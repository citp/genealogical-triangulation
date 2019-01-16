from collections import namedtuple, defaultdict
from itertools import chain, product, combinations
from os import listdir, makedirs
from os.path import join, exists, abspath, dirname
from shutil import rmtree
from warnings import warn
from random import shuffle, getstate, setstate, seed

from scipy.stats import gamma
import numpy as np
from progressbar import progressbar

from common_segments import common_segment_lengths
from data_logging import write_log
from population_genomes import generate_genomes
from population_statistics import all_ancestors_of
from gamma import fit_hurdle_gamma
from cm import centimorgan_data_from_directory
from shared_segment_detector import SharedSegmentDetector

ZERO_REPLACE = 1e-12

GammaParams = namedtuple("GammaParams", ["shape", "scale"])
HurdleGammaParams = namedtuple("HurdleGammaParams", ["shape", "scale", "zero_prob"])

DEFAULT_SMOOTHING = HurdleGammaParams(1.2088571040214136, 11686532.312642237, 0.9876864782229996)

class LengthClassifier:
    """
    Classifies based total length of shared segments
    """
    def __init__(self, distributions, labeled_nodes,
                 cryptic_distribution = None):
        self._distributions = distributions
        self._labeled_nodes = labeled_nodes
        self._by_unlabeled = None
        if cryptic_distribution is not None:
            self._cryptic_distribution = cryptic_distribution
        else:
            self._cryptic_distribution = DEFAULT_SMOOTHING

    @property
    def smoothing_parameters(self):
        smoothing = getattr(self, "_cryptic_distribution",
                            DEFAULT_SMOOTHING)
        if smoothing:
            return smoothing
        else:
            return DEFAULT_SMOOTHING

    @smoothing_parameters.setter
    def smoothing_parameters(self, params):
        self._cryptic_distribution = params

    @property
    def group_by_unlabeled(self):
        try:
            by_unlabeled = self._by_unlabeled
        except AttributeError:
            by_unlabeled = None
            
        if by_unlabeled is not None:
            return by_unlabeled

        by_unlabeled = defaultdict(set)
        for unlabeled_id, labeled_id in self._distributions.keys():
            by_unlabeled[unlabeled_id].add(labeled_id)

        by_unlabeled = dict(by_unlabeled)
        self._by_unlabeled = by_unlabeled
        return by_unlabeled
        

    def get_batch_smoothing_gamma(self, lengths):
        shape, scale, zero_prob = self.smoothing_parameters
        lengths = np.asarray(lengths, np.uint32)
        zero_i = (lengths == 0)
        nonzero_i = np.invert(zero_i)
        ret = np.empty_like(lengths, dtype = np.float64)
        gamma_probs = gamma.pdf(lengths[nonzero_i], a = shape, scale = scale)
        gamma_probs = gamma_probs * (1 - zero_prob)
        ret[nonzero_i] = gamma_probs
        ret[zero_i] = zero_prob
        ret[ret <= 0.0] = ZERO_REPLACE
        return ret

    def get_batch_pdf(self, lengths, query_nodes, labeled_nodes):
        assert len(lengths) == len(query_nodes) == len(labeled_nodes)
        assert len(lengths) > 0
        lengths = np.array(lengths, dtype = np.uint32)
        zero_i = (lengths == 0)
        nonzero_i = np.invert(zero_i)
        distributions = self._distributions
        params = (distributions[query_node, labeled_node]
                  for query_node, labeled_node
                  in zip(query_nodes, labeled_nodes))
        shape_scale_zero = list(zip(*params))
        shapes = np.array(shape_scale_zero[0], dtype = np.float64)
        scales = np.array(shape_scale_zero[1], dtype = np.float64)
        zero_prob = np.array(shape_scale_zero[2], dtype = np.float64)
        del shape_scale_zero
        ret = np.empty_like(lengths, dtype = np.float64)
        ret[zero_i] = zero_prob[zero_i]
        gamma_probs = gamma.pdf(lengths[nonzero_i],
                                a = shapes[nonzero_i],
                                scale = scales[nonzero_i])
        gamma_probs = gamma_probs * (1 - zero_prob[nonzero_i])
        ret[nonzero_i] = gamma_probs
        leq_zero = (ret <= 0.0)
        ret[leq_zero] = ZERO_REPLACE
        # ret[ret > 1.0] = 1.0
        return (ret, leq_zero)

    def __contains__(self, item):
        return item in self._distributions

def related_pairs(unlabeled_nodes, labeled_nodes, population, generations):
    """
    Given a population and labeled nodes, returns a list of pairs of nodes
    (unlabeled node, labeled node)
    where the labeled node and unlabeled node share at least 1 common ancestor
    going back generation generations from the latest generation.
    """
    if type(labeled_nodes) != set:
        labeled_nodes = set(labeled_nodes)
    ancestors = dict()
    ancestors = {node: all_ancestors_of(node)
                 for node
                 in chain(unlabeled_nodes, labeled_nodes)}

    return [(unlabeled, labeled) for unlabeled, labeled
            in product(unlabeled_nodes, labeled_nodes)
            if (unlabeled != labeled and
                not ancestors[unlabeled].isdisjoint(ancestors[labeled]))]


# TODO:
# At some point this should probably be turned into a "builder" class,
# too much state is getting passed along in this long parameter list.
def generate_classifier(population, labeled_nodes, genome_generator,
                        recombinators, directory, clobber = True,
                        iterations = 1000, generations_back_shared = 7,
                        min_segment_length = 0, non_paternity = 0.0):    
    if not exists(directory):
        makedirs(directory)
    elif clobber:
        rmtree(directory)
        makedirs(directory)
    num_generations = population.num_generations
    clear_index = max(num_generations - generations_back_shared, 0)
    to_clear = population.generations[clear_index].members
    for node in to_clear:
        node.suspected_mother = None
        node.suspected_mother_id = None
        node.suspected_father = None
        node.suspected_father_id = None
    if 0 < iterations:
        shared_to_directory(population, labeled_nodes, genome_generator,
                            recombinators, directory, clobber = clobber,
                            iterations = iterations,
                            min_segment_length = min_segment_length,
                            generations_back_shared = generations_back_shared,
                            non_paternity = non_paternity)


    
    print("Generating classifiers.")
    classifier = classifier_from_directory(directory, population.id_mapping)

    print("Calculating cryptic parameters")
    cryptic_params = cryptic_parameters(population.id_mapping,
                                        classifier._labeled_nodes,
                                        set(classifier._distributions.keys()))
    classifier._cryptic_distribution = cryptic_params

    return classifier

def cryptic_parameters(id_map, labeled_nodes, related_pairs):

    # TODO: Add argument for this
    recomb_dir = abspath(join(dirname(__file__),
                              "../data/recombination_rates/"))
    cm_data = centimorgan_data_from_directory(recomb_dir)
    ibd_detector = SharedSegmentDetector(0, 5, cm_data)

    # We try to only include the labeled nodes the analyst would have
    # access to This only works if we use the deterministic random
    # argument when we evaluate the classifier.
    labeled_copy = sorted(labeled_nodes)
    rand_state = getstate()
    seed(42)
    shuffle(labeled_copy)
    setstate(rand_state)
    
    labeled_node_pairs = set(combinations(labeled_copy[:1000], 2))
    related_pairs = set(related_pairs)
    cryptic_pairs = set(x for x in labeled_node_pairs if x not in related_pairs)
    lengths = []
    for node_a_id, node_b_id in cryptic_pairs:
        node_a = id_map[node_a_id]
        node_b = id_map[node_b_id]
        genome_a = node_a.genome
        genome_b = node_b.genome
        length = ibd_detector.shared_segment_length(genome_a, genome_b)
        lengths.append(length)
    np_lengths = np.array(lengths, dtype = np.uint64)
    params = fit_hurdle_gamma(np_lengths)
    assert all(x is not None for x in params)
    return params

def shared_to_directory(population, labeled_nodes, genome_generator,
                        recombinators, directory, min_segment_length = 0,
                        clobber = True, iterations = 1000,
                        generations_back_shared = 7,
                        non_paternity = 0.0):

    labeled_nodes = set(labeled_nodes)
    unlabeled_nodes = chain.from_iterable(generation.members
                                          for generation
                                          in population.generations[-3:])
    unlabeled_nodes = set(unlabeled_nodes) # - labeled_nodes
    print("Finding related pairs.")
    pairs = related_pairs(unlabeled_nodes, labeled_nodes, population,
                          generations_back_shared)
    print("{} related pairs.".format(len(pairs)))
    print("Opening file descriptors.")
    if clobber:
        mode = "w"
    else:
        mode = "a"
    fds = {node: open(join(directory, str(node._id)), mode)
           for node in labeled_nodes}
    print("Calculating shared lengths.")
    for i in range(iterations):
        print("iteration {}".format(i))
        print("Cleaning genomes.")
        population.clean_genomes()
        print("Generating genomes")
        generate_genomes(population, genome_generator, recombinators, 3,
                         true_genealogy = False)
        print("Calculating shared length")
        _calculate_shared_to_fds(pairs, fds, min_segment_length)
    for fd in fds.values():
        fd.close()
        
def _calculate_shared_to_fds(pairs, fds, min_segment_length):
    """
    Calculate the shared length between the pairs, and store the
    shared length in the given directory. Each labeled node has a file
    in the given directory. The files contain tab separated entries,
    where the first entry is the unlabeled node id, and the second
    entry is the amount of shared material.
    """
    shared_iter = ((unlabeled, labeled,
                    shared_segment_length_genomes(unlabeled.genome,
                                                  labeled.genome,
                                                  min_segment_length))
                   for unlabeled, labeled in pairs)
    for unlabeled, labeled, shared in shared_iter:
        fd = fds[labeled]
        fd.write("{}\t{}\n".format(unlabeled._id, shared))

def classifier_from_directory(directory, id_mapping):
    distributions = distributions_from_directory(directory, id_mapping)
    labeled_nodes = set(int(filename) for filename in listdir(directory))
    return LengthClassifier(distributions, labeled_nodes)

def distributions_from_directory(directory, id_mapping):
    """
    Calculate distributions from a directory created by
    calculate_shared_to_directory.
    """
    distributions = dict()
    for labeled_filename in progressbar(listdir(directory)):
        lengths = defaultdict(list)
        labeled = int(labeled_filename)
        with open(join(directory, labeled_filename), "r") as labeled_file:
            for line in labeled_file:
                # If the program crashed, the output can be left in an
                # inconsistent state.
                try:
                    unlabeled_id, shared_str = line.split("\t")
                except ValueError:
                    warn("Malformed line:\n{}".format(line), stacklevel = 0)
                    continue
                unlabeled = int(unlabeled_id)
                if unlabeled not in id_mapping:
                    error_string = "No such unlabeled node with id {}."
                    warn(error_string.format(unlabeled_id), stacklevel = 0)
                    continue
                try:
                    shared_float = float(shared_str)
                except ValueError:
                    error_string = "Error formatting value as float: {}."
                    warn(error_string.format(shared_str), stacklevel = 0)
                    continue
                lengths[unlabeled].append(shared_float)
        for unlabeled, lengths in lengths.items():
            shape, scale, zero_prob = fit_hurdle_gamma(np.array(lengths,
                                                                dtype = np.uint32))
            if shape is None:
                continue
            params = HurdleGammaParams(shape, scale, zero_prob)
            distributions[unlabeled, labeled] = params
    return distributions
    
def shared_segment_length_genomes(genome_a, genome_b, minimum_length):
    lengths = common_segment_lengths(genome_a, genome_b)
    seg_lengths = (x for x in lengths if x >= minimum_length)
    return sum(seg_lengths)

def shared_stats_genomes(genome_a, genome_b, minimum_length):
    lengths = [x for x in common_segment_lengths(genome_a, genome_b)
               if x >= minimum_length]
    num_segments = len(lengths)
    if num_segments == 0:
        return (0.0, 0)
    return (sum(lengths) / num_segments, num_segments)
    
def _shared_segment_length(node_a, node_b, minimum_length):
    return shared_segment_length_genomes(node_a.genome, node_b.genome,
                                          minimum_length)

from collections import namedtuple
from heapq import nlargest

# import pyximport; pyximport.install()
import numpy as np

from classify_relationship import LengthClassifier
from data_logging import write_log
from util import first_missing_ancestor, all_related

#TODO: Change this tuple, the last two values are not used.
ProbabilityData = namedtuple("ProbabilityData", ["start_i", "stop_i",
                                                 "cryptic_start_i",
                                                 "cryptic_stop_i"])

RawIdentified = namedtuple("RawIdentified", ["sibling_group",
                                             "ln_ratio",
                                             "identified_node"])
MINIMUM_LABELED_NODES = 5
INF = float("inf")
INF_REPLACE = 1.0

class BayesDeanonymize:
    def __init__(self, population, classifier = None, only_related = False,
                 search_generations_back = 7, cryptic_logging = False):
        self._population = population
        if classifier is None:
            self._length_classifier = LengthClassifier(population, 1000)
        else:
            self._length_classifier = classifier
        self._only_related = only_related
        if only_related:
            print("Only searching nodes related to labeled nodes within {} generations".format(search_generations_back))
            self._search_generations_back = search_generations_back
            self._compute_related()
        self._restrict_search_nodes = None
        self.cryptic_logging = cryptic_logging
        # self.__remove_erroneous_labeled()

    def __remove_erroneous_labeled(self):
        print("Removing erroneous labeled nodes")
        id_map = self._population.id_mapping
        labeled_nodes = [id_map[labeled_node_id] for labeled_node_id
                         in self._length_classifier._labeled_nodes]
        to_remove = set()
        for labeled_node in labeled_nodes:
            missing = first_missing_ancestor(labeled_node)
            if missing < 1:
                to_remove.add(labeled_node._id)
        new_labeled = set(self._length_classifier._labeled_nodes) - to_remove
        print("Removeing {} labeled nodes."
              " New size is {}.".format(len(to_remove), len(new_labeled)))
        self._length_classifier._labeled_nodes = list(new_labeled)

    def _to_search(self, shared_list, sex):
        labeled = set(self._length_classifier._labeled_nodes)
        genome_nodes = set(member for member in self._population.members
                           if (member.genome is not None and
                               member._id not in labeled and
                               member.sex == sex))
        if not self._only_related:
            return genome_nodes
        id_map = self._population.id_mapping
        potential_nodes = set()
        for labeled_node_id, shared in shared_list:
            if shared > 0:
                labeled_node = id_map[labeled_node_id]
                related = self._labeled_related[labeled_node]
                potential_nodes.update(related)
        if self._restrict_search_nodes is not None:
            return potential_nodes.intersection(self._restrict_search_nodes)
        else:
            return potential_nodes

    def _add_node_id_relatives(self, node_id, nodes):
        id_map = self._population.id_mapping
        labeled_node = id_map[node_id]
        related = all_related(labeled_node, True, self._search_generations_back)
        self._labeled_related[labeled_node] = related.intersection(nodes)

    def _compute_related(self):
        nodes = set(member for member in self._population.members
                    if member.genome is not None)
        self._labeled_related = dict()
        length_classifier = self._length_classifier
        for labeled_node_id in length_classifier._labeled_nodes:
            self._add_node_id_relatives(labeled_node_id, nodes)

    def add_labeled_node_id(self, node_id):
        self._length_classifier._labeled_nodes.append(node_id)
        if self._only_related:
            nodes = list(member for member in self._population.members
                         if member.genome is not None)
            self._add_node_id_relatives(node_id, nodes)

    def restrict_search(self, nodes):
        self._restrict_search_nodes = set(nodes)

    def identify(self, genome, actual_node, segment_detector):
        id_map = self._population.id_mapping
        length_classifier = self._length_classifier
        # TODO: Eliminated shared_list and use shared_dict everywhere
        shared_list = []
        for labeled_node_id in length_classifier._labeled_nodes:
            labeled_node = id_map[labeled_node_id]
            s = segment_detector.shared_segment_length(genome,
                                                       labeled_node.suspected_genome)
            shared_list.append((labeled_node_id, s))

        shared_dict = dict(shared_list)
        labeled_nodes = set(length_classifier._labeled_nodes)

        labeled_nodes_cryptic, all_lengths = list(zip(*shared_dict.items()))
        # We convert to python floats, as summing is faster.
        all_cryptic_possibilities = [float(x) for x
                                     in np.log(length_classifier.get_batch_smoothing_gamma(all_lengths))]
        # Maps labeled nodes to the log cryptic value of the IBD detected
        cryptic_lookup = dict(zip(labeled_nodes_cryptic,
                                  all_cryptic_possibilities))
        node_data = dict()
        batch_node_id = []
        batch_labeled_node_id = []
        batch_lengths = []
        # Keep for logging purposes
        # batch_cryptic_lengths = []
        node_cryptic_log_probs = dict()
        by_unlabeled = length_classifier.group_by_unlabeled
        nodes = self._to_search(shared_list, actual_node.sex)
        if len(nodes) == 0:
            # We have no idea which node it is
            return RawIdentified(set(), float("-inf"), None)
        
        empty = set()
        for node in nodes:
            node_start_i = len(batch_node_id)
            node_id = node._id
            cryptic_probability = 0
            node_cryptic_log_probs[node] = 0

            cryptic_nodes = labeled_nodes - by_unlabeled.get(node_id, empty)
            if len(cryptic_nodes) > 0:
                cryptic_probability = sum(cryptic_lookup[labeled_node_id]
                                          for labeled_node_id
                                          in cryptic_nodes)
                node_cryptic_log_probs[node] = cryptic_probability
        
            non_cryptic_nodes = list(labeled_nodes - cryptic_nodes)
            if len(non_cryptic_nodes) > 0:
                batch_node_id.extend([node_id] * len(non_cryptic_nodes))
                batch_labeled_node_id.extend(non_cryptic_nodes)
                batch_lengths.extend(shared_dict[labeled_node_id]
                                     for labeled_node_id in non_cryptic_nodes)
            
            node_stop_i = len(batch_node_id)
            #TODO: change this tuple, remove the last two values.
            node_data[node] = ProbabilityData(node_start_i, node_stop_i,
                                              -1, -1) 

        assert len(node_data) > 0
        if len(batch_lengths) > 0:
            pdf_vals = length_classifier.get_batch_pdf(batch_lengths,
                                                       batch_node_id,
                                                       batch_labeled_node_id)
            calc_prob, zero_replace = pdf_vals
        else:
            calc_prob = []

        node_probabilities = dict()
        for node, prob_data in node_data.items():
            start_i, stop_i, cryptic_start_i, cryptic_stop_i = prob_data
            node_calc = calc_prob[start_i:stop_i]
            log_prob = (np.sum(np.log(node_calc)) +
                        node_cryptic_log_probs[node])
            node_probabilities[node] = log_prob
        assert len(node_probabilities) > 0
        write_log("identify", {"node": actual_node._id,
                               "probs": {node._id: prob
                                         for node, prob
                                         in node_probabilities.items()}})
        # The value 8 is somewhat arbitrary. We are always able to
        # generate our confidence value with the top 8, as sibships
        # tend to be small. This number may need to be larger for
        # populations with large sibships.
        potential_nodes = nlargest(8, node_probabilities.items(),
                                   key = lambda x: x[1])
        top, top_log_prob = potential_nodes[0]
        sibling_group = get_suspected_sibling_group(top)
        for node, log_prob in potential_nodes[1:]:
            if node in sibling_group:
                continue
            next_node = node
            next_log_prob = log_prob
            break
        else:
            next_node, next_log_prob = potential_nodes[1]
                
        log_ratio  = top_log_prob - next_log_prob
        return RawIdentified(get_sibling_group(top), log_ratio, top)

def _get_logging_cryptic_lengths(shared_dict, cryptic_nodes, unique_lengths):
    lengths_iter = (shared_dict[labeled_node_id]
                    for labeled_node_id
                    in cryptic_nodes)
    temp_cryptic_lengths = np.fromiter(lengths_iter,
                                       dtype = np.uint64)
    lengths, counts = np.unique(temp_cryptic_lengths,
                                return_counts = True)
    store_counts = np.zeros(len(unique_lengths), dtype = np.uint32)
    i = np.searchsorted(unique_lengths, lengths)
    store_counts[i] = counts
    return store_counts

def get_sibling_group(node):
    """
    Returns the set containing node and all its full siblings
    """
    if node.mother is None or node.father is None:
        return set([node])
    return set(node.mother.children).intersection(node.father.children)

def get_suspected_sibling_group(node):
    if node.suspected_mother is None or node.suspected_father is None:
        return set([node])
    return set(node.suspected_mother.suspected_children).intersection(node.suspected_father.suspected_children)

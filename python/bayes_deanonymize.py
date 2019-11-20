from collections import namedtuple
from heapq import nlargest
from statistics import median

# import pyximport; pyximport.install()
import numpy as np

from classify_relationship import LengthClassifier
from calculate_probabilities import calculate_probabilities
from data_logging import write_log
from util import first_missing_ancestor, all_related

#TODO: Change this tuple, the last two values are not used.
ProbabilityData = namedtuple("ProbabilityData", ["node", "start_i", "stop_i", "cryptic_prob"])

RawIdentified = namedtuple("RawIdentified", ["sibling_group",
                                             "ln_ratio",
                                             "identified_node"])
MINIMUM_LABELED_NODES = 5
INF = float("inf")
INF_REPLACE = 1.0

class BayesDeanonymize:
    def __init__(self, population, classifier = None, only_related = False,
                 search_generations_back = 7, cryptic_logging = False,
                 probability_logging = True):
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
        self._exclude_search = set()
        self.exclude_anchors = set()
        self.probability_logging = probability_logging
        self._genome_nodes_cache = list(member for member in population.members
                                        if member.genome is not None)

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

    @property
    def exclude_from_search(self):
        return self._exclude_search

    @exclude_from_search.setter
    def exclude_from_search(self, to_exclude):
        self._exclude_search = set(to_exclude)


    def _to_search_basic(self, sex):
        genome_nodes = set(member for member in self._genome_nodes_cache
                           if member.sex == sex)
        id_map = self._population.id_mapping
        genome_nodes.difference_update(id_map[x] for x
                                       in self._length_classifier._labeled_nodes)
        return genome_nodes

    def _to_search_related(self, shared_list, sex):
        potential_nodes = set()
        id_map = self._population.id_mapping
        for labeled_node_id, shared in shared_list:
            if shared > 0:
                labeled_node = id_map[labeled_node_id]
                related = self._labeled_related[labeled_node]
                potential_nodes.update(node for node in related
                                       if node.sex == sex)
        if self._restrict_search_nodes is not None:
            return potential_nodes.intersection(self._restrict_search_nodes)
        else:
            return potential_nodes

    def _to_search(self, shared_list, sex):
        if self._only_related:
            ret = self._to_search_related(shared_list, sex)
        else:
            ret = self._to_search_basic(sex)

        if len(self._exclude_search) != 0:
            ret -= self._exclude_search
        return ret

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
        if node_id in self._length_classifier._labeled_nodes:
            return
        self._length_classifier._labeled_nodes.append(node_id)
        if self._only_related:
            nodes = list(member for member in self._population.members
                         if member.genome is not None)
            self._add_node_id_relatives(node_id, nodes)

    def remove_labeled_node_id(self, node_id):
        self._length_classifier._labeled_nodes.remove(node_id)

    def restrict_search(self, nodes):
        self._restrict_search_nodes = set(nodes)

    def identify(self, genome, actual_node, segment_detector):
        id_map = self._population.id_mapping
        length_classifier = self._length_classifier
        # TODO: Eliminated shared_list and use shared_dict everywhere
        shared_list = []
        anchors = set(length_classifier._labeled_nodes) - self.exclude_anchors
        sorted_labeled = sorted(anchors)
        np_sorted_labeled = np.array(sorted_labeled, dtype = np.uint32)
        sorted_shared = []
        for labeled_node_id in sorted_labeled:
            labeled_node = id_map[labeled_node_id]
            s = segment_detector.shared_segment_length(genome,
                                                       labeled_node.suspected_genome)
            shared_list.append((labeled_node_id, s))
            sorted_shared.append(s)

        write_log("positive ibd count", sum(0.0 < x for x in sorted_shared))
        #write_log("shared", sorted_shared)
        shared_dict = dict(shared_list)
        sorted_shared = np.array(sorted_shared, dtype = np.float64)

        labeled_nodes_cryptic, all_lengths = list(zip(*shared_dict.items()))
        np_cryptic = np.log(length_classifier.get_batch_smoothing_gamma(sorted_shared))

        node_data = []
        batch_shape = []
        batch_scale = []
        batch_zero_prob = []
        batch_lengths = []
        # Keep for logging purposes
        # batch_cryptic_lengths = []
        nodes = self._to_search(shared_list, actual_node.sex)
        if len(nodes) == 0:
            # We have no idea which node it is
            return RawIdentified(set(), float("-inf"), None)


        
        for node in nodes:
            node_start_i = len(batch_shape)
            node_id = node._id
            #node_cryptic_log_probs[node] = 0

            if node_id in length_classifier._distributions:
                labeled_ids, shape, scale, zero_prob = length_classifier._distributions[node_id]
            else:
                labeled_ids = np.array([], dtype = np.uint32)
                shape = scale = zero_prob = np.array([], dtype = np.float64)
            calc_data = calculate_probabilities(labeled_ids,
                                                shape, scale, zero_prob,
                                                sorted_shared,
                                                np_sorted_labeled,
                                                np_cryptic,
                                                node_id)
            cur_lengths, cur_shapes, cur_scales, cur_zero_prob, cur_cryptic = calc_data
            batch_lengths.extend(cur_lengths)
            batch_shape.extend(cur_shapes)
            batch_scale.extend(cur_scales)
            batch_zero_prob.extend(cur_zero_prob)
            
            node_stop_i = len(batch_shape)
            node_data.append(ProbabilityData(node, node_start_i, node_stop_i,
                                             cur_cryptic))


        assert len(node_data) > 0
        if len(batch_lengths) > 0:
            pdf_vals = length_classifier.batch_pdf_distributions(batch_lengths,
                                                                 batch_shape,
                                                                 batch_scale,
                                                                 batch_zero_prob)
            calc_prob, zero_replace = pdf_vals
        else:
            calc_prob = []

        log_calc_prob_cum = np.cumsum(np.log(calc_prob))
        del calc_prob
        log_calc_prob_cum = np.concatenate(([0.0], log_calc_prob_cum))
        node_probabilities = dict()
        for node, start_i, stop_i, cryptic_prob in node_data:
            log_prob = (log_calc_prob_cum[stop_i] - log_calc_prob_cum[start_i]) + cryptic_prob
            node_probabilities[node] = log_prob
        assert len(node_probabilities) > 0
        if self.probability_logging:
            write_log("identify", {"node": actual_node._id,
                                   "probs": {node._id: prob
                                             for node, prob
                                             in node_probabilities.items()}})

        if len(node_probabilities) == 0:
            return  RawIdentified(set(), -INF, None)
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
            if len(potential_nodes) > 1:
                next_node, next_log_prob = potential_nodes[1]

        if len(potential_nodes) > 1:
            log_ratio  = top_log_prob - next_log_prob
        else:
            log_ratio = -INF
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

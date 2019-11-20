from collections import Counter, defaultdict, namedtuple
from datetime import datetime
from pickle import dump
from sys import stdout
from math import sqrt
from random import shuffle, sample
from time import perf_counter

from scipy.stats import norm

from data_logging import write_log
from bayes_deanonymize import BayesDeanonymize
from shared_segment_detector import SharedSegmentDetector
from util import ancestor_roots, all_related


IdentifyResult = namedtuple("IdentifyResult", ["target_node",
                                               "sibling_group",
                                               "identified_node",
                                               "ln_ratio",
                                               "correct",
                                               "run_number"])
class Evaluation:
    def __init__(self, population, classifier, labeled_nodes = None,
                 ibd_detector = None, search_related = False,
                 smoothing_parameters = None, cryptic_logging = False,
                 out_of_genealogy = 0, randomize_labeled = False):
        self._population = population
        self._classifier = classifier
        if labeled_nodes is not None:
            self._original_labeled = classifier._labeled_nodes.copy()
            self.labeled_nodes = labeled_nodes
        if randomize_labeled:
            self._randomize_labeled = len(self.labeled_nodes)
        else:
            self._randomize_labeled = 0
        if search_related:
            print("Calculating related nodes")
            self._bayes = BayesDeanonymize(population, classifier,
                                           True, search_related,
                                           cryptic_logging = cryptic_logging)
        else:
            self._bayes = BayesDeanonymize(population, classifier, False,
                                           cryptic_logging = cryptic_logging)
        if ibd_detector is None:
            ibd_detector = SharedSegmentDetector()
        if smoothing_parameters:
            classifier.smoothing_parameters = smoothing_parameters
        assert isinstance(ibd_detector, SharedSegmentDetector)
        self._ibd_detector = ibd_detector
        self._run_number = 0
        self._out_of_genealogy = out_of_genealogy
        self.reset_metrics()

    @property
    def probability_logging(self):
        return self._bayes.probability_logging

    @probability_logging.setter
    def probability_logging(self, value):
        self._bayes.probability_logging = bool(value)

    @property
    def accuracy(self):
        total = self.correct + self.incorrect
        return self.correct / total

    def add_labeled_node_id(self, node):
        self._bayes.add_labeled_node_id(node)

    @property
    def labeled_nodes(self):
        return self._classifier._labeled_nodes.copy()

    @labeled_nodes.setter
    def labeled_nodes(self, labeled_nodes):
        # TODO: Handle the case when we are only searching x
        # generations back here.
        self._classifier._labeled_nodes = list(labeled_nodes)

    def restrict_search(self, nodes):
        assert len(nodes) > 0
        self._bayes.restrict_search(nodes)

    def print_metrics(self):
        total = self.correct + self.incorrect
        print("{} correct, {} incorrect, {} total.".format(self.correct,
                                                           self.incorrect,
                                                           total))
        stdout.flush()

        write_log("correct", self.correct)
        write_log("incorrect", self.incorrect)
        write_log("total", total)
        percent_accurate = self.accuracy
        plus_minus = confidence_interval(self.correct, total)
        print("{}±{:0.3} percent accurate.".format(percent_accurate, plus_minus))
        for generation, counter in self.generation_error.items():
            gen_correct = counter["correct"]
            gen_incorrect = counter["incorrect"]
            total = gen_correct + gen_incorrect
            format_string = "For generation {}: {} accuracy, {} total."
            print(format_string.format(generation, gen_correct / total, total))

    def _evaluate_node(self, node):
        raw_identified = self._bayes.identify(node.genome, node,
                                              self._ibd_detector)
        sibling_group, ln_ratio, identified_node = raw_identified
        node_generation = self._population.node_to_generation[node]
        print("Confidence score: {}".format(ln_ratio))
        if node in sibling_group:
            self.generation_error[node_generation]["correct"] += 1
            self.correct += 1
            print("correct")
        else:
            self.generation_error[node_generation]["incorrect"] += 1
            print("incorrect")
            self.incorrect += 1
            
        write_log("evaluate", {"target node": node._id,
                               "log ratio": ln_ratio,
                               "identified": set(x._id for x in sibling_group),
                               "run_number": self._run_number})
        stdout.flush()
        return IdentifyResult(node,
                              sibling_group,
                              identified_node,
                              ln_ratio,
                              node in sibling_group,
                              self._run_number)

    def _remove_from_classifier(self, nodes):
        node_ids = set(node._id for node in nodes)
        self._bayes.exclude_anchors = node_ids
        self._bayes.exclude_from_search = nodes

    def _to_remove(self, node):
        out_of_genealogy = self._out_of_genealogy
        ancestors = ancestor_roots(node, False, out_of_genealogy)
        ancestors.update(ancestor_roots(node, True, out_of_genealogy))
        generation_map = self._population.node_to_generation
        node_generation = generation_map[node]
        save_generation = node_generation - out_of_genealogy
        to_save = set(a for a in ancestors
                      if generation_map[a] == save_generation)
        related = all_related(node, False, out_of_genealogy)
        related.update(all_related(node, True, out_of_genealogy))
        return related - to_save


    def _remove_node_from_genealogy(self, node):
        to_remove = self._to_remove(node)
        assert node in to_remove
        self._remove_from_classifier(to_remove)

    def _reset_out_of_genealogy(self):
        self._bayes.exclude_from_search = set()
        self._bayes.exclude_anchors = set()

    def reset_metrics(self):
        # Maps generation -> counter with keys "correct" and "incorrect"
        self.generation_error = defaultdict(Counter)
        self.identify_results = []
        self.correct = 0
        self.incorrect = 0

    def run_expansion_round(self, identify_candidates, confidence_ratio,
                            expansion_data = None, expansion_filename = None,
                            revisit = True):
        print("Running expansion round.")
        if revisit:
            to_evaluate = list(identify_candidates)
        else:
            id_map = self._population.id_mapping
            labeled_genomes = set(id_map[x].suspected_genome for x
                              in self.labeled_nodes)
            to_evaluate = [x for x in identify_candidates
                           if x.suspected_genome not in labeled_genomes]
        shuffle(to_evaluate)
        added = []
        correct_add_count = 0
        new_added = 0
        self._bayes.probability_logging = False
        write_log("expansion_confidence_ratio", confidence_ratio)
        for i, node in enumerate(to_evaluate):
            self.run_evaluation([node], expansion = i)
            result = self.identify_results[-1]
            if result.ln_ratio > confidence_ratio:
                print("Adding node.")
                
                identified_node = result.identified_node
                added.append(result)
                if not expansion_data.identified_before(result.target_node):
                    new_added += 1
                if result.correct:
                    correct_add_count += 1
                prev_added = expansion_data.add_node(result)
                if prev_added != identified_node._id:
                    self._bayes.add_labeled_node_id(identified_node._id)
                    if prev_added is not None:
                        self._bayes.remove_labeled_node_id(prev_added._id)
                
            if 0 < i and i % 20 == 0:
                self.print_metrics()
                print("Nodes added this round: {}, New nodes added: {}".format(len(added), new_added))
                print("Correct nodes added: {}".format(correct_add_count))
            if expansion_data and expansion_filename and i % 500 == 0:
                remaining = set(node._id for node in to_evaluate[i:])
                expansion_data.remaining = remaining
                with open(expansion_filename, "wb") as expansion_file:
                    dump(expansion_data, expansion_file)
                write_log("expansion_data_written", {"current_node": node._id,
                                                     "complete": False})
        expansion_data.remaining = None
        write_log("expansion_round", {"added": len(added),
                                      "evaluated": len(to_evaluate),
                                      "correct_added": correct_add_count,
                                      "accuracy": self.accuracy})
        self.print_metrics()
        print("Added {} nodes this round.".format(len(added)))
        if len(added) == 0:
            input("No nodes added. Press enter to continue")
        return added

    def run_evaluation(self, unlabeled, expansion = None):
        # generation_map = population.node_to_generation
        # write_log("labeled_nodes", [node._id for node in labeled_nodes])
        # write_log("target_nodes", [node._id for node in unlabeled])
        if expansion is None:
            print("Attempting to identify {} random nodes.".format(len(unlabeled)),
                  flush = True)
        write_log("start time", datetime.now())
        for i, node in enumerate(unlabeled):
            if expansion:
                i = expansion
            print("Iteration: {}, actual node ID: {}".format(i + 1, node._id))
            start_time = perf_counter()
            if not expansion and self._randomize_labeled:
                potential = self._original_labeled.copy()
                if node._id in potential:
                    potential.remove(node._id)
                assert len(potential) > self._randomize_labeled
                new_anchors = sample(potential, self._randomize_labeled)
                self.labeled_nodes = new_anchors
            if self._out_of_genealogy > 0:
                self._reset_out_of_genealogy()
                self._remove_node_from_genealogy(node)
            self.identify_results.append(self._evaluate_node(node))
            end_time = perf_counter()
            print("It took {} seconds".format(end_time - start_time))

        write_log("end time", datetime.now())
        self._run_number += 1
        return self.identify_results

# from http://www.statsmodels.org/devel/_modules/statsmodels/stats/proportion.html#proportion_confint
def confidence_interval(correct, total, sig = 0.05):
    frac_correct = correct / total
    std = sqrt(frac_correct * (1 - frac_correct) / total)
    dist = norm.isf(sig / 2) * std
    return dist

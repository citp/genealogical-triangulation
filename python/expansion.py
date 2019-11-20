from evaluation import IdentifyResult

def flat_copy_identify_result(identify_result):
    target_node = identify_result.target_node._id
    assert isinstance(target_node, int)
    identified_node = identify_result.identified_node._id
    assert isinstance(identified_node, int)
    sibling_group = frozenset(node._id for node
                              in identify_result.sibling_group)
    return IdentifyResult(target_node, sibling_group, identified_node,
                          identify_result.ln_ratio, identify_result.correct,
                          identify_result.run_number)

class ExpansionData:
    def __init__(self, start_labeled):
        self.start_labeled = start_labeled
        self.added = []
        self._rounds = 0
        self._remaining = None
        self._original_pool = None
        self._target_identified = dict()
        self._identified_target = dict()

    def add_round(self):
        self._rounds += 1


    def _add_identification(self, correct, identified, target):
        old_identified = None
        if identified in self._identified_target:
            old_target = self._identified_target[identified]
            del self._target_identified[old_target]
        if target in self._target_identified:
            old_identified = self._target_identified[target]
            old_identified.suspected_genome = None
            del self._identified_target[old_identified]
        if correct:
            identified.suspected_genome = None
        else:
            identified.suspected_genome = target.genome
        self._target_identified[target] = identified
        self._identified_target[identified] = target
        return old_identified

    def identified_before(self, node):
        return node in self._target_identified

    def add_node(self, result):
        """
        Adds the new identification to the expansion data.
        If result.target_node has been identified before, returns the
        node that it was identified as, otherwise returns None
        """
        prev_identified = self._add_identification(result.correct,
                                                   result.identified_node,
                                                   result.target_node)
        self.added.append(flat_copy_identify_result(result))
        return prev_identified

    def update_remaining(self, remaining):
        self._remaining = remaining

    def adjust_genomes(self, population):
        self._deduplicate()
        id_mapping = population.id_mapping
        for result in self.added:
            identified = id_mapping[result.identified_node]
            target = id_mapping[result.target_node]
            self._add_identification(result.correct, identified, target)

    # TODO Consider callling this in __setstate__ so we don't have to
    # worry about when we call adjust_genomes.
    def _deduplicate(self):
        to_keep = []
        identified_targets = set()
        for identify in reversed(self.added):
            if identify.target_node in identified_targets:
                continue
            identified_targets.add(identify.target_node)
            to_keep.append(identify)
        to_keep.reverse()
        self.added = to_keep

    def __getstate__(self):
        state = self.__dict__.copy()
        state["_target_identified"] = dict()
        state["_identified_target"] = dict()
        return state

    @property
    def original_pool(self):
        return self._original_pool

    @original_pool.setter
    def original_pool(self, val):
        self._original_pool = val

    @property
    def remaining(self):
        return self._remaining

    @remaining.setter
    def remaining(self, remaining):
        self._remaining = remaining

    @property
    def labeled_nodes(self):
        ret = list(self.start_labeled)
        ret.extend(result.identified_node for result in self.added)
        return ret

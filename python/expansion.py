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
        self.added = set()
        self._rounds = 0
        self._remaining = None

    def add_round(self, to_add):
        self._rounds += 1
        self.extend_round(to_add)

    def extend_round(self, to_add, remaining = None):
        self.added.update(flat_copy_identify_result(res) for res in to_add)
        self._remaining = remaining

    def adjust_genomes(self, population):
        id_mapping = population.id_mapping
        for result in self.added:
            if result.correct:
                continue
            identified = id_mapping[result.identified_node]
            target = id_mapping[result.target_node]
            identified.suspected_genome = target.genome

    @property
    def remaining(self):
        return self._remaining

    @property
    def labeled_nodes(self):
        ret = list(self.start_labeled)
        ret.extend(result.identified_node for result in self.added)
        return ret

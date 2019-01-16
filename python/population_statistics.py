from itertools import combinations, chain
from random import sample
from collections import Counter

import numpy as np

from util import descendants_of, get_sample_of_cousins

def proportion_within_distance(population, distance, percent_labeled):
    members = population.generations[-1].members
    num_labeled = int(percent_labeled * len(members))
    assert num_labeled > 0
    labeled = sample(members, num_labeled)
    ancestors = list()
    to_search = [(node, distance) for node in labeled]
    while len(to_search) > 0:
        current_node, current_distance = to_search.pop()
        assert current_node.mother is not None
        assert current_node.father is not None
        if current_distance == 1:
            ancestors.append(current_node.mother)
            ancestors.append(current_node.father)
            continue
        to_search.append((current_node.mother, current_distance - 1))
        to_search.append((current_node.father, current_distance - 1))
    assert len(to_search) is 0
    descendant_set = set()
    to_search.extend(ancestors)
    while len(to_search) > 0:
        current_node = to_search.pop()
        if current_node in descendant_set:
            continue
        descendant_set.add(current_node)
        to_search.extend(current_node.children)
    return len(descendant_set.intersection(members)) / len(members)

def ancestors_of(node, distance, suspected = True):
    assert distance > 0
    to_visit = [(node, 0)]
    ancestors = set()
    while len(to_visit) > 0:
        current_node, current_distance = to_visit.pop()
        assert current_node is not None
        if current_distance == distance:
            ancestors.add(current_node)
            continue
        if suspected:
            mother = current_node.suspected_mother
            father = current_node.suspected_father
        else:
            mother = current_node.mother
            father = current_node.father
        if mother is not None:
            to_visit.append((mother, current_distance + 1))
        if father is not None:
            to_visit.append((father, current_distance + 1))
    return ancestors

def all_ancestors_of(node, suspected = True):
    """
    Return the set set of all ancestors, including the given node.
    """
    to_visit = [node]
    ancestors = set([node])
    while len(to_visit) > 0:
        current_node = to_visit.pop()
        if suspected:
            mother = current_node.suspected_mother
            father = current_node.suspected_father
        else:
            mother = current_node.mother
            father = current_node.father
        if mother is not None:
            ancestors.add(mother)
            to_visit.append(mother)
        if father is not None:
            ancestors.add(father)
            to_visit.append(father)
    return ancestors

def cdf_num_labeled_within_distance(population, distance, percent_labeled):
    # TODO: This function only works if there is strict monogamy
    members = population.generations[-1].members
    num_labeled = int(percent_labeled * len(members))
    assert num_labeled > 0
    labeled = sample(members, num_labeled)
    ancestors = set(chain.from_iterable(ancestors_of(labeled_node, distance)
                                        for labeled_node in labeled))

    descendants = chain.from_iterable(descendants_of(ancestor)
                                      for ancestor in ancestors)
    counter = Counter(descendants)
    
    # Remove everyone not in the last generation
    extras = set(counter.keys()) - set(members)
    for extra in extras:
        del counter[extra]

    for member in members:
        if member not in counter:
            counter[member] = 0
        else:
            # We double count ancestors,
            # eg if by brother is in the set, then I show up twice in the
            # count because I am added from both of my parent's sides. Fix
            # this when we move away from strict monogamy assumptions.
            counter[member] //= 2

    # key is the number of people labeled within distance, and value
    # is the number of people who have this many labeled people within
    # that distance.
    # eg if value_counts[3] is 5, then 5 people have relations of the
    # given distance with with 3 labeled individuals.
    # value_counts = Counter(counter.values())
    # assert len(members) == sum(value_counts.values()), \
    #     "Expected len(members) to equal sum(value_counts.values()), got {} and {} instead.".format(len(members), sum(value_co\unts.values()))
    #  total_count = len(members)
    # x = sorted(value_counts.keys())
    # y = np.cumsum([value_counts[key] for key in x]) / total_count
    # return (x, y)
    return np.array(list(counter.values()))

def shared_segment_distribution(population, distance):
    assert distance > 0
    cousin_pairs = get_sample_of_cousins(population, distance)
    lengths = []
    for p_1, p_2 in cousin_pairs:
            by_autosome = common_segment_lengths(p_1.genome, p_2.genome)
            lengths.extend(chain.from_iterable(by_autosome.values()))
    return lengths

from itertools import combinations, product, chain
from random import sample
from collections import deque

inf = float("inf")

# TODO: Deduplicate with version in population_statistics. Very confusing.
# Make sure to get suspected parameter correct in transition.
def all_ancestors(node, suspected = False):
    ancestors = set([node])
    current_generation = set([node])
    while len(current_generation) > 0:
        if suspected:
            mothers = set(node.suspected_mother for node in current_generation
                          if node.suspected_mother is not None)
            fathers = set(node.suspected_father for node in current_generation
                          if node.suspected_father is not None)
        else:
            mothers = set(node.mother for node in current_generation
                          if node.mother is not None)
            fathers = set(node.father for node in current_generation
                          if node.father is not None)
        current_generation = set(chain(mothers, fathers))
        ancestors.update(current_generation)
    return ancestors

def get_path_descendant(node_a, node_b, suspected = False):
    """
    Returns a list of nodes giving a path from node_a to node_b
    node_b must be an ancestor of node_a.
    """
    previous = dict()
    visited = set()
    queue = deque([node_a])
    while len(queue) > 0:
        current = queue.pop()
        if current in visited:
            continue
        if current == node_b:
            break
        if suspected:
            father = current.suspected_father
            mother = current.suspected_mother
        else:
            father = current.father
            mother = current.mother
        if father is not None:
            previous[father] = current
            queue.appendleft(father)
        if mother is not None:
            previous[mother] = current
            queue.appendleft(mother)
        visited.add(current)
    path = [node_b]
    while path[-1] in previous:
        path.append(previous[path[-1]])
    path.reverse()
    return path

def first_missing_ancestor(node, suspected = True):
    queue = deque([(node, 0)])
    while len(queue) > 0:
        current, distance = queue.pop()
        if suspected:
            father = current.suspected_father
            mother = current.suspected_mother
        else:
            father = current.father
            mother = current.mother
        if mother is None or father is None:
            return distance
        queue.appendleft((current.mother, distance + 1))
        queue.appendleft((current.father, distance + 1))
    return inf

def closest_error(node):
    visited = set()
    to_explore = deque([(0, node)])
    while len(to_explore) > 0:
        distance, node = to_explore.popleft()
        children = set(node.children)
        if (node.mother != node.suspected_mother
            or node.father != node.suspected_father
            or children != set(node.suspected_children)):
            return (distance, node)
        visited.add(node)
        if node.mother not in visited:
            to_explore.append((distance + 1, node.mother))
        if node.father not in visited:
            to_explore.append((distance + 1, node.father))
        to_explore.extend((distance + 1, n) for n in children - visited)
    return (None, None)

def closest_error_descendants(node):
    descendants = deque([(0, node)])
    while len(descendants) > 0:
        distance, node = descendants.popleft()
        children = set(node.children)
        if children != set(node.suspected_children):
            return (distance, node)
        descendants.extend((distance + 1, n) for n in children)
    return (None, None)

def closest_error_ancestors(node):
    ancestors = deque([(0, node)])
    while len(ancestors) > 0:
        distance, node = ancestors.popleft()
        if (node.mother != node.suspected_mother
            or node.father != node.suspected_father):
            return (distance, node)
        ancestors.append((distance + 1, node.mother))
        ancestors.append((distance + 1, node.father))
    return (None, None)

def closest_error_descendants_ancestors(node):
    """
    Return the closest error edge between the node and one of its
    ancestors or descendants.
    """
    ancestor_error = closest_error_ancestors(node)
    descendant_error = closest_error_descendants(node)

    if ancestor_error is None and descendant_error is None:
        return (None, None)
    if ancestor_error is None:
        return descendant_error
    if descendant_error is None:
        return ancestor_error

    if ancestor_error[0] < descendant_error[0]:
        return ancestor_error
    else:
        return descendant_error

def error_on_path(node_a, node_b, suspected = False):
    """
    Determine if there is an error on the path from node_a to node_b,
    where node_a is a descendant of node b.  Returns a list of
    distances along the path where errors occur.

    The path can be made by following the "suspected_mother" and
    father links or by following the ground truth data, and there is
    an error if these two paths are different.

    If suspected is "True" then we assume this path was built using
    the suspected parentage, otherwise the path was constructed using
    the ground truth parentage.
    """
    path = get_path_descendant(node_a, node_b, suspected)
    error_nodes = []
    for current_node, next_node in zip(path, path[1:]):
        if suspected:
            father = current_node.suspected_father
            mother = current_node.suspected_mother
        else:
            father = current_node.father
            mother = current_node.mother
        if father == next_node:
            if current_node.father != current_node.suspected_father:
                error_nodes.append(current_node)
        elif mother == next_node:
            if current_node.mother != current_node.suspected_mother:
                error_nodes.append(current_node)
        else:
            raise Exception("Found non child link in path")
    return [path.index(node) for node  in error_nodes]

def error_between_nodes(node_a, node_b, generation_map, suspected = True):
    """
    Determine the amount of error between node_a and node_b through
    their most recent common ancestor.

    Returns a pair, wit the first element being a list of distances
    where error occurs from node_a to the recent common ancestor, and
    the second being the same from node_b.
    """
    rca, dist = recent_common_ancestor(node_a, node_b, generation_map,
                                       suspected)
    if rca is None:
        Exception("Given nodes are not related. This error may occur if suspected parameter is incorrect")
    node_a_error = error_on_path(node_a, rca, suspected)
    node_b_error = error_on_path(node_b, rca, suspected)
    return (node_a_error, node_b_error)

def recent_common_ancestor(node_a, node_b, generation_map, suspected = False):
    assert node_a is not None
    assert node_b is not None
    if node_a == node_b:
        return (node_a, 0, 0)
    if not suspected: #TODO: Do something less verbose here
        if node_a.mother is not None and node_a.mother == node_b.mother:
            return (node_a.mother, 1, 2)
        if node_a.father is not None and node_a.father == node_b.father:
            return (node_a.father, 1, 2)
    a_ancestors = all_ancestors(node_a, suspected)
    b_ancestors = all_ancestors(node_b, suspected)
    common_ancestors = a_ancestors.intersection(b_ancestors)
    if len(common_ancestors) == 0:
        return (None, None, None)
    recent_generation = -1
    ancestor = None
    for common_ancestor in common_ancestors:
        ancestor_generation = generation_map[common_ancestor]
        if ancestor_generation > recent_generation:
            recent_generation = ancestor_generation
            ancestor = common_ancestor
    younger_pair = max(generation_map[node_a], generation_map[node_b])
    node_a_distance = abs(generation_map[node_a] - recent_generation)
    node_b_distance = abs(generation_map[node_b] - recent_generation)
    total_distance = node_a_distance + node_b_distance
    return (ancestor, abs(younger_pair - recent_generation), total_distance)

def get_sample_of_cousins(population, distance, percent_ancestors = 0.1,
                          percent_descendants = 0.1):
    """
    return a sample of pairs of individuals whos most recent common
    ancestor is exactly generations back.
    """
    assert 0 < distance < len(population.generations)
    assert 0 < percent_descendants <= 1
    assert 0 < percent_ancestors <= 1
    common_ancestors = population.generations[-(distance + 1)].members
    last_generation = set(population.generations[-1].members)
    ancestors_sample = sample(common_ancestors,
                              int(len(common_ancestors) * percent_ancestors))
    pairs = []
    for ancestor in ancestors_sample:
        temp_pairs = descendants_with_common_ancestor(ancestor, last_generation)
        temp_pairs = list(temp_pairs)
        pairs.extend(sample(temp_pairs,
                            int(len(temp_pairs) * percent_descendants)))
    return pairs

def ancestor_roots(node, suspected = False, generations_back = float("inf")):
    """"
    Returns all ancestor nodes that have no known ancestors themselves.
    """
    current_generation = [node]
    ancestors = set()
    i = 0
    while i < generations_back and len(current_generation) > 0:
        if suspected:
            mothers = (node.suspected_mother for node in current_generation
                       if node.suspected_mother is not None)
            fathers = (node.suspected_father for node in current_generation
                       if node.suspected_father is not None)
            new_ancestors = (node for node in current_generation
                             if node.suspected_father is None and
                             node.suspected_mother is None)
        else:
            mothers = (node.mother for node in current_generation
                       if node.mother is not None)
            fathers = (node.father for node in current_generation
                       if node.father is not None)
            new_ancestors = (node for node in current_generation
                             if node.father is None and
                             node.mother is None)
        ancestors.update(new_ancestors)
        current_generation = set(chain(mothers, fathers))
        i += 1
    ancestors.update(current_generation)
    return ancestors

def all_related(node, suspected = False, generations_back = float("inf")):
    """
    Returns all nodes related to the given node.
    """
    ancestors = ancestor_roots(node, suspected,
                               generations_back = generations_back)
    if len(ancestors) == 0:
        return set()
    descendant_sets = (descendants_of(ancestor, suspected)
                       for ancestor in ancestors)
    return set(chain.from_iterable(descendant_sets))

def descendants_of(node, suspected = False):
    descendants = set()
    to_visit = [node]
    while len(to_visit) > 0:
        ancestor = to_visit.pop()
        descendants.add(ancestor)
        if suspected:
            to_visit.extend(ancestor.suspected_children)
        else:
            to_visit.extend(ancestor.children)
    return descendants

def descendants_with_common_ancestor(ancestor, generation_members):
    """
    Returns pairs of individuals descendent from ancestor in the given
    generation who have ancestor as their most recent ancestor.
    """
    # Find the descendents of the children, remove the pairwise
    # intersection, and return pairs from different sets.
    ancestor_children = ancestor.children
    if len(ancestor_children) < 2:
        return []
    if generation_members.issuperset(ancestor_children):
        # Depth is only 1 generation, so return all combinations of children.
        return combinations(ancestor_children, 2)
    descendant_sets = [descendants_of(child).intersection(generation_members)
                       for child in ancestor_children]    
    pair_iterables = []
    for descendants_a , descendants_b in combinations(descendant_sets, 2):
        intersection = descendants_a.intersection(descendants_b)
        if len(intersection) > 0:
            # Remove individuals who have a more recent common ancestor
            descendants_a = descendants_a - intersection
            descendants_b = descendants_b - intersection
        pair_iterables.append(product(descendants_a, descendants_b))
    return chain.from_iterable(pair_iterables)


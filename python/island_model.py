class IslandNode():
    def __init__(self, switch_probability, members = None, is_leaf = False):
        """
        Create a new island node with constituent members.  If is_leaf
        is true, then members will be individuals (Node objects from
        node.py that represent people).
        """
        self._parent = None
        self._is_leaf = is_leaf
        self._switch_probability = switch_probability
        if members is None:
            members = set()
        # Don't use a generic "members" container, to potentially
        # catch mistakes where user think this node is a leaf when it
        # is not, or vice versa.
        if is_leaf:
            self._individuals = set(members)
        else:
            self._islands = set(members)
            for node in members:
                node._parent = self

    def _add_individual(self, individual):
        self._individuals.add(individual)

    def _remove_individual(self, individual):
        self._individuals.remove(individual)

    def _add_island(self, island):
        self._islands.add(island)
        island._parent = self

    def __in__(self, individual):
        return individual in self._individuals

    @property
    def is_leaf(self):
        return self._is_leaf

    @property
    def islands(self):
        return self._islands

    @property
    def individuals(self):
        return self._individuals

    @property
    def parent(self):
        return self._parent

    @property
    def switch_probability(self):
        return self._switch_probability


class IslandTree():
    """
    Class to encapsulate switching probabilities in Hierarchical
    island tree model.
    """
    def __init__(self, root_node):
        self._root = root_node
        self._leaves = []
        self._individual_island = {}
        nodes = []
        nodes.append(root_node)
        while len(nodes) > 0:
            current_node = nodes.pop()
            if current_node.is_leaf:
                self._leaves.append(current_node)
                self._individual_island.update((individual, current_node)
                                               for individual
                                               in current_node.individuals)
            else:
                nodes.extend(current_node.islands)

    def get_island(self, individual):
        return self._individual_island[individual]

    def add_individual(self, island, individual):
        island._add_individual(individual)
        self._individual_island[individual] = island

    def move_individual(self, individual, destination_island):
        """
        Remove an individual from their current island and move them
        to destination_island.
        """
        source_island = self._individual_island[individual]
        if source_island == destination_island:
            return
        source_island._remove_individual(individual)
        destination_island._add_individual(individual)
        self._individual_island[individual] = destination_island
        
                
    @property
    def root(self):
        return self._root

    @property
    def leaves(self):
        return self._leaves

    @property
    def individuals(self):
        return list(self._individual_island)

def tree_from_string(tree):
    if isinstance(tree, str):
        tree = tree.split("\n")
    nodes = []
    indentation = [-1]
    previous_node = None
    for line in tree:
        if line.lstrip().startswith("#"):
            # Ignore comments in tree file.
            continue
        line_indentation = len(line) - len(line.lstrip())
        current_node = (float(line.strip()), [])
        if line_indentation > indentation[-1]:
            if previous_node is None:
                previous_node = current_node
                continue
            nodes.append(previous_node)
            indentation.append(line_indentation)
        elif line_indentation < indentation[-1]:
            nodes.pop()
            indentation.pop()
        nodes[-1][1].append(current_node)
        previous_node = current_node
        
    root_node = _tree_from_string_helper(nodes[0])
    return IslandTree(root_node)

def _tree_from_string_helper(node):
    if len(node[1]) is 0:
        return IslandNode(node[0], is_leaf = True)
    else:
        children = [_tree_from_string_helper(child) for child in node[1]]
        return IslandNode(node[0], children)

def tree_from_file(f):
    with open(f, "r") as tree_file:
        return tree_from_string(tree_file.readlines())

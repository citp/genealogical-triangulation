class IslandNode():
    def __init__(self, switch_probability, members = None):
        """
        Create a new island node with constituent members.  If is_leaf
        is true, then members will be individuals (Node objects from
        node.py that represent people).
        """
        self._switch_probability = switch_probability
        if members is None:
            members = set()
        self._individuals = set(members)

    def _add_individual(self, individual):
        self._individuals.add(individual)

    def _remove_individual(self, individual):
        self._individuals.remove(individual)

    def __in__(self, individual):
        return individual in self._individuals

    @property
    def individuals(self):
        return self._individuals

    @property
    def switch_probability(self):
        return self._switch_probability


class IslandModel():
    """
    Class to encapsulate switching probabilities in Hierarchical
    island tree model.
    """
    def __init__(self, islands):
        assert islands is not None
        self._islands = islands
        self._individual_island = dict()
        for island in islands:
            for node in island.individuals:
                self._individual_island[node] = island

    def get_island(self, individual):
        return self._individual_island[individual]

    # TODO: Move this method to the island object
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
    def individuals(self):
        return list(self._individual_island)

    @property
    def islands(self):
        return self._islands

def islands_from_string(tree):
    if isinstance(tree, str):
        tree = tree.split("\n")
    islands = []
    for line in tree:
        line = line.strip()
        if line.startswith("#"):
            # Ignore comments in tree file.
            continue
        rate = float(line)
        islands.append(IslandNode(rate))
    return IslandModel(islands)

def islands_from_file(f):
    with open(f, "r") as tree_file:
        return islands_from_string(tree_file.readlines())

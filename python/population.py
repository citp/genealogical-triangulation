from random import shuffle, uniform, choice, sample
from itertools import chain
from types import GeneratorType
from pickle import Unpickler

from generation import Generation
from sex import Sex

class Population:
    def __init__(self, initial_generation = None):
        self._generations = []
        self._kinship_coefficients = None
        self._node_to_generation = (None, -1)
        # The last n generations with genomes defined
        self._generations_with_genomes = None
        if initial_generation is not None:
            self._generations.append(initial_generation)

    @property
    def id_mapping(self):
        """
        Returns the node_id -> node mapping.
        Assumes all nodes use the same mapping
        """
        # Return any nodes mapping, because they should all share the
        # same mapping
        return self._generations[0].members[0].mapping
        

    @property
    def generations(self):
        return list(self._generations)

    @property
    def members(self):
        return chain.from_iterable(generation.members for
                                   generation in self._generations)

    @property
    def size(self):
        return sum(generation.size for generation in self._generations)

    @property
    def node_to_generation(self):
        """
        Maps nodes in this population to their generation.
        Higher numbered generations are more recent.
        """
        # Cache the results of this function so it only needs to be
        # computed when new generations are added.
        cached = self._node_to_generation
        if cached[0] is not None and cached[1] == len(self._generations):
            return cached[0]        
        generations = [generation.members for generation in self._generations]
        node_to_generation = dict()
        for generation_num, members in enumerate(generations):
            for node in members:
                node_to_generation[node] = generation_num
        self._node_to_generation = (node_to_generation, len(self._generations))
        return node_to_generation

    def clean_genomes(self):
        for person in self.members:
            person.genome = None


    @property
    def num_generations(self):
        """
        Return the number of generations
        """
        return len(self._generations)

    def _symmetric_members(self):
        """
        Generator that yields 2 item tuples that contain members of
        the popluation. Generates all pairs of members, where members
        are visited in the first entry of the tuple first.
        TODO: Can this be replaced by itertools.combinations_with_replacement?
        """
        members = list(self.members)
        return ((members[y], members[x])
                for x in range(len(members))
                for y in range(x, len(members)))


class IslandPopulation(Population):
    """
    A population where mates are selected based on
    locality. Individuals exists on "islands", and will search for
    mates from a different island with a given switching probability.
    """
    def __init__(self, island_model, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._island_model = island_model
        # Assume existing members of the tree are a part of the founders
        if len(island_model.individuals) > 0:
            assert len(self._generations) is 0
            self._generations.append(Generation(island_model.individuals))


    def _pick_island(self, individual):
        current_island = self._island_model.get_island(individual)
        if uniform(0, 1) < current_island.switch_probability:
            candidate_islands = set(self._island_model.islands)
            candidate_islands.remove(current_island)
            return choice(tuple(candidate_islands))
        return current_island

    def _island_members(self, generation_members):
        """
        Returns a dictionary mapping islands to sets of individuals
        from generation_members.
        """
        if isinstance(generation_members, GeneratorType):
            generation_members = list(generation_members)
        members = dict()
        for leaf in self._island_model.islands:
            m = leaf.individuals.intersection(generation_members)
            members[leaf] = m
        return members

    def migrate_generation(self, generation):
        """
        Cause all members of the population to migrate.
        """
        for member in generation.members:
            island = self._pick_island(member)
            self._island_model.move_individual(member, island)
            

    def new_generation(self, node_generator, size = None,
                       non_paternity_rate = 0, adoption_rate = 0,
                       unknown_mother_rate = 0, unknown_father_rate =0,
                       monogamy_rate = 0.0):
        """
        Generates a new generation of individuals from the previous
        generation. If size is not passed, the new generation will be
        the same size as the previous generation.
        """
        assert 0.0 <= monogamy_rate <= 1.0
        if size is None:
            size = self._generations[-1].size
        self.migrate_generation(self._generations[-1])
        previous_generation = set(self._generations[-1].members)

        monogamous_map = dict()
        random_map = dict()
        for island in self._island_model.islands:
            members = list(previous_generation.intersection(island.individuals))
            shuffle(members)
            men = [x for x in members if x.sex == Sex.Male]
            women = [x for x in members if x.sex == Sex.Female]

            monogamy_cutoff = min(int(monogamy_rate * len(men)),
                                  int(monogamy_rate * len(women)))
            monogamous_men = men[:monogamy_cutoff]
            non_monogamous_men = men[monogamy_cutoff:]


            monogamous_women = women[:monogamy_cutoff]
            non_monogamous_women = women[monogamy_cutoff:]

            for male, female in zip(monogamous_men, monogamous_women):
                monogamous_map[male] = female
                monogamous_map[female] = male

            for male in non_monogamous_men:
                random_map[male] = non_monogamous_women

            for female in non_monogamous_women:
                random_map[female] = non_monogamous_men

        new_nodes = []
        num_twins = int(0.003 * size)
        pre_twins_size = size - num_twins
        previous_generation_list = list(previous_generation)
        while len(new_nodes) < pre_twins_size:
            parent_a = choice(previous_generation_list)
            if parent_a in monogamous_map:
                parent_b = monogamous_map[parent_a]
            else:
                candidate_parents = random_map[parent_a]
                if len(candidate_parents) == 0:
                    continue
                parent_b = choice(candidate_parents)

            father, mother = _sort_sex(parent_a, parent_b)
            child = node_generator.generate_node(father, mother)
            island = self._island_model.get_island(father)
            self._island_model.add_individual(island, child)
            new_nodes.append(child)
        
        apportioned = apportion(new_nodes, non_paternity_rate, adoption_rate)
        non_paternity, adopted = apportioned
        already_error = non_paternity.union(adopted)
        no_error = list(set(new_nodes) - already_error)

        # We don't use apportion here, because they are not mutually exclusive
        unknown_mother = sample(no_error,
                                int(unknown_mother_rate * len(new_nodes)))
        unknown_father = sample(no_error,
                                int(unknown_father_rate * len(new_nodes)))
        for node in unknown_mother:
            node.set_suspected_mother(None)
        for node in unknown_father:
            node.set_suspected_father(None)

        men_by_island = self._island_members(self._generations[-1].men)
        men_by_island = {island: list(men) for island, men
                         in men_by_island.items()}
        for node in non_paternity:
            child_island = self._island_model.get_island(node)
            suspected_father = choice(men_by_island[child_island])
            node.set_suspected_father(suspected_father)

        women_by_island = self._island_members(self._generations[-1].women)
        women_by_island = {island: list(women) for island, women
                           in women_by_island.items()}
        for node in adopted:
            child_island = self._island_model.get_island(node)
            suspected_father = choice(men_by_island[child_island])
            node.set_suspected_father(suspected_father)
            suspected_mother = choice(women_by_island[child_island])
            node.set_suspected_mother(suspected_mother)

        # Generate twins. Twins will essentially be copies of their sibling
        for template_node in sample(new_nodes, num_twins):
            twin = node_generator.twin_node(template_node)
            island = self._island_model.get_island(template_node)
            self._island_model.add_individual(island, twin)
            new_nodes.append(twin)


                
        SIZE_ERROR = "Generation generated is not correct size. Expected {}, got {}."
        assert len(new_nodes) == size, SIZE_ERROR.format(size, len(new_nodes))
        self._generations.append(Generation(new_nodes))

    @property
    def island_tree(self):
        return self._island_tree

def apportion(original, *rates):
    remainder = list(original)
    assigned = set()
    length = len(original)
    subsets = []
    for rate in rates:
        assert int(length * rate) < len(remainder)
        if rate == 0:
            current_group = set()
        else:
            current_group = set(sample(remainder, int(length * rate)))
        assigned.update(current_group)
        subsets.append(current_group)
        if 0 < len(current_group):
            remainder = [node for node in original if node not in assigned]
        
    return tuple(subsets)

def _sort_sex(a, b):
    if a.sex == Sex.Male:
        assert b.sex == Sex.Female
        return (a, b)
    else:
        assert b.sex == Sex.Male
        return (b, a)
        

def fix_twin_parents(population):
    """
    There was a bug in the inital implementation of twins such that
    their suspected parents would not have the both twins in the list
    of suspected children. This didn't impact performance, as the
    children property is mostly for easy of analysis and debugging.

    This method sets all members suspected children list to what is
    indicated by each nodes suspected mother/father property.
    """
    for member in population.members:
        member._suspected_children = []
        member.set_suspected_mother(member.suspected_mother)
        member.set_suspected_father(member.suspected_father)
                    

class PopulationUnpickler(Unpickler):
    
    def load(self):
        result = super().load()
        for member in result.members:
            member._resolve_parents()
        return result

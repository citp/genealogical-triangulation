from math import floor
from random import shuffle, uniform, choice, random, sample
from itertools import chain
from types import GeneratorType
from pickle import Unpickler
from collections import defaultdict

import json

from random_queue import RandomQueue
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


class HierarchicalIslandPopulation(Population):
    """
    A population where mates are selected based on
    locality. Individuals exists on "islands", and will search for
    mates from a different island with a given switching probability.
    """
    def __init__(self, island_tree, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._island_tree = island_tree
        # Assume existing members of the tree are a part of the founders
        if len(island_tree.individuals) > 0:
            assert len(self._generations) is 0
            self._generations.append(Generation(island_tree.individuals))


    def _pick_island(self, individual):
        island = self._island_tree.get_island(individual)
        # Traverse upward
        while uniform(0, 1) < island.switch_probability:
            island = island.parent
            if island.parent is None:
                break
        # Now drill down into final island.
        while not island.is_leaf:
            island = choice(list(island.islands))
        return island

    def _island_members(self, generation_members):
        """
        Returns a dictionary mapping islands to sets of individuals
        from generation_members.
        """
        if isinstance(generation_members, GeneratorType):
            # If generation_members is a generator, then all the
            # elements are "used up" on the first leaf. We need to
            # iterate over the elements of generation_members multiple
            # time.
            generation_members = list(generation_members)
        members = dict()
        for leaf in self._island_tree.leaves:
            m = leaf.individuals.intersection(generation_members)
            members[leaf] = m
        return members

    def migrate_generation(self, generation):
        """
        Cause all members of the population to migrate.
        """
        for member in generation.members:
            island = self._pick_island(member)
            self._island_tree.move_individual(member, island)
            

    def new_generation(self, size = None, non_paternity_rate = 0,
                       adoption_rate = 0, unknown_mother_rate = 0,
                       unknown_father_rate = 0, multi_partner_probs = None):
        """
        Generates a new generation of individuals from the previous
        generation. If size is not passed, the new generation will be
        the same size as the previous generation.
        """
        if size is None:
            size = self._generations[-1].size
        self.migrate_generation(self._generations[-1])
        previous_generation = list(self._generations[-1].members)
        shuffle(previous_generation)
        boundary = len(previous_generation) // 2
        seekers = RandomQueue(previous_generation[:boundary])
        mates = set(previous_generation[boundary:])
        mate_set = defaultdict(set)
        mate_attempts = defaultdict(int)
        available_mates = self._island_members(mates)
        pairs = []
        while len(seekers) > 0:
            seeker = seekers.dequeue()
            if sum(len(members) for members in available_mates.values()) is 0:
                break
            island = self._island_tree.get_island(seeker)
            potential_mates = [mate for mate in available_mates[island]
                               if mate.sex != seeker.sex]
            shuffle(potential_mates)
            for potential_mate in potential_mates:
                share_mother = (seeker.mother is not None and
                                potential_mate.mother == seeker.mother)
                share_father = (seeker.father is not None and
                                potential_mate.father == seeker.father)
                if not share_mother and not share_father:
                    # TODO insert these pairs with male node first.
                    pairs.append((seeker, potential_mate))
                    existing_mates = mate_set[seeker]
                    existing_mates.add(potential_mate)
                    num_mates = len(existing_mates)
                    if (multi_partner_probs is not None and
                        multi_partner_probs[num_mates] < random()):
                        seekers.enqueue(seeker)
                    else:
                        available_mates[island].remove(potential_mate)
                    break
            mate_attempts[seeker] += 1
        num_twins = int(0.003 * size)
        non_twin_size = size - num_twins
        min_children = floor(non_twin_size / len(pairs))
        # Number of families with 1 more than the min number of
        # children. Because only having 2 children per pair only works
        # if there is an exact 1:1 ratio of men to women.
        extra_child = non_twin_size - min_children * len(pairs)
        node_generator = pairs[0][0].node_generator
        new_nodes = []
        for i, (seeker, mate) in enumerate(pairs):
            if seeker.sex == Sex.Male:
                man = seeker
                woman = mate
            else:
                man = mate
                woman = seeker

            if i < extra_child:
                extra = 1
            else:
                extra = 0
                
            # Child will be based at mate's island
            island = self._island_tree.get_island(mate)
            for i in range(min_children + extra):
                child = node_generator.generate_node(man, woman)
                new_nodes.append(child)
                self.island_tree.add_individual(island, child)

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
            child_island = self._island_tree.get_island(node)
            suspected_father = choice(men_by_island[child_island])
            node.set_suspected_father(suspected_father)

        women_by_island = self._island_members(self._generations[-1].women)
        women_by_island = {island: list(women) for island, women
                           in women_by_island.items()}
        for node in adopted:
            child_island = self._island_tree.get_island(node)
            suspected_father = choice(men_by_island[child_island])
            node.set_suspected_father(suspected_father)
            suspected_mother = choice(women_by_island[child_island])
            node.set_suspected_mother(suspected_mother)

        # Generate twins. Twins will essentially be copies of their sibling
        for template_node in sample(new_nodes, num_twins):
            twin = node_generator.twin_node(template_node)
            island = self._island_tree.get_island(template_node)
            self._island_tree.add_individual(island, twin)
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


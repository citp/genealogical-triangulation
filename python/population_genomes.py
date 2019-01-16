from bisect import bisect_left
from collections import deque
from random import random

import numpy as np

from sex import Sex
from recomb_genome import RecombGenome, Diploid, CHROMOSOME_ORDER

def _pick_chroms_for_diploid(genome, recombinator):
    """
    Takes a genome and returns a diploid chromosome that is the result
    of recombination events and randomly picking a diploid for each
    autosome.
    """
    recomb_genome = recombinator.recombination(genome)
    mother = recomb_genome.mother
    father = recomb_genome.father
    starts = []
    founder = []
    offsets = recombinator._chrom_start_offset
    chrom_stop = 0
    for chrom_name in CHROMOSOME_ORDER:
        if random() < 0.5:
            tmp_diploid = mother
        else:
            tmp_diploid = father

        chrom_start = bisect_left(tmp_diploid.starts, offsets[chrom_name])
        if chrom_name != CHROMOSOME_ORDER[-1]:
            chrom_stop = bisect_left(tmp_diploid.starts,
                                     offsets[chrom_name + 1],
                                     chrom_start)
        else:
            chrom_stop = len(tmp_diploid.starts)
        starts.extend(tmp_diploid.starts[chrom_start:chrom_stop])
        founder.extend(tmp_diploid.founder[chrom_start:chrom_stop])
    return Diploid(np.array(starts, dtype = np.uint32),
                   mother.end,
                   np.array(founder, dtype = np.uint32))


def mate(mother, father, mother_recombinator, father_recombinator):
    """
    Takes a mother and father, and returns a genome for a child
    """
    assert mother is not None
    assert father is not None
    from_mother = _pick_chroms_for_diploid(mother, mother_recombinator)
    from_father = _pick_chroms_for_diploid(father, father_recombinator)
    return RecombGenome(from_mother, from_father)

def generate_genomes_ancestors(root_nodes, generator, recombinators):
    queue = deque(root_nodes)
    visited = set()
    while len(queue) > 0:
        person = queue.popleft()
        if person in visited:
            continue
        if person.mother is not None:
            assert person.father is not None
            if person.mother not in visited or person.father not in visited:
                continue
            person.genome =  mate(person.mother.genome, person.father.genome,
                                  recombinators[Sex.Female],
                                  recombinators[Sex.Male])
        else:
            person.genome = generator.generate()
        queue.extend(person.children)
        visited.add(person)

def generate_genomes(population, generator, recombinators, keep_last = None,
                     true_genealogy = True):
    assert keep_last is None or keep_last > 0
    for generation_num, generation in enumerate(population.generations):
        for person in generation.members:
            if person.genome is not None:
                continue
            if person.twin is not None and person.twin.genome is not None:
                person.genome = person.twin.genome
                continue
            if true_genealogy:
                mother = person.mother
                father = person.father
            else:
                mother = person.suspected_mother
                father = person.suspected_father
            
            if mother is None and father is None:
                person.genome = generator.generate()
                continue
            if mother is None:
                mother_genome = generator.generate()
            else:
                mother_genome = mother.genome
            if father is None:
                father_genome = generator.generate()
            else:
                father_genome = father.genome
                
            assert mother_genome is not None
            assert father_genome is not None
            person.genome = mate(mother_genome, father_genome,
                                 recombinators[Sex.Female],
                                 recombinators[Sex.Male])
        if keep_last is not None and keep_last <= generation_num:
            to_delete = population.generations[generation_num - keep_last]
            for person in to_delete.members:
                person.genome = None

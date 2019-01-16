from itertools import combinations
from random import sample
from time import perf_counter

import threading
import pdb

import numpy as np
# import pyximport; pyximport.install()

from population import PopulationUnpickler
from classify_relationship import shared_segment_length_genomes

print("Loading population")
with open("population_10000.pickle", "rb") as pickle_file:
    population = PopulationUnpickler(pickle_file).load()


def calculate_to_list(pairs, lengths):
    lengths.extend(shared_segment_length_genomes(node_a.genome, node_b.genome,
                                                 0)
                   for node_a, node_b in pairs)

                   
# print("Comparing pairs.")
# nodes = population.generations[-1].members
# nodes = sample(nodes, 1500)
# pairs = list(combinations(nodes, 2))
# lengths = []
# start = perf_counter()
# boundary = len(pairs) // 2
# thread_1 = threading.Thread(target=calculate_to_list,
#                             args = (pairs[:boundary], lengths))
# thread_2 = threading.Thread(target=calculate_to_list,
#                             args = (pairs[boundary:], lengths))
# thread_1.start()
# thread_2.start()
# thread_1.join()
# thread_2.join()
# stop = perf_counter()
# print(stop - start)

print("Comparing pairs.")
nodes = population.generations[-1].members
nodes = sample(nodes, 1500)
start = perf_counter()
lengths = [shared_segment_length_genomes(node_a.genome, node_b.genome, 0)
           for node_a, node_b in combinations(nodes, 2)]
stop = perf_counter()
print(stop - start)

# import pdb
# pdb.set_trace()
# shared = [len(np.flatnonzero(np.unpackbits(a.genome._founder_bits & b.genome._founder_bits)))
#           for a, b in combinations(nodes, 2)]

# print(np.average(shared))
# print(np.std(shared))
# print(max(lengths))

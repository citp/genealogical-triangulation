#!/usr/bin/env python3.5

import pickle

import matplotlib.pyplot as plt
import numpy as np

from population_statistics import proportion_within_distance, cdf_num_labeled_within_distance, shared_segment_distribution

# POPULATION_FILE = "population_no_genome_100000.pickle"
POPULATION_FILE = "population_50000.pickle"
RELATION_DISTANCE = 5

with open(POPULATION_FILE, "rb") as pickle_file:
    population = pickle.load(pickle_file)


lengths = np.array(shared_segment_distribution(population, 5))
weights = np.ones_like(lengths)/float(len(lengths))
# plt.hist(lengths, bins = 20, normed = True)
plt.hist(lengths, bins = 20, weights = weights)
# plt.xscale('log')
plt.show()

# counts = cdf_num_labeled_within_distance(population, RELATION_DISTANCE, 0.02)
# plt.hist(counts, bins = 2 ** (RELATION_DISTANCE - 1),
#          normed = True, cumulative = True)
# plt.xlabel("Number of labeled relatives within distance {}".format(RELATION_DISTANCE))
# plt.ylabel("CDF of people with labeled relatives.")
# plt.show()

# labeled_fraction = np.arange(0.001, 0.01, 0.0001)
# percent_covered = []
# for percent in labeled_fraction:
#     covered = np.average([proportion_within_distance(population,
#                                                      RELATION_DISTANCE,
#                                                      percent)
#                           for _ in range(5)])
#     percent_covered.append(covered)

# plt.plot(labeled_fraction, percent_covered)
# plt.xlabel("Fraction labeled")
# plt.ylabel("Fraction within {} generations. ({}th cousin or less)".format(RELATION_DISTANCE, RELATION_DISTANCE - 1))
# # axes = plt.gca()
# # axes.set_xlim([0.0, 0.015])
# # axes.set_ylim([0.0, 1.0])
# plt.show()


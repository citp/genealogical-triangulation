#!/usr/bin/env python3

import pdb

from pympler import asizeof, tracker, summary, muppy

from model import Generation, Population, generate_population

print("Generating population")
population = generate_population(100000)

print("Creating 10 generations.")

for _ in range(10):
    population.new_generation()


for i, generation in enumerate(population._generations):
    if i is 0:
        continue
    previous_generation_members = set(population._generations[i-1].members)
    for member in generation.members:
        assert member.mother in previous_generation_members
        assert member.father in previous_generation_members


# print("Calculating kinship coefficients")
# tr = tracker.SummaryTracker()
# kinship = population.kinship_coefficients
# tr.print_diff()

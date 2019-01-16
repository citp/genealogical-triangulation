#!/usr/bin/env python3

from random import choice

from island_model import tree_from_file
from population import HierarchicalIslandPopulation
from node import Node

# from pympler import tracker

print("Reading tree from file.")
tree = tree_from_file("test_tree")

num_leaves = len(tree.leaves)
assert num_leaves is 3, "Expected 3 leaves, got {}".format(num_leaves)

print("Creating population.")
# tr = tracker.SummaryTracker()

initial_population = [Node() for _ in range(1000000)]

leaves = tree.leaves
for person in initial_population:
    tree.add_individual(choice(leaves), person)

population = HierarchicalIslandPopulation(tree)

print("Creating 10 generations.")
for i in range(10):
    print("Creating generation {}".format(i))
    population.new_generation()
    

# tr.print_diff()

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

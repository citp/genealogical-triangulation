#!/usr/bin/env python3

from random import sample, seed
from pickle import load
from argparse import ArgumentParser

import pdb

from bayes_deanonymize import BayesDeanonymize
from population import PopulationUnpickler

parser = ArgumentParser(description = "Evaluate performance of classification.")
parser.add_argument("population")
parser.add_argument("classifier")
parser.add_argument("output_file")
parser.add_argument("--trials", type = int, default = 50)
args = parser.parse_args()


print("Loading population.")
with open(args.population, "rb") as pickle_file:
    population = PopulationUnpickler(pickle_file).load()

print("Loading classifier")
with open(args.classifier, "rb") as pickle_file:
    classifier = load(pickle_file)

nodes = set(member for member in population.members
             if member.genome is not None)

bayes = BayesDeanonymize(population, classifier)

id_mapping = population.id_mapping
labeled_nodes = set(id_mapping[node_id] for node_id
                    in classifier._labeled_nodes)
labeled_node_ids = list(classifier._labeled_nodes)

unlabeled = sample(list(nodes - labeled_nodes), args.trials)

for labeled_size in range(20, len(labeled_nodes) + 1, 5):
    classifier._labeled_nodes = labeled_node_ids[:labeled_size]
    correct = 0
    incorrect = 0
    for i, node in enumerate(unlabeled):
        identified = bayes.identify(node.genome, node, population, 0.03)
        if node in identified:
            correct += 1
        else:
            incorrect += 1
    with open(args.output_file, "a") as output_file:
        output_file.write("{}\t{}\n".format(labeled_size,
                                            (correct / len(unlabeled))))
    print("Calculated for {} labeled nodes.".format(labeled_size))

#!/usr/bin/env python3

from random import sample
from pickle import load
from argparse import ArgumentParser
from util import recent_common_ancestor

import pdb

from scipy import stats

from bayes_deanonymize import BayesDeanonymize
from population import PopulationUnpickler

parser = ArgumentParser(description = "Evaluate performance of classification.")
parser.add_argument("population")
parser.add_argument("classifier")
parser.add_argument("--num_node", "-n", type = int, default = 10)
parser.add_argument("--start", "-st", type = int, default = 5,
                    help = "Number of labeled nodes to start with.")
parser.add_argument("--stride", "-s", type = int, default = 5,
                    help = "Amount to increase labeled nodes by")
parser.add_argument("--stop", type = int, default = None,
                    help = "The maximum number of labeled nodes to identify with.")
args = parser.parse_args()


print("Loading population.")
with open(args.population, "rb") as pickle_file:
    population = PopulationUnpickler(pickle_file).load()

print("Loading classifier")
with open(args.classifier, "rb") as pickle_file:
    classifier = load(pickle_file)

nodes = set(member for member in population.members
             if member.genome is not None)


def evaluate(unlabeled, bayes):
    correct = 0
    incorrect = 0
    for i, node in enumerate(unlabeled):
        identified = bayes.identify(node.genome, node, population)
        if node in identified:
            correct += 1
        else:
            incorrect += 1
    return correct / len(unlabeled)

bayes = BayesDeanonymize(population, classifier)
id_mapping = population.id_mapping
labeled_nodes = set(id_mapping[node_id] for node_id
                    in classifier._labeled_nodes)
unlabeled = sample(list(nodes - labeled_nodes),
                   args.num_node)

labeled_copy = list(classifier._labeled_nodes)
if args.stop is None:
    stop = len(classifier._labeled_nodes)
else:
    stop = min(len(classifier._labeled_nodes), args.stop)

print("# labeled nodes\taccuracy")
for i in range(args.start, stop + args.stride, args.stride):
    classifier._labeled_nodes = labeled_copy[:i]
    accuracy = evaluate(unlabeled, bayes)
    print("{}\t{}".format(i, accuracy))

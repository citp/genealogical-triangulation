from argparse import ArgumentParser
from random import sample
from itertools import chain

from population import PopulationUnpickler
from classify_relationship import related_pairs
from to_json import to_json

parser = ArgumentParser(description = "Serialize a population for simulation in Rust. Anchor nodes and adversary horizon are picked before serialization")
parser.add_argument("population_file", help = "Pickled file with population")
parser.add_argument("output", help = "Name of file for json output")
parser.add_argument("--gen-back", "-g", type = int, default = 6,
                    help = "Ignore common ancestry more than the given number of generations back. This sets the analyst horizion.")
parser.add_argument("--num-anchor-nodes", "-n", type = int, default = 0,
                    help = "Number of nodes to include in the 'known' set. Nodes are picked randomly from the last generation. If no value is given, the number is the population size multiplied by 0.01. (Note that this is not the same as a random sample of 1%% of the entire population as only the last generation is used.)")
args = parser.parse_args()

print("Loading population")
with open(args.population_file, "rb") as pickle_file:
    population = PopulationUnpickler(pickle_file).load()

potentially_labeled = list(chain.from_iterable([generation.members
                                                for generation
                                                in population.generations[-3:]]))
if args.num_anchor_nodes <= 0:
    num_labeled_nodes = len(potentially_labeled) // 100
else:
    num_labeled_nodes = args.num_anchor_nodes
labeled_nodes = sample(potentially_labeled, num_labeled_nodes)

num_generations = population.num_generations
clear_index = max(num_generations - args.gen_back, 0)
to_clear = population.generations[clear_index].members
for node in to_clear:
    node.suspected_mother = None
    node.suspected_mother_id = None
    node.suspected_father = None
    node.suspected_father_id = None
unlabeled_nodes = set(potentially_labeled)
print("Computing related pairs.")
related_nodes = related_pairs(unlabeled_nodes, labeled_nodes, population,
                              args.gen_back)
print("Generating serialization.")
json = to_json(population, labeled_nodes, related_nodes)
print("Writing json file.")
with open(args.output, "w") as json_file:
    json_file.write(json)

#!/usr/bin/env python3

from argparse import ArgumentParser
from random import sample
from itertools import chain
from pickle import dump
from os import listdir
from os.path import dirname, join, abspath

from population import PopulationUnpickler, fix_twin_parents
from sex import Sex
from classify_relationship import generate_classifier, related_pairs
from recomb_genome import recombinators_from_directory, RecombGenomeGenerator
from cm import centimorgan_data_from_directory
from to_json import to_json

parser = ArgumentParser(description = "Generate a classifier which can (hopefully) identify individuals in a population.")
parser.add_argument("population_file", help = "Pickled file with population")
parser.add_argument("work_dir",
                    help = "Directory to put shared length calculations in.")
parser.add_argument("num_iterations", type = int, default = 1000,
                    help = "Number of samples to collect from empirical distributions")
parser.add_argument("--recover", "-r", default = False,
                    action="store_true",
                    help = "work directory from interrupted run.")
parser.add_argument("--gen_back", "-g", type = int, default = 7,
                    help = "Ignore common ancestry more than the given number of generations back.")
parser.add_argument("--num_labeled_nodes", "-n", type = int, default = 0,
                    help = "Number of nodes to include in the 'known' set. Nodes are picked randomly from the last generation. If no value is given, the number is the population size multiplied by 0.01. (Note that this is not the same as a random sample of 1%% of the entire population as only the last generation is used.)")
parser.add_argument("--output_pickle", default = "distributions.pickle",
                    help = "File to store distributions in. Pickle format will be used. Default is 'distributions.pickle'")
parser.add_argument("--non_paternity", "-np", type = float, default = 0.0,
                    help = "Non paternity rate for the adversary to assume.")
parser.add_argument("--to_json", default = None,
                    help = "If this flag is present, will instead store the population as json for faster computation in another language")

args = parser.parse_args()

print("Loading population")
with open(args.population_file, "rb") as pickle_file:
    population = PopulationUnpickler(pickle_file).load()
fix_twin_parents(population)

if not args.recover:
    potentially_labeled = list(chain.from_iterable([generation.members
                                                    for generation
                                                    in population.generations[-3:]]))
    if args.num_labeled_nodes <= 0:
        num_labeled_nodes = population.size // 100
    else:
        num_labeled_nodes = args.num_labeled_nodes
    labeled_nodes = sample(potentially_labeled, num_labeled_nodes)
else:
    print("Recovering run")
    labeled_nodes = [population.id_mapping[int(filename)]
                     for filename in listdir(args.work_dir)]

if args.to_json:
    num_generations = population.num_generations
    clear_index = max(num_generations - args.gen_back, 0)
    to_clear = population.generations[clear_index].members
    for node in to_clear:
        node.suspected_mother = None
        node.suspected_mother_id = None
        node.suspected_father = None
        node.suspected_father_id = None
    unlabeled_nodes = set(chain.from_iterable(generation.members
                                          for generation
                                          in population.generations[-3:]))
    related_nodes = related_pairs(unlabeled_nodes, labeled_nodes, population,
                                  args.gen_back)
    json = to_json(population, labeled_nodes, related_nodes)
    with open(args.to_json, "w") as json_file:
        json_file.write(json)
    exit()

# TODO: Add a command line option for this
recomb_dir = abspath(join(dirname(__file__), "../data/recombination_rates/"))
print("Loading recombination data.")
recombinators = recombinators_from_directory(recomb_dir)
chrom_sizes = recombinators[Sex.Male]._num_bases
genome_generator = RecombGenomeGenerator(chrom_sizes)


print("Populating length classifier.")

clobber = not (args.recover or args.num_iterations == 0)

classifier = generate_classifier(population, labeled_nodes,
                                 genome_generator, recombinators,
                                 args.work_dir,
                                 iterations = args.num_iterations,
                                 clobber = clobber,
                                 generations_back_shared = args.gen_back,
                                 min_segment_length = 5000000,
                                 non_paternity = args.non_paternity)

del recombinators
del labeled_nodes
del genome_generator
del population

# import pdb
# pdb.set_trace()

print("Pickling classifier")
with open(args.output_pickle, "wb") as pickle_file:
    dump(classifier, pickle_file)

print("Pickling complete")
# import pdb
# pdb.set_trace()

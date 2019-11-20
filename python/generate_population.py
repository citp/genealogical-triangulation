#!/usr/bin/env python3

from argparse import ArgumentParser
from random import choice
from pickle import dump, HIGHEST_PROTOCOL
from os.path import exists
from sys import exit

from population import IslandPopulation
from population_genomes import generate_genomes
from node import NodeGenerator
from recomb_genome import recombinators_from_directory, RecombGenomeGenerator
from island_model import islands_from_file
from sex import Sex

def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def generation_sizes(final_size, avg_children, num_generations):
    cur = int(final_size)
    ret = []
    for _ in range(num_generations):
        ret.append(cur)
        cur = round(cur * (2 / avg_children))
    ret.reverse()
    return ret

parser = ArgumentParser(description = "Generate a population and its associated genomes.")
parser.add_argument("island_file",
                    help = "Describes island model.")
parser.add_argument("recombination_dir",
                    help = "Directory containing Hapmap and decode data.")
parser.add_argument("--generation-size", help = "Number of individuals in the final generation. If avg-children is not used, this is the size of every generation.",
                    type = int, default = 100000)
parser.add_argument("--num_generations", type = int, default = 10)
parser.add_argument("--no_genomes", action="store_true", default = False,
                    help = "Don't generate genomes for the individuals in the population.")
parser.add_argument("--non_paternity", "-p", type = float, default = 0,
                    help = "Rate with which the suspected father is not the true father.")
parser.add_argument("--adoption", "-a", type = float, default = 0,
                    help = "Rate with which the suspected mother and father are not the true mother and father.")
parser.add_argument("--missing_mother", "-mm", type = float, default = 0,
                    help = "Rath with which the mother is unknown.")
parser.add_argument("--missing_father", "-mf", type = float, default = 0,
                    help = "Rate with which the father is unknown.")
parser.add_argument("--monogamy-rate", type = float, default = 0.0,
                    help = "Fraction of the population that is monogamous.")
parser.add_argument("--avg-children", type = float, default = None,
                    help = "Average number of children per individual in each generation. Default results in constant sized generations.")


parser.add_argument("--output_file", default = "population.pickle",
                    help = "Outputs a pickle file containing a Population object to this file. This file will be clobbered if it exists.")

args = parser.parse_args()
if args.num_generations < 1:
    parser.error("num_generations must be >= 1")

if not 0 <= args.non_paternity <= 1:
    parser.error("Non-paternity rate must be in the range [0, 1]")
    
if not 0 <= args.adoption <= 1:
    parser.error("adoption rate must be in the range [0, 1]")

if args.avg_children is not None:
    sizes = generation_sizes(args.generation_size, args.avg_children,
                             args.num_generations)
else:
    sizes = [args.generation_size for _ in range(args.num_generations)]

node_generator = NodeGenerator()
print("Generating founders")
founders = [node_generator.generate_node() for _ in range(sizes[0])]

island_model = islands_from_file(args.island_file)
islands = island_model.islands
for person in founders:
    island_model.add_individual(choice(islands), person)
population = IslandPopulation(island_model)

print("Adding more generations")
for generation_size in sizes[1:]:
    population.new_generation(node_generator, size = generation_size,
                              non_paternity_rate = args.non_paternity,
                              adoption_rate = args.adoption,
                              unknown_mother_rate = args.missing_mother,
                              unknown_father_rate = args.missing_father,
                              monogamy_rate = args.monogamy_rate)

if not args.no_genomes:
    print("Loading recombination rates")
    recombinators = recombinators_from_directory(args.recombination_dir)
    chrom_sizes = recombinators[Sex.Male]._num_bases
    genome_generator = RecombGenomeGenerator(chrom_sizes)
    print("Generating genomes")
    generate_genomes(population, genome_generator, recombinators, 3)


if args.output_file:
    if exists(args.output_file):
        print("Output file already exists.")
        exit(1)

    with open(args.output_file, "wb") as pickle_file:
        # Trees cause deep recursion in the pickle module, so we need
        # to raise the recursion limit. This is the stack depth for
        # python functions, you may need to increase the native stack
        # depth using ulimit -s
        # https://docs.python.org/3.4/library/pickle.html#what-can-be-pickled-and-unpickled
        print("Saving population file to {}".format(args.output_file))
        dump(population, pickle_file, protocol = HIGHEST_PROTOCOL)

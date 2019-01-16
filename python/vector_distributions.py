from collections import deque, OrderedDict
from argparse import ArgumentParser
from os.path import basename
import pdb


import numpy as np
from matplotlib import pyplot

from node import NodeGenerator
from sex import Sex
from population_genomes import generate_genomes_ancestors
from recomb_genome import recombinators_from_directory, RecombGenomeGenerator
from classify_relationship import shared_segment_length_genomes

node_generator = NodeGenerator()
name_map = None

def parse_genealogy_file(filename):
    genders = dict()
    parents = OrderedDict()
    reading_names = True
    with open(filename, "r") as genealogy_file:
        for line in genealogy_file:
            line = line.strip()
            if len(line) is 0:
                continue
            if line == "#":
                reading_names = False
                continue
            line = line.split()
            name = line[0]
            if reading_names:
                if len(line) > 1:
                    genders[name] = Sex[line[1].lower().capitalize()]
                else:
                    genders[name] = None
            else:
                father = line[1]
                mother = line[2]
                assert father in genders
                assert mother in genders
                parents[name] = (father, mother)
    founder_names = set(genders.keys()) - parents.keys()
    founders = {name : node_generator.generate_node(sex = genders[name])
                for name in founder_names}
    all_nodes = dict(founders)
    for name, (father_name, mother_name) in parents.items():
        all_nodes[name] = node_generator.generate_node(all_nodes[father_name],
                                                       all_nodes[mother_name],
                                                       sex = genders[name])
    return (all_nodes, set(founders.values()))

def simulate_sharing(founders, pair, genome_generator, recombinators,
                     iterations = 10000):
    sharing = []
    for i in range(iterations):
        generate_genomes_ancestors(founders, genome_generator, recombinators)
        shared = shared_segment_length_genomes(pair[0].genome, pair[1].genome,
                                               0)
        sharing.append(shared)
    return sharing
        
parser = ArgumentParser(description = "Examine the distributions of various relationships.")
parser.add_argument("relationship_file", nargs = 2,
                    help = "File to describe the genealogy.")
parser.add_argument("pair_1", nargs = 2)
parser.add_argument("pair_2", nargs = 2)

args = parser.parse_args()

all_nodes_0, founders_0 = parse_genealogy_file(args.relationship_file[0])
name_map = {node: node_name for node_name, node in all_nodes_0.items()}
all_nodes_1, founders_1 = parse_genealogy_file(args.relationship_file[1])

recombinators = recombinators_from_directory("../data/recombination_rates")
chrom_sizes = recombinators[Sex.Male]._num_bases
genome_generator = RecombGenomeGenerator(chrom_sizes)

pair_0 = [all_nodes_0[name] for name in args.pair_1]
sharing_0 = simulate_sharing(founders_0, pair_0, genome_generator,
                             recombinators)

pair_1 = [all_nodes_1[name] for name in args.pair_2]
sharing_1 = simulate_sharing(founders_1, pair_1, genome_generator,
                             recombinators)

weights_0 = np.ones_like(sharing_0)/float(len(sharing_0))
weights_1 = np.ones_like(sharing_1)/float(len(sharing_1))
# pyplot.hist(sharing_0, bins = 30, alpha = 0.5, normed = False,
#             label = basename(args.relationship_file[0]))
# pyplot.hist(sharing_1, bins = 30, alpha = 0.5, normed = False,
#             label = basename(args.relationship_file[1]))
pyplot.hist(sharing_0, bins = 30, normed = False,
            weights = weights_0,
            label = basename(args.relationship_file[0]))
pyplot.hist(sharing_1, bins = 30, normed = False,
            weights = weights_1,
            label = basename(args.relationship_file[1]))
pyplot.ylabel("Simliarity function probability", size = 20)
pyplot.xlabel("Similarity function value", size = 20)
pyplot.legend(loc='upper right')
pyplot.show()
print(abs(np.mean(sharing_0) - np.mean(sharing_1)))
print(abs(np.median(sharing_0) - np.median(sharing_1)))



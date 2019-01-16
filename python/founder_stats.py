from argparse import ArgumentParser
from collections import defaultdict
from json import dump

from population_genomes import generate_genomes
from population import PopulationUnpickler
from sex import Sex
from recomb_genome import recombinators_from_directory, RecombGenomeGenerator

def diploid_founders(diploid, founder_map):
    starts = diploid.starts
    lengths = (starts[1:] - starts[:-1]).tolist()
    for founder, length in zip(diploid.founder[:-1], lengths):
        founder_map[founder].append(length)
    last_length = diploid.end - lengths[-1]
    last_founder = diploid.founder[-1]
    founder_map[last_founder].append(last_length)
    
def simulate_founder_stats(population, genome_generator, recombinators,
                           iterations, output_file):
    node_stats = defaultdict(lambda: defaultdict(list))
    for i in range(iterations):
        genome_generator.reset()
        print("iteration {}".format(i))
        print("Cleaning genomes.")
        population.clean_genomes()
        print("Generating genomes")
        generate_genomes(population, genome_generator, recombinators, 3,
                         true_genealogy = False)
        print("Calculating founder stats")
        nodes = (node for node in population.members if node.genome is not None)
        for node in nodes:
            mother = node.genome.mother
            father = node.genome.father
            temp_founder_lengths = defaultdict(list)
            diploid_founders(mother, temp_founder_lengths)
            diploid_founders(father, temp_founder_lengths)
            founder_map = node_stats[str(node._id)]
            for founder, lengths in temp_founder_lengths.items():
                founder_map[str(founder)].append((sum(lengths), len(lengths)))


    with open(output_file, "w") as output_json:
        dump(node_stats, output_json)
            
parser = ArgumentParser(description = "Generate a classifier which can (hopefully) identify individuals in a population.")
parser.add_argument("population_file", help = "Pickled file with population")
parser.add_argument("num_iterations", type = int, default = 1000,
                    help = "Number of samples to collect from empirical distributions")
parser.add_argument("output")
parser.add_argument("--gen_back", "-g", type = int, default = 7,
                    help = "Ignore common ancestry more than the given number of generations back.")

args = parser.parse_args()

print("Loading population")
with open(args.population_file, "rb") as pickle_file:
    population = PopulationUnpickler(pickle_file).load()

print("Loading recombination data.")
recombinators = recombinators_from_directory("../data/recombination_rates/")
chrom_sizes = recombinators[Sex.Male]._num_bases
genome_generator = RecombGenomeGenerator(chrom_sizes)


print("Calculating founder stats.")

num_generations = population.num_generations
clear_index = num_generations - args.gen_back
to_clear = population.generations[clear_index].members
for node in to_clear:
    node.suspected_mother = None
    node.suspected_mother_id = None
    node.suspected_father = None
    node.suspected_father_id = None

simulate_founder_stats(population, genome_generator, recombinators,
                       args.num_iterations, args.output)

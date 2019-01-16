from argparse import ArgumentParser
from itertools import combinations
from os.path import realpath, split, join
from pickle import load

import numpy as np
from progressbar import progressbar

from shared_segment_detector import SharedSegmentDetector
from gamma import fit_hurdle_gamma
from population import PopulationUnpickler, fix_twin_parents
from cm import centimorgan_data_from_directory

parser = ArgumentParser(description = "Generate hurdle gamma params for a given population and classifier.")
parser.add_argument("population")
parser.add_argument("classifier")
parser.add_argument("--cm-ibd-threshold", type = float, default = 0.0,
                    help = "IBD segments smaller than length in cM will "
                    "go undetected")

args = parser.parse_args()

print("Loading population.", flush = True)
with open(args.population, "rb") as pickle_file:
    population = PopulationUnpickler(pickle_file).load()
fix_twin_parents(population)

print("Loading classifier", flush = True)
with open(args.classifier, "rb") as pickle_file:
    classifier = load(pickle_file)

if args.cm_ibd_threshold > 0:
    cur_path = realpath(__file__)
    parent = split(split(cur_path)[0])[0]
    rates_dir = join(parent, "data", "recombination_rates")
    print("Loading recombination data for centimorgan cutoff.", flush = True)
    recomb_data = centimorgan_data_from_directory(rates_dir)
    ibd_detector = SharedSegmentDetector(0, 5, recomb_data)
else:
    ibd_detector = SharedSegmentDetector(5000000)

labeled_nodes = classifier._labeled_nodes
labeled_node_pairs = set(combinations(labeled_nodes, 2))
related_pairs = set(classifier._distributions.keys())
cryptic_pairs = set(x for x in combinations(labeled_nodes, 2)
                    if x not in related_pairs)

print("Calculating IBD for pairs.")
lengths = []
id_map = population.id_mapping
for node_a_id, node_b_id in progressbar(cryptic_pairs):
    node_a = id_map[node_a_id]
    node_b = id_map[node_b_id]
    genome_a = node_a.genome
    genome_b = node_b.genome
    length = ibd_detector.shared_segment_length(genome_a, genome_b)
    lengths.append(length)

np_lengths = np.array(lengths, dtype = np.uint64)
print(fit_hurdle_gamma(np_lengths))

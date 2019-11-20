#!/usr/bin/env python3

from random import shuffle, getstate, setstate, seed, sample
from pickle import load, dump
from argparse import ArgumentParser
from os.path import exists, realpath, split, join

from evaluation import Evaluation
from shared_segment_detector import SharedSegmentDetector
from expansion import ExpansionData
from population import PopulationUnpickler, fix_twin_parents
from population_genomes import generate_genomes
from sex import Sex
from recomb_genome import recombinators_from_directory, RecombGenomeGenerator
from cm import centimorgan_data_from_directory
from data_logging import (write_log, change_logfile_name,
                          stop_logging, start_logging)

parser = ArgumentParser(description = "Evaluate performance of classification.")
parser.add_argument("population")
parser.add_argument("classifier")
parser.add_argument("--data-logfile",
                    help = "Filename to write log data to.")
parser.add_argument("--num_node", "-n", type = int, default = 10)
parser.add_argument("--test_node", "-t", type = int, action = "append")
parser.add_argument("--test-node-file")
parser.add_argument("--subset_labeled", "-s", type = int, default = None,
                    help = "Chose a random subset of s nodes from the set of anchor nodes. If using expansion rounds, this is the size of the initial set of anchor nodes.")
parser.add_argument("--deterministic-labeled", "-ds", action = "store_true",
                    help = "Seed the random number generator to ensure anchor node subset is deterministic.")
parser.add_argument("--rerandomize-anchors", action = "store_true",
                    default = False,
                    help = "Select a new random set of anchors for each identification.")
parser.add_argument("--anchor-node-file", help = "File specifying the ids of the anchor nodes to use.")
parser.add_argument("--ibd-threshold", type = int, default = 5000000,
                    help = "IBD segments smaller than this value will "
                    "go undetected")
parser.add_argument("--cm-ibd-threshold", type = float, default = 0.0,
                    help = "IBD segments smaller than length in cM will "
                    "go undetected")
parser.add_argument("--deterministic_random", "-d", action = "store_true",
                    help = "Seed the random number generator such that the same labeled nodes will be chosen on runs with the same number of nodes.")
parser.add_argument("--search-related", type = int, default = False,
                    help = "Search only nodes that are related to labeled nodes for which there is nonzero ibd.")
parser.add_argument("--expansion-rounds-data",
                    help = "Pickle file with data from expansion rounds.")
parser.add_argument("--smoothing-parameters",
                    help = "File with smoothing parameters for cryptic IBD. Format is a space separated list of cutoff, above_cutoff, below_cutoff, minus_eps.")
parser.add_argument("--expansion-ratio", "-r", type = float, default = 9.0,
                    help = "Confidence value required to add a node for snowball identification.")
parser.add_argument("--expansion-pool", type = int,
                    help = "Number of individuals per round in expansion rounds..")
parser.add_argument("--recombination_dir",
                    help = "Directory containing Hapmap and decode data. If this is specified, new genomes will be generated.")
parser.add_argument("--disable-probability-logging", action = "store_true", default = False, help = "This option will disable the logging of individual probabilties to the log file. Much less disk space is used when this option is specified.")
parser.add_argument("--out-of-genealogy", default = 0, type = int,
                    help = "All nodes to be identified will be erased from analyst view. This should result in always inferring the incorrect individual, as the analyst is restricted to guessing individuals in its genealogy.")

args = parser.parse_args()

if args.test_node and args.test_node_file:
    parser.error("Cannot specify both test nodes and a test node file.")

if args.anchor_node_file and args.subset_labeled:
    parser.error("Cannot specify both anchor nodes subset size and a anchor node file.")

if args.expansion_rounds_data:
    expansion_file_exists = exists(args.expansion_rounds_data)
    if not expansion_file_exists and args.subset_labeled is None:
        parser.error("A subset of labeled nodes is necessary for expansion rounds when expansion rounds data file does not already exist.")
    if expansion_file_exists:
        with open(args.expansion_rounds_data, "rb") as expansion_file:
            expansion_data = load(expansion_file)
    else:
        expansion_data = None


if args.data_logfile:
    change_logfile_name(args.data_logfile)
    start_logging()
else:
    stop_logging()

write_log("args", args)


print("Loading population.", flush = True)
with open(args.population, "rb") as pickle_file:
    population = PopulationUnpickler(pickle_file).load()
fix_twin_parents(population)

if args.recombination_dir:
    print("Generating new genomes for population.")
    print("Loading recombination rates")
    recombinators = recombinators_from_directory(args.recombination_dir)
    chrom_sizes = recombinators[Sex.Male]._num_bases
    genome_generator = RecombGenomeGenerator(chrom_sizes)
    population.clean_genomes()
    print("Generating genomes")
    generate_genomes(population, genome_generator, recombinators, 3)

print("Loading classifier", flush = True)
with open(args.classifier, "rb") as pickle_file:
    classifier = load(pickle_file)

cur_path = realpath(__file__)
parent = split(split(cur_path)[0])[0]
rates_dir = join(parent, "data", "recombination_rates")
print("Loading recombination data.", flush = True)
recomb_data = centimorgan_data_from_directory(rates_dir)
ibd_detector = SharedSegmentDetector(recomb_data, args.ibd_threshold,
                                     args.cm_ibd_threshold)    

if args.smoothing_parameters:
    print("Loading smoothing parameters from file.")
    with open(args.smoothing_parameters, "r") as params_file:
        params_lines = params_file.readlines()
    smoothing_params = [float(x) for x in params_lines[0].strip().split()]
    assert len(smoothing_params) == 4
    # Cutoff is the first param, which is an int
    smoothing_params[0] = int(smoothing_params[0])
else:
    smoothing_params = None


original_labeled = set(classifier._labeled_nodes)
if args.expansion_rounds_data and expansion_data is not None:
    print("Loading {} nodes from expansion data".format(len(expansion_data.labeled_nodes)))
    expansion_data.adjust_genomes(population)
    run_labeled = expansion_data.labeled_nodes
elif args.anchor_node_file:
    with open(args.anchor_node_file, "r") as anchor_file:
        anchor_lines = anchor_file.readlines()
    anchors = [int(anchor.strip()) for anchor in anchor_lines]
    run_labeled = anchors
elif args.subset_labeled:
    # we want the labeled nodes to be chosen randomly, but the same
    # random nodes chosen every time if the same number of labeled
    # nodes is chosen.
    sorted_labeled = list(classifier._labeled_nodes)
    sorted_labeled.sort()
    if args.deterministic_random or args.deterministic_labeled:
        rand_state = getstate()
        seed(42)
        shuffle(sorted_labeled)
        setstate(rand_state)
    else:
        shuffle(sorted_labeled)
    
    run_labeled = sorted_labeled[:args.subset_labeled]
else:
    # Use all the labeled nodes
    run_labeled = None

evaluation = Evaluation(population, classifier, run_labeled,
                        ibd_detector = ibd_detector,
                        search_related = args.search_related,
                        smoothing_parameters = smoothing_params,
                        out_of_genealogy = args.out_of_genealogy,
                        randomize_labeled = args.rerandomize_anchors)

if args.disable_probability_logging:
    evaluation.probability_logging = False

if args.expansion_rounds_data and expansion_data is None:
    expansion_data = ExpansionData(evaluation.labeled_nodes)

write_log("labeled nodes", evaluation.labeled_nodes)

id_mapping = population.id_mapping
nodes = set(member for member in population.members
             if member.genome is not None)
labeled_nodes = set(id_mapping[node_id] for node_id
                    in evaluation.labeled_nodes)
if args.test_node is not None and len(args.test_node) > 0:
    unlabeled = [id_mapping[node_id] for node_id in args.test_node]
elif args.test_node_file is not None:
    with open(args.test_node_file, "r") as test_node_file:
        node_ids = [int(node_id.strip()) for node_id
                    in test_node_file.readlines()]
    unlabeled = [id_mapping[node_id] for node_id in node_ids]
else:
    all_unlabeled = list(nodes - labeled_nodes)
    all_unlabeled.sort(key = lambda node: node._id)
    if args.deterministic_random:
        rand_state = getstate()
        seed(43)
        shuffle(all_unlabeled)
        setstate(rand_state)
    else:
        shuffle(all_unlabeled)
    unlabeled = all_unlabeled[:args.num_node]

if not args.expansion_rounds_data:
    write_log("to identify", [node._id for node in unlabeled])
    evaluation.run_evaluation(unlabeled)
    evaluation.print_metrics()
else:
    potential = set(id_mapping[node] for node in original_labeled)
    if expansion_data.remaining and len(expansion_data.remaining) > 0:
        print("Recovering identify candidates.")
        identify_candidates = set(id_mapping[node] for node
                                  in expansion_data.remaining)
    else:
        if expansion_data.original_pool:
            potential_ids = expansion_data.original_pool
            identify_candidates = set(id_mapping[node_id] for node_id
                                      in potential_ids)
            original_labeled = set(id_mapping[node_id] for node_id
                                   in expansion_data.original_labeled)
        else:
            to_sample = potential - labeled_nodes
            if args.expansion_pool:
                identify_candidates = sample(to_sample, args.expansion_pool)
            else:
                identify_candidates = to_sample
            expansion_data.original_pool = [x._id for x
                                            in identify_candidates]
            expansion_data.original_labeled = [x._id for x in labeled_nodes]
            original_labeled = labeled_nodes

    write_log("to identify", [node._id for node in identify_candidates])
    evaluation.restrict_search(potential - original_labeled)
    added = evaluation.run_expansion_round(identify_candidates,
                                           args.expansion_ratio,
                                           expansion_data,
                                           args.expansion_rounds_data)
    expansion_data.add_round()
    with open(args.expansion_rounds_data, "wb") as expansion_file:
        dump(expansion_data, expansion_file)
    write_log("expansion_data_written", {"current_node": None,
                                         "complete": True})

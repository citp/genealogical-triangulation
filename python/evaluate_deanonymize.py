#!/usr/bin/env python3

from random import shuffle, getstate, setstate, seed
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
parser.add_argument("--test_node_file")
parser.add_argument("--subset_labeled", "-s", type = int, default = None,
                    help = "Chose a random subset of s nodes from the set of labeled nodes. If using expansion rounds, this is the size of the initial set of labeled nodes.")
parser.add_argument("--ibd-threshold", type = int, default = 5000000,
                    help = "IBD segments smaller than this value will "
                    "go undetected")
parser.add_argument("--cm-ibd-threshold", type = float, default = 0.0,
                    help = "IBD segments smaller than length in cM will "
                    "go undetected")
parser.add_argument("--deterministic_random", "-d", action = "store_true",
                    help = "Seed the random number generator such that the same labeled nodes will be chosen on runs with the same number of nodes.")
parser.add_argument("--deterministic_labeled", "-ds", action = "store_true",
                    help = "Seed the random number generator to ensure labeled node subset is deterministic.")
parser.add_argument("--search-related", type = int, default = False,
                    help = "Search only nodes that are related to labeled nodes for which there is nonzero ibd.")
parser.add_argument("--expansion-rounds-data",
                    help = "Pickle file with data from expansion rounds.")
parser.add_argument("--smoothing-parameters",
                    help = "File with smoothing parameters for cryptic IBD. Format is a space separated list of cutoff, above_cutoff, below_cutoff, minus_eps.")
parser.add_argument("--expansion-ratio", "-r", type = float, default = 9.0,
                    help = "Confidence value required to add a node for snowball identification.")
parser.add_argument("--recombination_dir",
                    help = "Directory containing Hapmap and decode data. If this is specified, new genomes will be generated.")
parser.add_argument("--log-cryptic", action = "store_true", default = False,
                    help = "Log cryptic lengths in order to fit new cryptic parameters.")

args = parser.parse_args()

if args.test_node and args.test_node_file:
    parser.error("Cannot specify both test nodes and a test node file.")

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

if args.cm_ibd_threshold > 0:
    cur_path = realpath(__file__)
    parent = split(split(cur_path)[0])[0]
    rates_dir = join(parent, "data", "recombination_rates")
    print("Loading recombination data for centimorgan cutoff.", flush = True)
    recomb_data = centimorgan_data_from_directory(rates_dir)
    ibd_detector = SharedSegmentDetector(args.ibd_threshold,
                                         args.cm_ibd_threshold,
                                         recomb_data)
else:
    ibd_detector = SharedSegmentDetector(args.ibd_threshold)
    

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

evaluation = Evaluation(population, classifier,
                        ibd_detector = ibd_detector,
                        search_related = args.search_related,
                        smoothing_parameters = smoothing_params,
                        cryptic_logging = args.log_cryptic)
original_labeled = set(evaluation.labeled_nodes)
if args.expansion_rounds_data and expansion_data is not None:
    print("Loading {} nodes from expansion data".format(len(expansion_data.labeled_nodes)))
    evaluation.labeled_nodes = expansion_data.labeled_nodes
    expansion_data.adjust_genomes(population)
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
    
    evaluation.labeled_nodes = sorted_labeled[:args.subset_labeled]

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
    # potential is the set of nodes the adversary believes it has not
    # identified, based on the fact that it has not added the nodes to
    # its labeled set.
    potential = set(id_mapping[node] for node
                    in original_labeled - set(evaluation.labeled_nodes))
    # These two identify_candidate possibilities may be different if
    # we have misidentified a node. Eg identified node x as node y,
    # which means y may in fact come up later.
    if expansion_data.remaining and len(expansion_data.remaining) > 0:
        print("Recovering identify candidates.")
        identify_candidates = set(id_mapping[node] for node
                                  in expansion_data.remaining)
    else:
        identify_candidates = potential
    write_log("to identify", [node._id for node in identify_candidates])
    evaluation.restrict_search(potential)
    added = evaluation.run_expansion_round(identify_candidates,
                                           args.expansion_ratio,
                                           expansion_data,
                                           args.expansion_rounds_data)
    expansion_data.add_round(added)
    with open(args.expansion_rounds_data, "wb") as expansion_file:
        dump(expansion_data, expansion_file)
    write_log("expansion_data_written", {"current_node": None,
                                         "complete": True})

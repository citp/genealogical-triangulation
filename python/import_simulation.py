from argparse import ArgumentParser
from os.path import isfile
from pickle import dump, HIGHEST_PROTOCOL

from population import PopulationUnpickler
from classify_relationship import classifier_from_directory, classifier_from_file

parser = ArgumentParser(description = "Import simulated data.")
parser.add_argument("population_file", help = "Pickled file with population")
parser.add_argument("work_dir",
                    help = "Directory to put shared length calculations in.")
parser.add_argument("--output-pickle", default = "distributions.pickle",
                    help = "File to store distributions in. Pickle format will be used. Default is 'distributions.pickle'")
args = parser.parse_args()

print("Loading population")
with open(args.population_file, "rb") as pickle_file:
    population = PopulationUnpickler(pickle_file).load()

print("Importing simulated data")
if isfile(args.work_dir):
    classifier = classifier_from_file(args.work_dir, population.id_mapping)
else:
    classifier = classifier_from_directory(args.work_dir, population.id_mapping)

if args.output_pickle:
    print("Pickling classifier")
    with open(args.output_pickle, "wb") as pickle_file:
        dump(classifier, pickle_file, protocol = HIGHEST_PROTOCOL)

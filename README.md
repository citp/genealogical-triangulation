Getting Started
===============

This project requires both [Rust](https://www.rust-lang.org/) and [Python](https://www.python.org/). Populations are generated in Python, then exported to a json file. Rust reads the population file, then runs simulations to generate IBD data, outputting the data to a directory. The output is then read in on the python side, fit to hurdle-gamma distributions and serialized. Identification tasks are performed in Python using the hurdle-gamma distributions.

Identification will require 16+GB of RAM for populations with generations of 100k individuals.

Python Requirements
-----------

* Python 3.4+ or Python 3.3+ with [enum](https://pypi.python.org/pypi/enum34) package
* scipy/numpy
* Cython
* Statsmodels
* progressbar2



Setup
-----

* run `fetch.sh` in `data/recombination_rates` directory. This will
  fetch the appropriate recombination data from the HapMap project.
* Install the python dependencies listed above.
* run `python3 setup.py build_ext --inplace` in the predict
  directory. This will need to be run whenever any of the .pyx files
  are modified.

Running
=======

Some of the script names aren't descriptive of what is happening, as code changed over the course of the project.

Work flow
---------

The typical work flow is a three step process

1. Generate population - `python3 generate_population.py --help`
2. Export population to Rust - `python3 run_classify_relationship.py --help`
4. Run simulations in rust - `cargo run --release --bin simulate -- --help`
5. Import output from simulations - `python3 run_classify_relationship.py --help` (see the `--recover` option to recover the data from the simulation)
6. Identify - `python3 evaluate_deanonymize.py --help`

Commands
--------

### Generate population

To generate a population use the `generate_population.py` script. For
example if you cd into the `python` directory and run: `python3
generate_population.py ../data/your_tree_file ../data/recombination_rates/
--generation_size 1000 --num_generations 10 --output
population.pickle` A population with 10 generations each with 1000
members will be generated and saved to population.pickle with Python's
pickle format.


### Simulate population to generate distributions

To run experiments in rust first convert the population to a format rust can understand.: `python3
run_classify_relationship.py population.pickle work_dir 100
--num_labeled_nodes 150 --to_json file_for_rust.json`


This command will pick 150 nodes from the last generation and mark
them as "labeled" (anchor nodes in the terminology of our paper). 
It will then output a json file for the rust simulator to read.
Then run `cargo run --release --bin simulate file_for_rust.json recombination_directory output_directory`
Then it will perform 1000 experiments to sample
from the simulated empirical distributions from the (labeled,
unlabeled) pairs. A directory `output_directory` will be created with the simulation data. This tends to take on the order of days to run. If it interrupted, it can be resumed using the same `cargo` command.

To fit the hurdle gamma distributions, run `run_classify_relationship.py population.pickle output_directory 0
--recover --output_pickle distributions.pickle`. The distributions will be saved to `distributions.pickle`.

If you provide the recover option, the `--num_labeled_nodes` option will be ignored, as
the labeled nodes will be determined by `work_dir`. Recovering will
try to do `num_iteration` new iterations, on top of what may already
be in the `work_dir`. If `num_iterations` is 0, no experiments will be
run, but rather the distributions will be calculated immediately.


### Identify individuals

The final step is identifying unlabeled individuals.

Running `python3 evaluate_deanonymize.py population.pickle
distributions.pickle -n 10` will try to identify 10 random unlabeled
individuals in the population.

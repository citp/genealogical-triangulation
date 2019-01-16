extern crate genetic_privacy;
extern crate bincode;
extern crate rand;
#[macro_use]
extern crate clap;

// NOTE: The identification code in Rust is not used at this time. The python version of this code is used.

use std::fs::File;
use std::io::prelude::*;
use std::collections::{HashMap, HashSet};
use std::rc::Rc;
use std::cell::RefCell;
use std::time::Instant;

//TODO: Update random library
#[allow(deprecated)]
use rand::{thread_rng, sample};
use bincode::deserialize;
use clap::{Arg, App};

use genetic_privacy::population::node::{import_json, Node, Population};
use genetic_privacy::identify::distributions::Distribution;
use genetic_privacy::genome::recombinator::recombinators_from_directory;
use genetic_privacy::genome::genome::RecombGenomeGenerator;
use genetic_privacy::genome::population_genomes::generate_genomes;
use genetic_privacy::identify::bayes::BayesDeanonymize;

fn main() {
    let matches = App::new("Run identification")
        .about("Identify nodes.")
        .arg(Arg::with_name("POPULATION")
            .required(true))
        .arg(Arg::with_name("DISTRIBUTIONS")
            .required(true))
        .arg(Arg::with_name("RECOMBINATORS")
            .required(true))
        .arg(Arg::with_name("num-nodes")
            .takes_value(true)
            .default_value("100")
            .long("num-nodes")
            .short("n"))
        .get_matches();
    println!("Loading population");
    let imported_population = import_json(matches.value_of("POPULATION").unwrap());
    let population = imported_population.node_generator.nodes;

    let mut file = File::open(matches.value_of("DISTRIBUTIONS").unwrap()).unwrap();
    let mut distribution_data = Vec::new();
    file.read_to_end(&mut distribution_data).unwrap();
    println!("Loading distributions");
    let distribution: Distribution = deserialize(&distribution_data).unwrap();

    println!("Loading recombination data");
    let recombinators = recombinators_from_directory(matches.value_of("RECOMBINATORS").unwrap());
    let mut chrom_lengths = HashMap::new();
    for (chrom, data) in &recombinators.male.chromosomes {
        chrom_lengths.insert(*chrom, data.num_bases);
    }
    
    let genome_generator = RecombGenomeGenerator::new(&chrom_lengths);

    println!("Generating fresh genomes");
    generate_genomes(&population, &genome_generator, &recombinators, true, 3);

    
    println!("Loading bayes identifier");
    let bayes = BayesDeanonymize::new(&population, distribution);
    
    let can_identify = with_genomes(&population);
    let mut rng = thread_rng();
    let num_to_identify = value_t!(matches.value_of("num-nodes"), usize).unwrap_or_else(|e| e.exit());
    #[allow(deprecated)]
    let to_identify = sample(&mut rng, can_identify, num_to_identify);

    println!("Identifying {} nodes", to_identify.len());
    let correct = identify_individuals(&to_identify, &bayes);
    println!("{}% correctly identified",
             correct as f64 / to_identify.len() as f64);
}

fn with_genomes(population: &Population) -> Vec<Rc<RefCell<Node>>> {
    let mut ret = Vec::new();
    for node_ref in &population.members {
        let node = node_ref.borrow();
        if node.genome.is_some() {
            ret.push(node_ref.clone());
        }
    }
    ret
}

fn identify_individuals(targets: &[Rc<RefCell<Node>>],
                        bayes: &BayesDeanonymize) -> u32 {
    let mut correct = 0;
    for (i, target) in targets.iter().enumerate() {
        let target_node = target.borrow();
        println!("Identifying node {} with ID: {}", i, target_node.id);
        let target_genome = target_node.genome.as_ref().unwrap();
        let identify_time = Instant::now();
        let sibling_group = bayes.identify(target_genome, 5000000);

        let ids = group_ids(&sibling_group);
        println!("Identified group: {:?}", ids);
        if ids.contains(&target_node.id) {
            correct += 1;
            println!("Correct");
        } else {
            println!("Incorrect");
        }
        println!("It took {} seconds to identify an individual.",
                 identify_time.elapsed().as_secs());
    }
    correct
}

fn group_ids(group: &[Rc<RefCell<Node>>]) -> HashSet<u32> {
    return group.iter().map(|node| node.borrow().id).collect();
}

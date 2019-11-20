extern crate genetic_privacy;
extern crate bincode;
extern crate fnv;
extern crate clap;


// NOTE: This code isn't used, the python version is used.

use std::fs::{File, read_dir};
use std::io::prelude::*;
use std::collections::HashMap;

use bincode::{serialize, Infinite};
use fnv::FnvHashMap;
use clap::{Arg, App};

use genetic_privacy::population::node::{import_json};
use genetic_privacy::statistics::fit;
use genetic_privacy::identify::distributions::Distribution;
use genetic_privacy::genome::recombinator::recombinators_from_directory;
use genetic_privacy::genome::genome::RecombGenomeGenerator;
use genetic_privacy::genome::population_genomes::generate_genomes;
use genetic_privacy::population::node_pair::NodeIdPair;

fn main() {
    let matches = App::new("Fit distributions")
        .about("Generates classifier fit.")
        .arg(Arg::with_name("POPULATION")
            .required(true))
        .arg(Arg::with_name("RECOMBINATORS")
            .required(true))
        .arg(Arg::with_name("WORK_DIR")
            .required(true))
        .arg(Arg::with_name("OUTFILE")
            .required(true))
        .get_matches();
    let directory = matches.value_of("WORK_DIR").unwrap();
    let labeled_nodes = get_labeled_nodes(directory);
    println!("Fitting distributions");
    let params = fit::fit_directory(directory);
    let fnv_params = convert_hashmap(params);

    println!("Loading population");
    let imported_population = import_json(matches.value_of("POPULATION").unwrap());
    let population = imported_population.node_generator.nodes;

    println!("Loading recombination data");
    let recombinators = recombinators_from_directory(matches.value_of("RECOMBINATORS").unwrap());
    let mut chrom_lengths = HashMap::new();
    for (chrom, data) in &recombinators.male.chromosomes {
        chrom_lengths.insert(*chrom, data.num_bases);
    }

    let genome_generator = RecombGenomeGenerator::new(&chrom_lengths);

    println!("Generating fresh genomes");
    generate_genomes(&population, &genome_generator, &recombinators, true, 3);

    let cryptic_params = fit::fit_cryptic_labeled(&imported_population.related, &imported_population.labeled);
    let distribution = Distribution::new(fnv_params, cryptic_params, &labeled_nodes);
    println!("Serializing distributions");
    let encoded = serialize(&distribution, Infinite).unwrap();
    let mut file = File::create(matches.value_of("OUTFILE").unwrap()).unwrap();
    println!("Writing distributions to file");
    file.write_all(&encoded).unwrap();

//    println!("Writing json");
//    let distributions_vec: Vec<_> = distribution.distributions.iter().collect();
//    let json_distributions = serde_json::to_string(&distributions_vec).unwrap();
//    File::create("/home/paul/Data/genetic-privacy/distributions_2.json").unwrap().write_all(json_distributions.as_bytes()).unwrap();
}


fn get_labeled_nodes(directory: &str) -> Vec<u32> {
    let mut ret = Vec::new();
    let paths = read_dir(directory).unwrap();
    for path in paths {
        let file = path.unwrap();
        let labeled_id = file.file_name().to_str().unwrap().parse::<u32>().unwrap();
        ret.push(labeled_id);
    }
    ret
}

fn convert_hashmap<V>(mut input_map: HashMap<NodeIdPair, V>) -> FnvHashMap<NodeIdPair, V> {
    let mut map = FnvHashMap::default();
    for (key, value) in input_map.drain() {
        map.insert(key, value);
    }
    map
}

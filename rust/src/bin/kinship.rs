extern crate genetic_privacy;
extern crate clap;

use std::fs::File;
use std::io::prelude::*;

use clap::{Arg, App};

use genetic_privacy::population::node::import_json;
use genetic_privacy::population::kinship_coefficient::calculate_kinship_to_file_recursive_subset;


fn main() {
    let matches = App::new("Genetics simulations")
        .about("Calculate kinship for genealogy.")
        .arg(Arg::with_name("POPULATION")
            .required(true))
        .arg(Arg::with_name("OUT")
            .required(true))
        .arg(Arg::with_name("subset")
            .takes_value(true))
        .get_matches();
    println!("Importing population");
    let imported = import_json(matches.value_of("POPULATION").unwrap());
    let population = &imported.node_generator.nodes;

    let keep_subset: Vec<u32>;
    let subset = match matches.value_of("subset") {
        Some(filename) => {
            println!("Reading subset file.");
            let mut file = File::open(filename).unwrap();
            let mut contents = String::new();
            file.read_to_string(&mut contents).unwrap();
            keep_subset = contents.split_whitespace().map(|str_id| str_id.parse().unwrap()).collect();
            Some(keep_subset.as_slice())
        },
        None => None
    };

    println!("Calculating kinship.");
    calculate_kinship_to_file_recursive_subset(&population, subset, 3, matches.value_of("OUT").unwrap());
}

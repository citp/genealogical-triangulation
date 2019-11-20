extern crate genetic_privacy;
#[macro_use]
extern crate clap;

use std::collections::HashMap;

use clap::{Arg, App};

use genetic_privacy::population::node::import_json;
use genetic_privacy::genome::recombinator::recombinators_from_directory;
use genetic_privacy::genome::genome::RecombGenomeGenerator;
use genetic_privacy::simulate::simulate::simulate;
use genetic_privacy::genome::cm::CmConverter;
use genetic_privacy::simulate::serializer::SimulationSerializer;

fn main() {
    let matches = App::new("Genetics simulations")
        .about("Simulates population geneomes.")
        .arg(Arg::with_name("POPULATION")
            .required(true))
        .arg(Arg::with_name("RECOMBINATORS")
            .required(true))
        .arg(Arg::with_name("WORK_DIR")
            .required(true))
        .arg(Arg::with_name("iterations")
            .takes_value(true)
            .default_value("1000")
            .short("i")
            .long("iterations"))
        .arg(Arg::with_name("clobber")
            .short("c")
            .long("clobber"))
        .get_matches();
    let iterations: u32 = value_t!(matches.value_of("iterations"), u32).unwrap_or_else(|e| e.exit());
    let imported = import_json(matches.value_of("POPULATION").unwrap());
    let population = &imported.node_generator.nodes;
    //let labeled = &imported.labeled;
    let related = &imported.related;

    let recombinators = recombinators_from_directory(matches.value_of("RECOMBINATORS").unwrap());
    let cm_converter = CmConverter::new(matches.value_of("RECOMBINATORS").unwrap());
    let mut chrom_lengths = HashMap::new();
    for (chrom, data) in &recombinators.male.chromosomes {
        chrom_lengths.insert(*chrom, data.num_bases);
    }
    let genome_generator = RecombGenomeGenerator::new(&chrom_lengths);
    println!("Running simulation.");
    // simulate(population, labeled, related, &genome_generator, &recombinators,
    //          &args[3], 2000, true, 5000000);
    let clobber = matches.is_present("clobber");
    let serializer = SimulationSerializer::new(matches.value_of("WORK_DIR").unwrap(), clobber, related);
    simulate(population, related, &genome_generator, &recombinators, &cm_converter,
             serializer, iterations);
}

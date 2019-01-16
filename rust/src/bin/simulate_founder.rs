extern crate genetic_privacy;
use std::env;
use std::collections::HashMap;

use genetic_privacy::population::node::import_json;
use genetic_privacy::genome::recombinator::recombinators_from_directory;
use genetic_privacy::genome::genome::RecombGenomeGenerator;
use genetic_privacy::simulate::simulate_founder_stats::simulate_founder_stats;

fn main() {
    let args: Vec<String> = env::args().collect();
    let imported = import_json(&args[1]);
    let population = &imported.node_generator.nodes;

    let recombinators = recombinators_from_directory(&args[2]);
    let mut chrom_lengths = HashMap::new();
    for (chrom, data) in &recombinators.male.chromosomes {
        chrom_lengths.insert(*chrom, data.num_bases);
    }
    let genome_generator = RecombGenomeGenerator::new(&chrom_lengths);
    println!("Running simulation.");
    simulate_founder_stats(population, &genome_generator, &recombinators,
                           1000, &args[3]);
}

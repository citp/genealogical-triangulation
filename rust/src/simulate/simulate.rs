
use std::time::Instant;

use rayon::prelude::*;

use population::node::{Population};
use population::node_pair::RelatedPair;
use genome::genome::{RecombGenomeGenerator, Genome};
use genome::common_segments::shared_segment_length_genomes;
use genome::recombinator::{RecombinatorPair};
use genome::population_genomes::{generate_genomes, clean_genomes};
use genome::cm::CmConverter;
use simulate::serializer::SimulationSerializer;


pub fn simulate(population: &Population,
                related_pairs: &[RelatedPair],
                genome_generator: &RecombGenomeGenerator,
                recombinators: &RecombinatorPair,
                cm_converter: &CmConverter,
                mut serializer: SimulationSerializer, num_iterations: u32) {

    for i in 0..num_iterations {
        println!("On iteration {} of {}", i, num_iterations);
        genome_generator.reset();
        let start = Instant::now();
        clean_genomes(population);
        let genome_time = Instant::now();
        generate_genomes(population, genome_generator, recombinators, false, 3);
        println!("It took {} seconds to generate genomes.",
                 genome_time.elapsed().as_secs());
        let shared_time = Instant::now();
        calculate_shared_to_serializer(&related_pairs, &mut serializer, cm_converter);
        //calculate_shared_to_serializer(related_pairs, &genomes, &mut serializer, cm_converter);
        // calculate_stats_to_fds(related_pairs, &mut fds, min_segment_length);
        println!("It took {} seconds to calculate shared length.",
                 shared_time.elapsed().as_secs());
        serializer.flush();
        //flush_buffers(&mut fds);
        println!("Iteration {} took {} seconds.", i, start.elapsed().as_secs());
    }
}



fn calculate_shared_to_serializer(pairs: &[RelatedPair], //genomes: &[(&Genome, &Genome)],
                                  serializer: &mut SimulationSerializer,
                                  cm_converter: &CmConverter) {

    let borrowed_pairs: Vec<_> = pairs.iter()
        .map(|pair| (pair.unlabeled.borrow(), pair.labeled.borrow()))
        .collect();
    let genomes: Vec<(&Genome, &Genome)> = borrowed_pairs.iter()
        .map(|(a, b)| (a.genome.as_ref().unwrap().as_ref(), b.genome.as_ref().unwrap().as_ref()))
        .collect();
    let shared: Vec<_> = genomes.par_iter().map(|(a, b)| shared_segment_length_genomes(a, b, cm_converter)).collect();

    for (pair, shared_length) in pairs.iter().zip(shared) {
        let labeled = pair.labeled.borrow();
        let unlabeled = pair.unlabeled.borrow();
       // let labeled_genome = labeled.genome.as_ref().unwrap();
        //let unlabeled_genome = unlabeled.genome.as_ref().unwrap();
        //let shared_length = shared_segment_length_genomes(labeled_genome, unlabeled_genome, &cm_converter);
        serializer.insert(labeled.id, unlabeled.id, shared_length);
    }
}



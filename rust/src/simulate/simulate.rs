use std::rc::{Rc};
use std::cell::{RefCell};
use std::hash::BuildHasher;
use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::fs::{create_dir_all, remove_file, OpenOptions};
use std::path::Path;
use std::time::Instant;

use fnv::FnvHashMap;

use population::node::{Population, Node};
use population::node_pair::RelatedPair;
use genome::genome::{RecombGenomeGenerator, Genome};
use genome::common_segments::{shared_segment_length_genomes,
                              common_segment_lengths};
use genome::recombinator::{RecombinatorPair};
use genome::population_genomes::{generate_genomes, clean_genomes};
    
pub fn simulate(population: &Population, labeled_nodes: &[Rc<RefCell<Node>>],
                related_pairs: &[RelatedPair],
                genome_generator: &RecombGenomeGenerator,
                recombinators: &RecombinatorPair,
                work_dir: &str, num_iterations: u32, clobber: bool,
                min_segment_length: u32) {
    create_dir_all(work_dir).unwrap();
    let mut fds = FnvHashMap::default();
    for labeled_node in labeled_nodes {
        let node_id = labeled_node.borrow().id;
        let path = Path::new(work_dir).join(node_id.to_string());
        let file;
        if clobber {
            if path.exists() {
                remove_file(&path).unwrap();
            }
            file = OpenOptions::new().create(true).write(true)
                .open(path).unwrap();
        } else {
            file = OpenOptions::new().append(true).open(path).unwrap();
        }
        let buffer = BufWriter::new(file);
        fds.insert(node_id, buffer);
    }
    for i in 0..num_iterations {
        println!("On iteration {}", i);
        genome_generator.reset();
        let start = Instant::now();
        clean_genomes(population);
        let genome_time = Instant::now();
        generate_genomes(population, genome_generator, recombinators, false, 3);
        println!("It took {} seconds to generate genomes.",
                 genome_time.elapsed().as_secs());
        let shared_time = Instant::now();
        calculate_shared_to_fds(related_pairs, &mut fds, min_segment_length);
        // calculate_stats_to_fds(related_pairs, &mut fds, min_segment_length);
        println!("It took {} seconds to calculate shared length.",
                 shared_time.elapsed().as_secs());
        flush_buffers(&mut fds);
        println!("Iteration {} took {} seconds.", i, start.elapsed().as_secs());
    }
}

fn flush_buffers<T: Write, U: BuildHasher>(bufs: &mut HashMap<u32, BufWriter<T>, U>) {
    for buf in bufs.values_mut() {
        buf.flush().unwrap();
    }
}

fn shared_segment_avglength_count(a: &Genome, b: &Genome, min_length: u32) -> (f64, u32) {
    let lengths: Vec<_> = common_segment_lengths(a, b).iter()
        .filter(|&length| *length >= min_length).cloned().collect();
    if lengths.len() == 0 {
        return (0.0, 0);
    }
    let sum: f64 = lengths.iter().map(|&x| x as f64).sum();
    let average_length = sum / (lengths.len() as f64);
    (average_length, lengths.len() as u32)
}

fn calculate_shared_to_fds<T: Write, U: BuildHasher>(pairs: &[RelatedPair],
                                     fds: &mut HashMap<u32, BufWriter<T>, U>,
                                     min_segment_length: u32) {
    for pair in pairs {
        let labeled = pair.labeled.borrow();
        let unlabeled = pair.unlabeled.borrow();
        let labeled_genome = labeled.genome.as_ref().unwrap();
        let unlabeled_genome = unlabeled.genome.as_ref().unwrap();
        let shared_length = shared_segment_length_genomes(labeled_genome,
                                                          unlabeled_genome,
                                                          min_segment_length);
        let fd = fds.get_mut(&labeled.id).unwrap();
        write!(fd, "{}\t{}\n", unlabeled.id, shared_length).unwrap();
    }
}

fn calculate_stats_to_fds<T: Write>(pairs: &[RelatedPair],
                                     fds: &mut HashMap<u32, BufWriter<T>>,
                                     min_segment_length: u32) {
    for pair in pairs {
        let labeled = pair.labeled.borrow();
        let unlabeled = pair.unlabeled.borrow();
        let labeled_genome = labeled.genome.as_ref().unwrap();
        let unlabeled_genome = unlabeled.genome.as_ref().unwrap();
        let avg_count = shared_segment_avglength_count(labeled_genome,
                                                       unlabeled_genome,
                                                       min_segment_length);
        let avg = avg_count.0;
        let count = avg_count.1;
        let fd = fds.get_mut(&labeled.id).unwrap();
        write!(fd, "{}\t{}\t{}\n", unlabeled.id, avg, count).unwrap();
    }
}




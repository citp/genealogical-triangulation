use std::rc::{Rc};
use std::cell::{RefCell};
use std::collections::HashMap;
use std::time::Instant;
use std::io::prelude::*;
use std::fs::File;

use serde_json;

use population::node::{Population, Node};
use genome::genome::{RecombGenomeGenerator, Genome};
use genome::diploid::Diploid;
use genome::recombinator::{RecombinatorPair};
use genome::population_genomes::{generate_genomes, clean_genomes};

#[derive(Clone, Debug, Copy, Serialize, Deserialize)]
pub struct IterationStats {
    pub count: u32,
    pub total_length: u64,
}

impl IterationStats {
    pub fn new(count: u32, length: u64) -> IterationStats {
        IterationStats {count: count, total_length: length}
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FounderStats {
    pub founder_stats: HashMap<u32, Vec<IterationStats>>        
}

impl FounderStats {
    pub fn new() -> FounderStats {
        FounderStats {founder_stats: HashMap::new()}
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NodeStats {
    pub node_map: HashMap<u32, FounderStats>
}

impl NodeStats {
    pub fn new() -> NodeStats {
        NodeStats {node_map: HashMap::new()}
    }

    pub fn update(&mut self, node_id: u32, genome: &Genome) {
        let mut temp_founder_lengths = HashMap::new();
        self.update_diploid(&genome.mother, &mut temp_founder_lengths);
        self.update_diploid(&genome.father, &mut temp_founder_lengths);
        for (founder_id, lengths) in &temp_founder_lengths {
            let total_length: u64 = lengths.iter().map(|x| *x as u64).sum();
            let count = lengths.len() as u32;
            self.add_founder_count_length(node_id, *founder_id,
                                          count, total_length);
        }
    }

    fn add_founder_count_length(&mut self, node: u32, founder: u32,
                           count: u32, length: u64) {
        let entry = self.node_map.entry(node)
            .or_insert_with(FounderStats::new);
        let lengths_vec = entry.founder_stats.entry(founder)
            .or_insert_with(Vec::new);
        lengths_vec.push(IterationStats::new(count, length));
    }

    fn update_diploid(&mut self, diploid: &Diploid,
                      temp_founder_lengths: &mut HashMap<u32, Vec<u32>>) {
        let starts = &diploid.starts;
        let founder = &diploid.founder;
        for i in 0..starts.len() - 1 {
            let start = starts[i];
            let stop = starts[i + 1];
            let founder = founder[i];
            let len = stop - start;
            let lengths = temp_founder_lengths.entry(founder)
                .or_insert_with(Vec::new);
            lengths.push(len);
            // self.add_founder_segment(node_id, founder, len);
        }
        let last_len = diploid.end - starts.last().unwrap();
        let last_founder = founder.last().unwrap();
        let lengths = temp_founder_lengths.entry(*last_founder)
                .or_insert_with(Vec::new);
        lengths.push(last_len);
        
    }
}

fn calculate_founder_stats(nodes: &[Rc<RefCell<Node>>],
                           stats: &mut NodeStats) {
    for node in nodes {
        let borrow_node = node.borrow();
        if let Some(x) = borrow_node.genome.as_ref() {
            stats.update(borrow_node.id, x);
        }
    }
}


pub fn simulate_founder_stats(population: &Population,
                              genome_generator: &RecombGenomeGenerator,
                              recombinators: &RecombinatorPair,
                              num_iterations: u32,
                              output_file: &str) {
    let mut node_stats = NodeStats::new();
    for i in 0..num_iterations {
        println!("On iteration {}", i);
        genome_generator.reset();
        let start = Instant::now();
        clean_genomes(population);
        let genome_time = Instant::now();
        generate_genomes(population, genome_generator, recombinators, false, 3);
        println!("It took {} seconds to generate genomes.",
                 genome_time.elapsed().as_secs());

        calculate_founder_stats(&population.members, &mut node_stats);
        let shared_time = Instant::now();
        println!("It took {} seconds to calculate founder stats.",
                 shared_time.elapsed().as_secs());
        println!("Iteration {} took {} seconds.", i, start.elapsed().as_secs());
    }
    let mut lens = Vec::new();
    for (_, founder_stats) in &node_stats.node_map {
        for (_, stats_vec) in &founder_stats.founder_stats {
            lens.push(stats_vec.len() as u32);
        }
    }

    println!("Max stat entries {}", lens.iter().max().unwrap());
    let json = serde_json::to_string(&node_stats).unwrap();
    let mut file = File::create(output_file).unwrap();
    file.write_all(json.as_bytes()).unwrap();
    file.flush().unwrap();
}

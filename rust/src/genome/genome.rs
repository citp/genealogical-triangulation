use std::collections::HashMap;
use std::cell::Cell;

use genome::diploid::Diploid;
use genome::recombinator::cum_sum;

pub const CHROMOSOMES: &'static [u32] = &[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22];
// const NUM_CHROMOSOMES: u32 = 22;

#[derive(Clone)]
pub struct Genome {
    pub mother: Diploid,
    pub father: Diploid
}

pub struct RecombGenomeGenerator {
    pub chromosome_lengths: HashMap<u32, u32>,
    pub total_length: u32,
    pub chrom_start_offset: HashMap<u32, u32>,
    genome_id: Cell<u32>
}

impl RecombGenomeGenerator {
    pub fn generate(&self) -> Genome {
        let mut starts: Vec<u32> = vec![];
        for chrom in CHROMOSOMES {
            let offset = self.chrom_start_offset.get(chrom).unwrap();
            starts.push(*offset);
        }
        let mother_founder = vec![self.genome_id.get(); CHROMOSOMES.len()];
        self.genome_id.set(self.genome_id.get() + 1);
        let father_founder = vec![self.genome_id.get(); CHROMOSOMES.len()];
        self.genome_id.set(self.genome_id.get() + 1);
        let mother = Diploid{starts: starts.clone(), founder: mother_founder,
                              end: self.total_length };
        let father = Diploid{starts: starts, founder: father_founder,
                             end: self.total_length };
        Genome {mother: mother, father: father}
    }

    pub fn new(chromosome_lengths: &HashMap<u32, u32>) 
               -> RecombGenomeGenerator {
        let lengths_copy = chromosome_lengths.clone();
        let total_length = lengths_copy.values().fold(0, |acc, &x| acc + x);
        let lengths_vec: Vec<u32> = CHROMOSOMES.iter()
            .map(|chrom| lengths_copy[chrom]).collect();
        let cum_lengths = cum_sum(&lengths_vec);
        let mut start_offsets: HashMap<u32, u32> = HashMap::new();
        start_offsets.insert(1, 0);
        for (i, cum_length) in cum_lengths.iter().enumerate() {
            start_offsets.insert((i + 2) as u32, *cum_length);
        }

        RecombGenomeGenerator { chromosome_lengths: lengths_copy,
                                total_length: total_length,
                                chrom_start_offset: start_offsets,
                                genome_id: Cell::new(0) }
    }

    pub fn reset(&self) {
        self.genome_id.set(0);
    }
}



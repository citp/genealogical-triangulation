use std::rc::Rc;

use rand::{Rng, thread_rng};
// use rand::distributions::{IndependentSample, Range};

use population::node::Population;
use genome::genome::{Genome, CHROMOSOMES, RecombGenomeGenerator};
use genome::recombinator::{Recombinator, RecombinatorPair};
use genome::diploid::Diploid;

pub fn generate_genomes(population: &Population,
                        generator: &RecombGenomeGenerator,
                        recombinators: &RecombinatorPair,
                        true_genealogy: bool,
                        keep_last: u32) {
    let generations = &population.generations;
    for (generation_num, generation) in generations.iter().enumerate() {
        for node_ref in &generation.members {
            let mut node = node_ref.borrow_mut();
            if node.genome.is_some() {
                continue;
            }
            let mother_option;
            let father_option;
            if true_genealogy {
                mother_option = node.mother.upgrade();
                father_option = node.father.upgrade();
            } else {
                mother_option = node.suspected_mother.upgrade();
                father_option = node.suspected_father.upgrade();
            }
            if mother_option.is_none() && father_option.is_none() {
                node.genome = Some(Rc::new(generator.generate()));
                continue;
            }
            if let Some(twin) = node.twin.upgrade() {
                let twin_borrow = twin.borrow();
                let twin_mother_option;
                let twin_father_option;
                if true_genealogy {
                    twin_mother_option = twin_borrow.mother.upgrade();
                    twin_father_option = twin_borrow.father.upgrade();
                } else {
                    twin_mother_option = twin_borrow.suspected_mother.upgrade();
                    twin_father_option = twin_borrow.suspected_father.upgrade();
                }
                if twin_borrow.genome.is_some() && mother_option == twin_mother_option && father_option == twin_father_option {
                    node.genome = twin_borrow.genome.clone();
                    continue;
                }
            }
            let mother_genome;
            if mother_option.is_none() {
                mother_genome = Rc::new(generator.generate());
            } else {
                let mother = mother_option.unwrap();
                let mother_ref = mother.borrow();
                mother_genome = mother_ref.genome.as_ref().unwrap().clone();
            }
            let father_genome;
            if father_option.is_none() {
                father_genome = Rc::new(generator.generate());
            } else {
                let father = father_option.unwrap();
                let father_ref = father.borrow();
                father_genome = father_ref.genome.as_ref().unwrap().clone();
            }
            node.genome = Some(Rc::new(mate(mother_genome.as_ref(),
                                            father_genome.as_ref(),
                                            &recombinators)));
        }
        if keep_last <= generation_num as u32 {
            let to_clean_gen_num = generation_num - keep_last as usize;
            let to_clean = &population.generations[to_clean_gen_num];
            for node_ref in &to_clean.members {
                let mut node = node_ref.borrow_mut();
                node.genome = None;
            }
        }
    }
    
}

pub fn clean_genomes(population: &Population) {
    for node_ref in &population.members {
        let mut node = node_ref.borrow_mut();
        node.genome = None;
    }
}

fn pick_chroms_for_diploid(genome: &Genome, recombinator: &Recombinator)
                           -> Diploid {
    let recomb_genome = recombinator.recombination(genome);
    let mut starts = Vec::new();
    let mut founder = Vec::new();
    let mut rng = thread_rng();
    for &chrom in CHROMOSOMES {
        let tmp_diploid;
        if rng.gen() {
            tmp_diploid = &recomb_genome.mother;
        } else {
            tmp_diploid = &recomb_genome.father;
        }
        let start = recombinator.chromosomes[&chrom].chrom_start_offset;
        let start_i = match tmp_diploid.starts.binary_search(&start) {
            Ok(i) => i,
            Err(i) => i
        };
        let stop_i;
        if chrom != CHROMOSOMES[CHROMOSOMES.len() - 1] {
            let stop = recombinator.chromosomes[&(chrom + 1)].chrom_start_offset;
            stop_i = match tmp_diploid.starts.binary_search(&stop) {
                Ok(i) => i,
                Err(i) => i
            };
        } else {
            stop_i = tmp_diploid.starts.len();
        }
        starts.extend_from_slice(&tmp_diploid.starts[start_i..stop_i]);
        founder.extend_from_slice(&tmp_diploid.founder[start_i..stop_i]);
    }
    Diploid {starts: starts, founder: founder, end: recomb_genome.mother.end}
}

fn mate(mother: &Genome, father: &Genome,
        recombinators: &RecombinatorPair) -> Genome {
    // mother_recombinator: &Recombinator,
    // father_recombinator: &Recombinator) -> Genome {
    let from_mother = pick_chroms_for_diploid(mother, &recombinators.female);
    let from_father = pick_chroms_for_diploid(father, &recombinators.male);
    Genome {mother: from_mother, father: from_father}
}

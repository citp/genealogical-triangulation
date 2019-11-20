use std::collections::HashMap;
use std::rc::Rc;
use std::cell::RefCell;
use std::f64::NEG_INFINITY;
use std::ops::Deref;

use rgsl::error::set_error_handler;
use rgsl::Value;

use genome::genome::Genome;
use population::node::{Population, Node};
use population::node_pair::NodeIdPair;
use identify::distributions::{Distribution, DistributionAlt, get_shape, get_scale, get_length};
use genome::common_segments::shared_segment_length_genomes;
use genome::cm::CmConverter;

pub struct BayesDeanonymize<'a> {
    population: &'a Population,
    distribution: Distribution,
    alt_distribution: DistributionAlt,
    id_map: HashMap<u32, Rc<RefCell<Node>>>,
    cm_converter: CmConverter
}

fn error_handling(error_str: &str, file: &str, line: u32, error_value: Value) {
    println!("[{:?}] '{}:{}': {}", error_value, file, line, error_str);
    println!("Shape {}, scale {}, length {}", get_shape(), get_scale(), get_length());
    panic!("Error in gsl encountered");
}

#[derive(Clone, Copy)]
struct Shared(u32, f64);

impl<'a> BayesDeanonymize<'a> {
    pub fn new (population: &'a Population,
            distribution: Distribution, cm_converter: CmConverter) -> BayesDeanonymize {
        let mut id_map = HashMap::new();
        for node_ref in &population.members {
            let node = node_ref.borrow();
            id_map.insert(node.id, node_ref.clone());
        }
//        println!("Distributions has {} entries",
//                 distribution.distributions.len());
        println!("Converting distributions");
        let alt = (&distribution).into();
        BayesDeanonymize {population: population, distribution: distribution, alt_distribution: alt,
                          id_map: id_map, cm_converter: cm_converter}
    }
    pub fn identify(&self, genome: &Genome) -> Vec<Rc<RefCell<Node>>> {
        let num_labeled = self.distribution.labeled_nodes.len();
        let mut shared_list = Vec::with_capacity(num_labeled);
        for labeled_node_id in &self.distribution.labeled_nodes {
            let labeled = self.id_map[labeled_node_id].borrow();
            let labeled_genome = labeled.genome.as_ref().unwrap();
            let shared = shared_segment_length_genomes(genome, labeled_genome, &self.cm_converter);
            shared_list.push(Shared(*labeled_node_id, shared));
        }

        let _old_handler = set_error_handler(Some(error_handling));

        //let mut node_probabilities = FnvHashMap::with_capacity_and_hasher(self.population.members.len(), Default::default());
        let mut temp_probs = Vec::with_capacity(shared_list.len());
        let mut node_probabilities = Vec::with_capacity(self.population.members.len());
        for node_ref in &self.population.members {
            let node = node_ref.borrow();
            //let mut prob_accum: f64 = 0.0;
            //probs.truncate(0);
            for &Shared(labeled_node_id, shared) in &shared_list {
                let pair = NodeIdPair { unlabeled: node.id,
                                        labeled: labeled_node_id};
                //prob_accum += self.distribution.get_probability(shared, &pair).ln();
                //temp_probs.push(self.distribution.get_probability(shared, &pair));
                temp_probs.push(self.distribution.get_log_probability(shared, &pair));
            }
            //let prob_accum: f64 = temp_probs.drain(..).map(|x| x.ln()).sum();
            let prob_accum: f64 = temp_probs.drain(..).sum();
            //node_probabilities.insert(node.id, prob_accum);
            node_probabilities.push((node.id, prob_accum));
        }

        let mut max_pair = (0, NEG_INFINITY);
        for (node_id, log_prob) in node_probabilities.iter() {
            if max_pair.1 < *log_prob {
                max_pair = (*node_id, *log_prob)
            }
        }
        let max_node_rc = &self.id_map[&max_pair.0];
        println!("Identified {}", max_pair.0);
        let sibling_group = get_sibling_group(max_node_rc);
        sibling_group
    }

    pub fn identify_alt(&self, genome: &Genome) -> Vec<Rc<RefCell<Node>>> {
        let num_labeled = self.distribution.labeled_nodes.len();
        let mut shared_list = Vec::with_capacity(num_labeled);
        for labeled_node_id in &self.distribution.labeled_nodes {
            let labeled = self.id_map[labeled_node_id].borrow();
            let labeled_genome = labeled.genome.as_ref().unwrap();
            let shared = shared_segment_length_genomes(genome, labeled_genome, &self.cm_converter);
            shared_list.push((*labeled_node_id, shared));
        }

        shared_list.sort_unstable_by_key(|x| x.0);

        let _old_handler = set_error_handler(Some(error_handling));

        let mut node_probabilities = Vec::with_capacity(self.population.members.len());
        for node_ref in &self.population.members {
            let node = node_ref.borrow();
            let mut temp_probs = self.alt_distribution.get_log_probabilities(node.id, &shared_list);
            let prob_accum: f64 = temp_probs.drain(..).sum();
            node_probabilities.push((node.id, prob_accum));
        }

        let mut max_pair = (0, NEG_INFINITY);
        for (node_id, log_prob) in node_probabilities.iter() {
            if max_pair.1 < *log_prob {
                max_pair = (*node_id, *log_prob)
            }
        }
        let max_node_rc = &self.id_map[&max_pair.0];
        println!("Identified {}", max_pair.0);
        let sibling_group = get_sibling_group(max_node_rc);
        sibling_group
    }
}

fn get_sibling_group(node: &Rc<RefCell<Node>>) -> Vec<Rc<RefCell<Node>>> {
    let mut ret = Vec::new();
    let borrowed = node.borrow();
    let n: &Node = borrowed.deref();
    let father_ref = match n.father.upgrade() {
        Some(node_ref) => node_ref,
        None => {
            ret.push(node.clone());
            return ret;
        }
    };
    let mother_ref = match n.mother.upgrade() {
        Some(node_ref) => node_ref,
        None => {
            ret.push(node.clone());
            return ret;
        }
    };
    for father_child in &father_ref.borrow().children {
        for mother_child in &mother_ref.borrow().children {
            if father_child.borrow().id == mother_child.borrow().id {
                ret.push(father_child.clone());
            }
        }
    }
    assert!(ret.len() > 0);
    ret
}

use std::collections::HashMap;
use std::rc::Rc;
use std::cell::RefCell;
use std::f64::NEG_INFINITY;
use std::ops::Deref;

use genome::genome::Genome;
use population::node::{Population, Node};
use population::node_pair::NodeIdPair;
use identify::distributions::Distribution;
use genome::common_segments::shared_segment_length_genomes;

pub struct BayesDeanonymize<'a> {
    population: &'a Population,
    distribution: Distribution,
    id_map: HashMap<u32, Rc<RefCell<Node>>>
}

#[derive(Clone, Copy)]
struct Shared(u32, u64);

impl<'a> BayesDeanonymize<'a> {
    pub fn new (population: &'a Population,
            distribution: Distribution) -> BayesDeanonymize {
        let mut id_map = HashMap::new();
        for node_ref in &population.members {
            let node = node_ref.borrow();
            id_map.insert(node.id, node_ref.clone());
        }
        println!("Distributions has {} entries",
                 distribution.distributions.len());
        BayesDeanonymize {population: population, distribution: distribution,
                          id_map: id_map}
    }
    pub fn identify(&self, genome: &Genome, min_length: u32) -> Vec<Rc<RefCell<Node>>> {
        let num_labeled = self.distribution.labeled_nodes.len();
        let mut shared_list = Vec::with_capacity(num_labeled);
        for labeled_node_id in &self.distribution.labeled_nodes {
            let labeled = self.id_map[labeled_node_id].borrow();
            let labeled_genome = labeled.genome.as_ref().unwrap();
            let shared = shared_segment_length_genomes(genome, labeled_genome,
                                                       min_length);
            shared_list.push(Shared(*labeled_node_id, shared));
        }

        let mut node_probabilities: HashMap<u32, f64> = HashMap::new();
        for node_ref in &self.population.members {
            let node = node_ref.borrow();
            let mut prob_accum: f64 = 0.0;
            for &Shared(labeled_node_id, shared) in &shared_list {
                let pair = NodeIdPair { labeled: node.id,
                                        unlabeled: labeled_node_id};
                prob_accum += self.distribution.get_probability(shared,
                                                                &pair).ln();
            }
            node_probabilities.insert(node.id, prob_accum);
        }

        let mut max_pair = (0, NEG_INFINITY);
        for (&node_id, &log_prob) in node_probabilities.iter() {
            if max_pair.1 < log_prob {
                max_pair = (node_id, log_prob)
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

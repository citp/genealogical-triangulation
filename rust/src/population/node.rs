use std::rc::{Weak, Rc};
use std::cell::{Cell, RefCell};
use std::hash::{Hasher, Hash};
// use std::collections::HashMap;
use std::io::prelude::*;
use std::fs::File;
use std::cmp;

extern crate serde_json;

use genome::genome::Genome;
use population::node_pair::RelatedPair;

#[derive(Serialize, Deserialize, Clone, Copy, Debug, Hash, Eq, PartialEq)]
pub enum Sex {
    Male,
    Female
}

pub struct Node {
    pub id: u32,
    pub sex: Sex,
    pub generation: u32,
    pub father: Weak<RefCell<Node>>,
    pub mother: Weak<RefCell<Node>>,
    pub suspected_father: Weak<RefCell<Node>>,
    pub suspected_mother: Weak<RefCell<Node>>,
    pub children: Vec<Rc<RefCell<Node>>>,
    pub suspected_children: Vec<Rc<RefCell<Node>>>,
    pub twin: Weak<RefCell<Node>>,
    pub genome: Option<Rc<Genome>>
}

impl Hash for Node {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.id.hash(state);
        self.sex.hash(state);
        self.generation.hash(state);
    }
}

impl PartialEq for Node {
    fn eq(&self, other: &Node) -> bool {
        self.id == other.id && self.sex == other.sex && self.generation == other.generation
    }
}

impl Eq for Node {}

pub struct Generation {
    pub members: Vec<Rc<RefCell<Node>>>
}

impl Generation {
    fn new() -> Generation {
        Generation { members: Vec::new() }
    }
}

pub struct Population {
    pub members: Vec<Rc<RefCell<Node>>>,
    pub generations: Vec<Generation>
}

impl Population {
    fn new(members: Vec<Rc<RefCell<Node>>>) -> Population {
        let max_generation = members.iter()
            .map(|node| node.borrow().generation)
            .fold(0, cmp::max) + 1;
        let mut generations = Vec::new();
        for _ in 0..max_generation {
            generations.push(Generation::new());
        }
        for node in &members {
            let node_generation = node.borrow().generation;
            let mut generation = &mut generations[node_generation as usize];
            generation.members.push(node.clone());
        }
        Population { members: members,
                     generations: generations }
    }
}

pub struct NodeGenerator {
    id: Cell<u32>,
    pub nodes: Population
    // id_mapping: HashMap<u32, Rc<Node>>
}

impl Node {
    fn twin_node(&self) -> Node {
        Node { id: self.id, father: self.father.clone(), sex: self.sex,
               mother: self.mother.clone(), generation: self.generation,
               suspected_father: self.suspected_father.clone(),
               suspected_mother: self.suspected_mother.clone(),
               children: Vec::new(), suspected_children: Vec::new(),
               // TODO: The twin property here is wrong.
               twin: self.twin.clone(), genome: self.genome.clone() }
    }
}

#[derive(Serialize, Deserialize)]
pub struct JsonNode {
    pub id: u32,
    pub sex: Sex,
    pub generation: u32,
    pub father: Option<u32>,
    pub mother: Option<u32>,
    pub suspected_father: Option<u32>,
    pub suspected_mother: Option<u32>,
    pub twin: Option<u32>
}

#[derive(Serialize, Deserialize)]
pub struct JsonRelated {
    labeled_node: u32,
    unlabeled_node: u32,
}

#[derive(Serialize, Deserialize)]
pub struct JsonPopulation {
    pub nodes: Vec<JsonNode>,
    pub related: Vec<JsonRelated>,
    pub labeled: Vec<u32>
}

pub struct ImportedPopulation {
    pub node_generator: NodeGenerator,
    pub related: Vec<RelatedPair>,
    pub labeled: Vec<Rc<RefCell<Node>>>
}

fn import_node(json_node: &JsonNode, nodes: &[Rc<RefCell<Node>>]) -> Node {
    Node {
        id: json_node.id,
        sex: json_node.sex,
        generation: json_node.generation,
        father: match json_node.father {
            Some(father_id) => Rc::downgrade(&nodes[father_id as usize]),
            None => Weak::new()
        },
        mother: match json_node.mother {
            Some(mother_id) => Rc::downgrade(&nodes[mother_id as usize]),
            None => Weak::new()
        },
        suspected_father: match json_node.suspected_father {
            Some(suspected_father_id) => Rc::downgrade(&nodes[suspected_father_id as usize]),
            None => Weak::new()
        },
        suspected_mother: match json_node.suspected_mother {
            Some(suspected_mother_id) => Rc::downgrade(&nodes[suspected_mother_id as usize]),
            None => Weak::new()
        },
        children: Vec::new(),
        suspected_children: Vec::new(),
        twin: Weak::new(),
        genome: None
    }
}

pub fn import_json(json_file: &str) -> ImportedPopulation {
    let mut file = File::open(json_file).unwrap();
    let mut contents = String::new();
    file.read_to_string(&mut contents).unwrap();
    let mut json_population: JsonPopulation = serde_json::from_str(&contents).unwrap();
    let num_nodes = json_population.nodes.len();
    json_population.nodes.sort_by_key(|n| n.id);
    let mut nodes = Vec::with_capacity(num_nodes);
    for (i, json_node) in json_population.nodes.iter().enumerate() {
        assert!(i as u32 == json_node.id);
        let new_node = import_node(&json_node, &nodes);

        let node_rc = Rc::new(RefCell::new(new_node));

        if let Some(father_id) = json_node.father {
            let mut father = nodes[father_id as usize].borrow_mut();
            father.children.push(node_rc.clone());
        }
        if let Some(mother_id) = json_node.mother {
            let mut mother = nodes[mother_id as usize].borrow_mut();
            mother.children.push(node_rc.clone());
        }
        if let Some(father_id) = json_node.suspected_father {
            let mut suspected_father = nodes[father_id as usize].borrow_mut();
            suspected_father.suspected_children.push(node_rc.clone());
        }
        if let Some(mother_id) = json_node.suspected_mother {
            let mut suspected_mother = nodes[mother_id as usize].borrow_mut();
            suspected_mother.suspected_children.push(node_rc.clone());
        }

        nodes.push(node_rc.clone());
    }

    for json_node in &json_population.nodes {
        if let Some(twin_id) = json_node.twin {
            let twin = Rc::downgrade(&nodes[twin_id as usize]);
            let node_rc = &nodes[json_node.id as usize];
            let mut actual_node = node_rc.borrow_mut();
            actual_node.twin = twin;
        }
    }

    let labeled_nodes = json_population.labeled.iter()
        .map(|&node_id| nodes[node_id as usize].clone())
        .collect();

    let mut related_pairs = Vec::with_capacity(json_population.related.len());
    for pair in json_population.related {
        let labeled_node = nodes[pair.labeled_node as usize].clone();
        let unlabeled_node = nodes[pair.unlabeled_node as usize].clone();
        related_pairs.push(RelatedPair {labeled: labeled_node,
                                        unlabeled: unlabeled_node});
    }

    let node_generator = NodeGenerator {id: Cell::new(num_nodes as u32),
                                        nodes: Population::new(nodes) };
    ImportedPopulation {node_generator: node_generator,
                        related: related_pairs,
                        labeled: labeled_nodes}
}

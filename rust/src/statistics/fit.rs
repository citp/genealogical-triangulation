use std::fs::{File, read_dir};
use std::path::Path;
use std::io::BufRead;
use std::io::BufReader;
use std::collections::HashMap;
use std::rc::Rc;
use std::cell::RefCell;
use std::io::Read;

use serde_json::{Result, Value};

use population::node_pair::NodeIdPair;

use super::gamma::{HurdleGammaParams, fit_hurdle_gamma};
use population::node::Node;
use population::node_pair::RelatedPair;
//use genome::common_segments::shared_segment_length_genomes;

#[derive(Clone, Copy, Debug)]
pub struct NodeParams(u32, HurdleGammaParams);

pub fn fit_file<T: AsRef<Path>>(filename: T) -> Vec<NodeParams> {
    let file = File::open(filename).unwrap();
    let reader = BufReader::new(file);
    let mut data: HashMap<u32, Vec<u64>> = HashMap::new();
    for line in reader.lines() {
        let unwrapped_line = line.unwrap();
        let parts : Vec<_> = unwrapped_line.split_whitespace().collect();
        if parts.len() != 2 {
            continue;
        }
        let node_id = match parts[0].parse::<u32>() {
            Ok(x) => x,
            Err(_) => {
                println!("Error parsing node value {}", parts[0]);
                continue;
            }
                
        };
        let shared = match parts[1].parse::<u64>() {
            Ok(x) => x,
            Err(_) => {
                println!("Error parsing IBD value {}", parts[1]);
                continue;
            }
        };
        let mut current_vec = data.entry(node_id).or_insert_with(Vec::new);
        current_vec.push(shared);
    }
    let mut results = Vec::with_capacity(data.len());
    for (node_id, data_vec) in data.iter() {
        if let Some(params) = fit_hurdle_gamma(data_vec.as_slice()) {
            assert!(!params.has_nan_parameters());
            results.push(NodeParams(*node_id, params));
        }
    }
    results
}

pub fn fit_directory<T: AsRef<Path>> (directory: T)
                                      -> HashMap<NodeIdPair,
                                                 HurdleGammaParams> {
    let paths = read_dir(directory).unwrap();
    let mut result = HashMap::new();
    for path in paths {
        let file = path.unwrap();
        if !file.metadata().unwrap().is_file() {
            continue;
        }
        let labeled_id = file.file_name().to_str().unwrap().parse::<u32>().unwrap();
        let file_params = fit_file(file.path());
        for distribution in file_params {
            let unlabeled_id = distribution.0;
            let params = distribution.1;
            let pair = NodeIdPair { labeled: labeled_id,
                                    unlabeled: unlabeled_id };
            result.insert(pair, params);
        }
        
    }
    result
}

pub fn fit_cryptic_labeled(_related: &[RelatedPair], _labeled: &[Rc<RefCell<Node>>]) -> HurdleGammaParams {
    //TODO Finish this
//    let labeled_refs: Vec<_> = labeled.iter().map(|x| x.borrow()).collect();
//    for (i, labeled_a) in labeled_refs.iter().enumerate() {
//        for labeled_b in labeled_refs[i..].iter() {
//
//        }
//
//    }
    HurdleGammaParams {shape: 1.1573974490526806, scale: 12642827.473324005, zero_prob: 0.9876864782229996}
}


pub fn load_python_json() -> Result<HashMap<NodeIdPair, HurdleGammaParams>> {
    let mut data = String::new();
    File::open("/home/paul/Data/genetic-privacy/python_distributions.json").unwrap().read_to_string(&mut data).unwrap();
    let v: Value = serde_json::from_str(&data)?;
    let arr = v.as_array().unwrap();
    let mut ret = HashMap::new();
    for value in arr {
        let nodes = value[0].as_array().unwrap();
        let unlabeled = nodes[0].as_u64().unwrap() as u32;
        let labeled = nodes[1].as_u64().unwrap() as u32;

        let params = value[1].as_array().unwrap();
        let shape = params[0].as_f64().unwrap();
        let scale = params[1].as_f64().unwrap();
        let zero_prob = params[2].as_f64().unwrap();

        let key = NodeIdPair {labeled, unlabeled};
        let params = HurdleGammaParams {shape, scale, zero_prob };
        ret.insert(key, params);
    }
    Ok(ret)
}
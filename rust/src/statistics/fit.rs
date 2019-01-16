use std::fs::{File, read_dir};
use std::path::Path;
use std::io::BufRead;
use std::io::BufReader;
use std::collections::HashMap;

use population::node_pair::NodeIdPair;

use super::gamma::{HurdleGammaParams, fit_hurdle_gamma};

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
        let params = fit_hurdle_gamma(data_vec.as_slice());
        results.push(NodeParams(*node_id, params));
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

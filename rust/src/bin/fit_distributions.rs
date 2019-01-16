extern crate genetic_privacy;
extern crate bincode;
extern crate fnv;

// NOTE: This code isn't used, the python version is used.

use std::fs::{File, read_dir};
use std::io::prelude::*;
use std::env;
use std::collections::HashMap;

use bincode::{serialize, Infinite};
use fnv::FnvHashMap;

use genetic_privacy::statistics::fit;
use genetic_privacy::identify::distributions::Distribution;
use genetic_privacy::population::node_pair::NodeIdPair;
use genetic_privacy::statistics::gamma::HurdleGammaParams;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        panic!("Incorrect number of arguments.");
    }
    let directory = &args[1];
    let labeled_nodes = get_labeled_nodes(directory);
    println!("Fitting distributions");
    let _params = fit::fit_directory(directory);
    let fnv_params = convert_hashmap(_params);
    // TODO: Actually calculate this
    let default_cryptic = HurdleGammaParams {shape: 0.0, scale: 0.0, zero_prob: 0.0};
    let distribution = Distribution {distributions: fnv_params,
        cryptic_distribution:default_cryptic,
        labeled_nodes: labeled_nodes};
    println!("Serializing distributions");
    let encoded = serialize(&distribution, Infinite).unwrap();
    let mut file = File::create(&args[2]).unwrap();
    println!("Writing distributions to file");
    file.write_all(&encoded).unwrap();
}

fn get_labeled_nodes(directory: &str) -> Vec<u32> {
    let mut ret = Vec::new();
    let paths = read_dir(directory).unwrap();
    for path in paths {
        let file = path.unwrap();
        let labeled_id = file.file_name().to_str().unwrap().parse::<u32>().unwrap();
        ret.push(labeled_id);
    }
    ret
}

fn convert_hashmap<V>(mut input_map: HashMap<NodeIdPair, V>) -> FnvHashMap<NodeIdPair, V> {
    let mut map = FnvHashMap::default();
    for (key, value) in input_map.drain() {
        map.insert(key, value);
    }
    map
}

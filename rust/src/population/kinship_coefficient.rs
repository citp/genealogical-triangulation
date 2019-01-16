use fnv::{FnvHashMap, FnvHashSet};
use std::time::Duration;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::cell::Ref;
use pbr::ProgressBar;

use population::node::{Node, Population};
use population::founders::founders_id;
use util::set_intersection::nonzero_intersection_simd;

#[derive(PartialEq, Eq, PartialOrd, Ord, Hash, Copy, Clone)]
pub struct KinshipKey {
    lower: u32,
    higher: u32
}

impl KinshipKey {
    pub fn new(id1: u32, id2: u32) -> KinshipKey {
        if id1 < id2 {
            KinshipKey {lower: id1, higher: id2}
        } else {
            KinshipKey {lower: id2, higher: id1}
        }
    }

    pub fn values(&self) -> (u32, u32) {
        (self.lower, self.higher)
    }
}


fn clear_entry<T: Write>(map: &mut FnvHashMap<KinshipKey, f64>, key: &KinshipKey, node_generation: u32, keep_last: u32, last_generation: u32, outfile: &mut BufWriter<T>) {
    if last_generation - node_generation < keep_last {
        let kinship_coef = map.get(&key).unwrap();
        let (lower, upper) = key.values();
        write!(outfile, "{}\t{}\t{}\n", lower, upper, kinship_coef).unwrap();
    }
    map.remove(&key);
}

fn get_partners(population_members: &Vec<Ref<'_, Node>>) -> FnvHashSet<KinshipKey> {
    let mut ret = FnvHashSet::default();
    for node in population_members {
        if let (Some(mother), Some(father)) = (node.mother.upgrade(), node.father.upgrade()) {
            let key = KinshipKey::new(mother.borrow().id, father.borrow().id);
            ret.insert(key);

        }
    }
    ret
}

pub fn calculate_kinship_to_file(population: &Population, keep_last: u32, filename: &str) {
    assert!(keep_last > 2);
    let mut ret = FnvHashMap::default();

    let borrowed: Vec<_> = population.members.iter().map(|x| x.borrow()).collect();
    for i in 0..borrowed.len() - 1 {
        assert!(borrowed[i].generation <= borrowed[i+1].generation);
        assert!(borrowed[i].id < borrowed[i+1].id);
    }
    let last_generation = borrowed.last().unwrap().generation;
    let partners = get_partners(&borrowed);

    let outfile = File::create(filename).unwrap();
    let mut writer = BufWriter::new(outfile);

    let mut pb = ProgressBar::new(((borrowed.len() * (borrowed.len() + 1)) / 2) as u64);
    pb.set_max_refresh_rate(Some(Duration::from_secs(2)));

    for (i, node_2) in borrowed.iter().enumerate() {
        let self_key = KinshipKey::new(node_2.id, node_2.id);

        if let (Some(mother), Some(father)) = (node_2.mother.upgrade(), node_2.father.upgrade()) {
            let parent_key = KinshipKey::new(mother.borrow().id, father.borrow().id);
            let parent_kinship = *ret.get(&parent_key).unwrap_or(&0.0);

            let coeff = 0.5 + 0.5 * parent_kinship;
            ret.insert(self_key, coeff);
        } else {
            ret.insert(self_key, 0.5);
        }
        for node_1 in borrowed[i + 1..].iter() {
            if node_2.generation + 2 < node_1.generation {
                break;
            }

            if let (Some(mother), Some(father)) = (node_1.mother.upgrade(), node_1.father.upgrade()) {
                let mother_ref = mother.borrow();
                let mother_key = KinshipKey::new(mother_ref.id, node_2.id);
                let mother_coeff = *ret.get(&mother_key).unwrap_or(&0.0);

                let father_ref = father.borrow();
                let father_key = KinshipKey::new(father.borrow().id, node_2.id);
                let father_coeff = *ret.get(&father_key).unwrap_or(&0.0);

                let coeff = 0.5 * (mother_coeff + father_coeff);
                if coeff > 0.0 {
                    let key = KinshipKey::new(node_1.id, node_2.id);
                    ret.insert(key, coeff);
                }


                if mother_coeff != 0.0 && mother_ref.id < node_2.id  && !partners.contains(&mother_key){
                    let max_mother_child_id = mother_ref.children.iter().map(|x| x.borrow().id).max().unwrap();
                    if max_mother_child_id == node_1.id {
                        clear_entry(&mut ret, &mother_key, mother_ref.generation, keep_last, last_generation, &mut writer);
                    }
                }

                if father_coeff != 0.0 && father_ref.id < node_2.id && !partners.contains(&father_key) {
                    let max_father_child_id = father_ref.children.iter().map(|x| x.borrow().id).max().unwrap();
                    if max_father_child_id == node_1.id {
                        clear_entry(&mut ret, &father_key, father_ref.generation, keep_last, last_generation, &mut writer);
                    }
                }
            } // If the parents on node_1 are none, we assume the nodes are unrelated
        }
        pb.add((borrowed.len() - i) as u64);
    }
    pb.finish_print("Done calculating kinship.");

    let mut write_pb = ProgressBar::new(ret.len() as u64);
    write_pb.set_max_refresh_rate(Some(Duration::from_secs(2)));
    println!("Writing remaining entries to disk.");
    for (i, (key, kinship_coef)) in ret.iter().enumerate() {
        let (lower, upper) = key.values();
        write!(writer, "{}\t{}\t{}\n", lower, upper, kinship_coef).unwrap();
        if 0 < i && i % 10000 == 0 {
            write_pb.add(10000);
        }
    }
    write_pb.finish();
}

fn clear_kinships(to_delete: &mut Vec<KinshipKey>, kinships: &mut FnvHashMap<KinshipKey, f64>) {
    let count = to_delete.len();
    for key in to_delete.iter() {
        assert!(kinships.remove(key).is_some());
    }
    println!("Deleted {} old entries. Currently {} entries in kinship dictionary.", count, kinships.len());
    to_delete.clear();
    to_delete.shrink_to_fit();
}


pub fn calculate_kinship(population: &Population, keep_last: u32) -> FnvHashMap<KinshipKey, f64> {
    assert!(keep_last > 2);
    let mut ret = FnvHashMap::default();
    let borrowed: Vec<_> = population.members.iter().map(|x| x.borrow()).collect();
    for i in 0..borrowed.len() - 1 {
        assert!(borrowed[i].generation <= borrowed[i+1].generation);
    }
    let last_generation = borrowed.last().unwrap().generation;
    let mut pb = ProgressBar::new(((borrowed.len() * (borrowed.len() + 1)) / 2) as u64);
    pb.set_max_refresh_rate(Some(Duration::from_secs(2)));
    let mut to_delete: Vec<Vec<_>> = Vec::new();
    to_delete.resize((last_generation + 1) as usize, Vec::new());

    for (i, node_2) in borrowed.iter().enumerate() {
        let current_generation = node_2.generation;
        let keep = last_generation < current_generation + keep_last;
        if current_generation > 1 && to_delete[(current_generation - 2) as usize].len() > 0 {
            clear_kinships(&mut to_delete[(current_generation - 2) as usize], &mut ret);
        }

        for node_1 in borrowed[i..].iter() {
            if node_2.generation + 2 < node_1.generation {
                break;
            }

            if node_1.id == node_2.id {
                let key = KinshipKey::new(node_1.id, node_2.id);

                if let Some(mother) = node_1.mother.upgrade() {
                    let father = node_1.father.upgrade().unwrap();
                    let parent_key = KinshipKey::new(mother.borrow().id, father.borrow().id);
                    let parent_kinship = *ret.get(&parent_key).unwrap_or(&0.0);

                    let coeff = 0.5 + 0.5 * parent_kinship;
                    ret.insert(key, coeff);
                } else {
                    ret.insert(key, 0.5);
                }
                if !keep {
                    to_delete[current_generation as usize].push(key);
                }
            } else {
                if let Some(mother) = node_1.mother.upgrade() {

                    let key1 = KinshipKey::new(mother.borrow().id, node_2.id);
                    let mother_coeff = *ret.get(&key1).unwrap_or(&0.0);

                    let father = node_1.father.upgrade().unwrap();
                    let key2 = KinshipKey::new(father.borrow().id, node_2.id);
                    let father_coeff = *ret.get(&key2).unwrap_or(&0.0);

                    let coeff = 0.5 * (mother_coeff + father_coeff);
                    if coeff > 0.0 {
                        let key = KinshipKey::new(node_1.id, node_2.id);
                        ret.insert(key, coeff);
                        if !keep {
                            to_delete[current_generation as usize].push(key);
                        }
                    }
                } // If the parents on node_1 are none, we assume the nodes are unrelated
            }
        }
        pb.add((borrowed.len() - i) as u64);
    }
    for mut delete_vec in to_delete {
        if delete_vec.len() > 0 {
            clear_kinships(&mut delete_vec, &mut ret)
        }
    }
    ret
}


pub fn calculate_kinship_to_file_recursive(population: &Population, keep_last: u32, filename: &str) {
    let borrowed: Vec<_> = population.members.iter().map(|x| x.borrow()).collect();
    for i in 0..borrowed.len() - 1 {
        assert!(borrowed[i].generation <= borrowed[i+1].generation);
        assert!(borrowed[i].id < borrowed[i+1].id);
    }

    let outfile = File::create(filename).unwrap();
    let mut writer = BufWriter::new(outfile);

    let last_generation = borrowed.last().unwrap().generation;
    let to_calculate: Vec<_> = borrowed.iter().filter(|x| last_generation - x.generation < keep_last).collect();

    let mut pb = ProgressBar::new(((to_calculate.len() * (to_calculate.len() + 1)) / 2) as u64);
    pb.set_max_refresh_rate(Some(Duration::from_secs(2)));
    for (i, node_1) in to_calculate.iter().enumerate() {
        for node_2 in to_calculate[i..].iter() {
            let kinship = recursive_kinship(node_1, node_2);
            if 0.0 < kinship {
                write!(writer, "{}\t{}\t{}\n", node_1.id, node_2.id, kinship).unwrap();
            }
            pb.inc();
        }
    }
    pb.finish();

}

pub fn calculate_kinship_to_file_recursive_subset(population: &Population, keep_ids: Option<&[u32]>, keep_last: u32, filename: &str) {
    let borrowed: Vec<_> = population.members.iter().map(|x| x.borrow()).collect();
    for i in 0..borrowed.len() - 1 {
        assert!(borrowed[i].generation <= borrowed[i+1].generation);
        assert!(borrowed[i].id < borrowed[i+1].id);
        assert!(borrowed[i].id == i as u32);
    }


    //let test: () = borrowed;

    println!("Calculating founders.");
    let mut founder_pb = ProgressBar::new(borrowed.len() as u64);
    founder_pb.set_max_refresh_rate(Some(Duration::from_secs(1)));
    let founders: Vec<_> = borrowed.iter().map(|node| {
        let mut f = founders_id(node);
        f.sort();
        assert!(f.len() > 0);
        founder_pb.inc();
        f
    }).collect();
    founder_pb.finish();

    let outfile = File::create(filename).unwrap();
    let mut writer = BufWriter::new(outfile);

    let last_generation = borrowed.last().unwrap().generation;
    let potential_nodes: Vec<_> = borrowed.iter().filter(|x| last_generation - x.generation < keep_last).collect();

    let keep_nodes = match keep_ids {
        Some(keep_ids) => keep_ids.iter().map(|id| &borrowed[*id as usize]).collect(),
        None => potential_nodes.clone()
     };

    let keep_nodes_specified = keep_ids.is_some();

    let mut pb = if keep_nodes_specified {
        ProgressBar::new((keep_nodes.len() * potential_nodes.len()) as u64)
    } else {
        ProgressBar::new(((potential_nodes.len() * (potential_nodes.len() + 1)) / 2) as u64)
    };
    pb.set_max_refresh_rate(Some(Duration::from_secs(2)));

    for  (i, node_1) in keep_nodes.iter().enumerate() {
        let to_calculate = if keep_nodes_specified {
            potential_nodes.iter()
        } else {
            potential_nodes[i..].iter()
        };
        for node_2 in to_calculate {
            let kinship = recursive_kinship_with_founders(node_1, node_2, &founders);
            if 0.0 < kinship {
                write!(writer, "{}\t{}\t{}\n", node_1.id, node_2.id, kinship).unwrap();
            }
            pb.inc();
        }
    }
    pb.finish();

}

fn recursive_kinship(node_1: &Node, node_2: &Node) -> f64 {
    if node_1.id == node_2.id {
        if let (Some(mother), Some(father)) = (node_1.mother.upgrade(), node_1.father.upgrade()) {
            return 0.5 + 0.5 * recursive_kinship(&mother.borrow(), &father.borrow());
        } else {
            return 0.5
        }
    }

    let (lower, upper) = if node_1.id < node_2.id {
        (node_1, node_2)
    } else {
        (node_2, node_1)
    };

    if let (Some(mother), Some(father)) = (upper.mother.upgrade(), upper.father.upgrade()) {
        0.5 * recursive_kinship(lower, &mother.borrow()) + 0.5 * recursive_kinship(lower, &father.borrow())
    } else {
        0.0
    }
}

fn recursive_kinship_with_founders(node_1: &Node, node_2: &Node, founders: &[Vec<u32>]) -> f64 {
    let founders_1 = &founders[node_1.id as usize];
    let founders_2 = &founders[node_2.id as usize];
    if !nonzero_intersection_simd(founders_1, founders_2) {
        return 0.0;
    }
    if node_1.id == node_2.id {
        if let (Some(mother), Some(father)) = (node_1.mother.upgrade(), node_1.father.upgrade()) {
            return 0.5 + 0.5 * recursive_kinship(&mother.borrow(), &father.borrow());
        } else {
            return 0.5
        }
    }

    if let Some(twin) = node_1.twin.upgrade() {
        if twin.borrow().id == node_2.id {
            return recursive_kinship_with_founders(node_1, node_1, founders);
        }
    }

    let (lower, upper) = if node_1.id < node_2.id {
        (node_1, node_2)
    } else {
        (node_2, node_1)
    };

    if let (Some(mother), Some(father)) = (upper.mother.upgrade(), upper.father.upgrade()) {
        0.5 * recursive_kinship_with_founders(lower, &mother.borrow(), founders) + 0.5 * recursive_kinship_with_founders(lower, &father.borrow(), founders)
    } else {
        // Founder intersection should handle this case
        unreachable!();
    }
}
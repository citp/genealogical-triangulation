use fnv::FnvHashSet;

use population::node::Node;

pub fn founders_id(node: &Node) -> Vec<u32> {

    let mut to_visit;
    if let (Some(mother), Some(father)) = (node.mother.upgrade(), node.father.upgrade()) {
        to_visit = Vec::new();
        to_visit.push(mother.clone());
        to_visit.push(father.clone());
    } else {
        return vec![node.id]
    }

    let mut founders = FnvHashSet::default();
    while let Some(current) = to_visit.pop() {
        let current_ref = current.borrow();
        if let (Some(mother), Some(father)) = (current_ref.mother.upgrade(), current_ref.father.upgrade()) {
            to_visit.push(mother.clone());
            to_visit.push(father.clone());
        } else {
            founders.insert(current_ref.id);
        }
    }

    founders.iter().map(|x| *x).collect()
}
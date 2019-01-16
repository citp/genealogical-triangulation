use std::collections::VecDeque;
use std::rc::Rc;
use std::cell::RefCell;
use std::hash::{Hasher, Hash};

use fnv::FnvHashSet;

use population::node::Node;

struct NodeKey {
    node: Rc<RefCell<Node>>
}

impl NodeKey {
    fn new(node: Rc<RefCell<Node>>) -> NodeKey {
        NodeKey {node}
    }
}

impl From<Rc<RefCell<Node>>> for NodeKey {
    fn from(node: Rc<RefCell<Node>>) -> NodeKey {
        NodeKey::new(node)
    }
}

impl Hash for NodeKey {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.node.borrow().hash(state);
    }
}

impl PartialEq for NodeKey {
    fn eq(&self, other: &NodeKey) -> bool {
        self.node.borrow().eq(&other.node.borrow())
    }
}

impl Eq for NodeKey {}

pub fn descendants(node: &Node) -> Vec<Rc<RefCell<Node>>> {
    let mut ret: FnvHashSet<NodeKey> = FnvHashSet::default();
    let mut to_visit:VecDeque<Rc<RefCell<Node>>> = VecDeque::new();
    for child in node.children.iter() {
        to_visit.push_back(child.clone());
    }
    while let Some(current) = to_visit.pop_front() {
        ret.insert(current.clone().into());
        let current_ref = current.borrow();
        for child in current_ref.children.iter() {
            to_visit.push_back(child.clone());
        }
    }
    ret.into_iter().map(|x| x.node).collect()
}

pub fn descendants_id(node: &Node) -> Vec<u32> {
    descendants(node).iter().map(|x| x.borrow().id).collect()
}
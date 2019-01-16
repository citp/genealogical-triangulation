use std::rc::Rc;
use std::cell::RefCell;

use population::node::Node;

#[derive(PartialEq, Eq, Hash, Clone, Copy, Debug, Serialize, Deserialize)]
pub struct NodeIdPair {
    pub labeled: u32,
    pub unlabeled: u32
}

pub struct RelatedPair {
    pub labeled: Rc<RefCell<Node>>,
    pub unlabeled: Rc<RefCell<Node>>
}
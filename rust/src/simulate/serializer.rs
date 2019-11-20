use std::path::Path;
use std::io::{BufWriter, Write};
use std::fs::{remove_file, OpenOptions, File};

use byteorder::{LittleEndian, WriteBytesExt};

use fnv::FnvHashMap;
use std::time::Instant;
use population::node_pair::RelatedPair;
use std::mem;

pub struct SimulationSerializer {
    buffers: FnvHashMap<u32, Vec<(u32, f64)>>,
    writer: BufWriter<File>
}

impl SimulationSerializer {
    pub fn new(outfile: &str, clobber: bool, related_pairs: &[RelatedPair]) -> Self {
        let write_header = !Path::new(outfile).exists() || clobber;
        let writer = set_up_file(outfile, clobber);
        let mut ret = SimulationSerializer {buffers: FnvHashMap::default(), writer: writer};
        if write_header {
            ret.write_header(related_pairs);
        }
        ret
    }

    pub fn insert(&mut self, anchor: u32, unlabeled: u32, shared: f64) {
        let buffer = self.buffers.entry(anchor).or_insert_with(Vec::new);
        buffer.push((unlabeled, shared));
    }

    fn write_header(&mut self, related_pairs: &[RelatedPair]) {
        let mut id_pairs: Vec<_> = related_pairs.iter()
            .map(|pair| (pair.labeled.borrow().id, pair.unlabeled.borrow().id))
            .collect();
        id_pairs.sort_unstable();
        let mut pair_bytes = Vec::with_capacity(id_pairs.len() * 2 * mem::size_of::<u32>());
        for (anchor, unlabeled) in id_pairs {
            pair_bytes.write_u32::<LittleEndian>(anchor).unwrap();
            pair_bytes.write_u32::<LittleEndian>(unlabeled).unwrap();
        }
        self.writer.write_u64::<LittleEndian>(pair_bytes.len() as u64).unwrap();
        self.writer.write_all(&pair_bytes).unwrap();
        self.writer.flush().unwrap();
    }

    pub fn flush(&mut self) {
        let mut buffers: Vec<_> = self.buffers.drain().collect();
        buffers.sort_unstable_by_key(|(anchor, _)| *anchor);
        let now = Instant::now();
        for (_anchor, unlabeled_shared) in buffers.iter_mut() {
            unlabeled_shared.sort_unstable_by_key(|(unlabeled, _shared)| *unlabeled);
            let shared: Vec<_> = unlabeled_shared.iter().map(|(_, shared)| *shared).collect();
            for s in shared {
                self.writer.write_f64::<LittleEndian>(s).unwrap();
            }
        }
        self.writer.flush().unwrap();
        println!("It took {} seconds to write to file.", now.elapsed().as_secs());
    }
}

fn set_up_file(filename: &str, clobber: bool) -> BufWriter<File> {
    let path = Path::new(filename);
    if clobber && path.exists() {
        remove_file(&path).unwrap();
    }
    let file;
    if path.exists() {
        file = OpenOptions::new().append(true).open(path).unwrap();
    } else {
        file = OpenOptions::new().create(true).write(true).open(path).unwrap();
    }
    BufWriter::new(file)
}
extern crate genetic_privacy;
extern crate rand;

use rand::thread_rng;
// TODO: Fix depricated
#[allow(deprecated)]
use rand::distributions::{IndependentSample, Range};

use genetic_privacy::genome::diploid::Diploid;
use genetic_privacy::genome::common_segments::common_homolog_segments;

fn main() {
    let founders = Range::new(0, 2);
    for _ in 0..1000 {
        let a = vector_sequence(100000, &founders);
        let b = vector_sequence(100000, &founders);
        let diploid_a = vector_sequence_to_diploid(&a);
        let diploid_b = vector_sequence_to_diploid(&b);
        let vec_equal = compare(&a, &b);
        let common_segments = common_homolog_segments(&diploid_a,
                                                      &diploid_b);
        assert_eq!(vec_equal, common_segments);
    }
}

fn compare(a: &[u8], b: &[u8]) -> Vec<(u32, u32)> {
    let equal = a.iter().zip(b).map(|(&a, &b)| a == b);
    let mut start = 0;
    let mut run = false;
    let mut res = Vec::new();
    for (i, current_equal) in equal.enumerate() {
        if !run && current_equal {
            run = true;
            start = i;
        } else if !current_equal && run {
            run = false;
            res.push((start as u32, i as u32));
        }
    }
    if run {
        res.push((start as u32, a.len() as u32));
    }
    res
}

fn vector_sequence(length: usize, founder_range: &Range<u8>) -> Vec<u8> {
    let mut res = Vec::with_capacity(length);
    let mut rng = thread_rng();
    for _ in 0..length {
        #[allow(deprecated)]
        res.push(founder_range.ind_sample(&mut rng))
    }
    res
}

fn vector_sequence_to_diploid(sequence: &[u8]) -> Diploid {
    let end = sequence.len() as u32;
    let mut starts = Vec::new();
    let mut founder = Vec::new();
    let mut current = sequence[0] as u32;
    starts.push(0);
    founder.push(current);
    for (i, &val) in sequence.iter().enumerate() {
        if val as u32 != current {
            current = val as u32;
            starts.push(i as u32);
            founder.push(current);
        }
    }
    Diploid {starts: starts, founder: founder, end: end}
}

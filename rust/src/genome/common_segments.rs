use genome::genome::{Genome};
use genome::diploid::Diploid;

pub fn common_segment_lengths(a: &Genome, b: &Genome) -> Vec<u32> {
    let mut ret = Vec::new();
    ret.extend_from_slice(&lengths(&common_homolog_segments(&a.mother,
                                                            &b.mother)));
    ret.extend_from_slice(&lengths(&common_homolog_segments(&a.father,
                                                            &b.mother)));
    ret.extend_from_slice(&lengths(&common_homolog_segments(&a.mother,
                                                            &b.father)));
    ret.extend_from_slice(&lengths(&common_homolog_segments(&a.father,
                                                            &b.father)));
    ret
}

pub fn shared_segment_length_genomes(a: &Genome, b: &Genome, min_length: u32)
                                 -> u64 {
    let lengths = common_segment_lengths(a, b);
    lengths.iter().filter(|&length| *length >= min_length)
        .map(|&x| x as u64).sum()
}

pub fn common_homolog_segments(homolog_a: &Diploid, homolog_b: &Diploid) -> Vec<(u32, u32)> {
    let mut shared_segments = Vec::new();
    let mut index_a = 0;
    let mut index_b = 0;
    let mut start;
    let mut stop;
    let starts_a = &homolog_a.starts;
    let starts_b = &homolog_b.starts;
    let founder_a = &homolog_a.founder;
    let founder_b = &homolog_b.founder;
    while index_a < starts_a.len() && index_b < starts_b.len() {
        let a_start = starts_a[index_a];
        let a_stop;
        if index_a + 1 < starts_a.len() {
            a_stop = starts_a[index_a + 1]
        } else {
            a_stop = homolog_a.end
        }
        let a_id = founder_a[index_a];

        let b_start = starts_b[index_b];
        let b_stop;
        if index_b + 1 < starts_b.len() {
            b_stop = starts_b[index_b + 1]
        } else {
            b_stop = homolog_b.end
        }
        let b_id = founder_b[index_b];

        if a_id == b_id {
            if a_start > b_start {
                start = a_start
            } else {
                start = b_start
            }
            if a_stop < b_stop {
                stop = a_stop;
            } else {
                stop = b_stop
            }
            shared_segments.push((start, stop));
        }
        if a_stop == b_stop {
            index_a += 1;
            index_b += 1;
        } else if a_stop > b_stop {
            index_b += 1;
        } else {
            index_a += 1;
        }
    }
    if shared_segments.len() <= 1 {
        return shared_segments;
    }
    consolidate_sequence(&shared_segments)
}

fn lengths(segments: &[(u32, u32)]) -> Vec<u32> {
    segments.iter().map(|&(a, b)| b - a).collect()
}

pub fn consolidate_sequence(sequence: &[(u32, u32)]) -> Vec<(u32, u32)> {
    assert!(sequence.len() > 1);
    let mut consolidated = Vec::with_capacity(sequence.len());
    let mut i = 0;
    let mut j = 1;
    while j < sequence.len() {
        if sequence[j - 1].1 != sequence[j].0 {
            consolidated.push((sequence[i].0, sequence[j - 1].1));
            i = j;
        }
        j += 1
    }
    consolidated.push((sequence[i].0, sequence[j-1].1));
    consolidated
}

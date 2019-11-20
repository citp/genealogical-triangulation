use std::cmp::{Ordering, max, min};

use genome::genome::{Genome};
use genome::diploid::Diploid;
use genome::cm::CmConverter;

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

pub fn common_segments_inbreeding(a: &Genome, b: &Genome) -> Vec<(u32, u32)> {
    let ibd_a_mother_b_mother = common_homolog_segments(&a.mother,&b.mother);
    let ibd_a_father_b_mother = common_homolog_segments(&a.father,&b.mother);
    let mut ibd_b_mother = merge_overlaps(ibd_a_mother_b_mother, ibd_a_father_b_mother);

    let ibd_a_mother_b_father = common_homolog_segments(&a.mother,&b.father);
    let ibd_a_father_b_father = common_homolog_segments(&a.father,&b.father);
    let mut ibd_b_father = merge_overlaps(ibd_a_mother_b_father, ibd_a_father_b_father);

    let mut b_inbreed =common_homolog_segments(&b.mother, &b.father);
    if 0 < b_inbreed.len() {
        let a_inbreed = common_homolog_segments(&a.mother, &a.father);
        subtract_regions(&mut b_inbreed, &a_inbreed);
        remove_inbreeding(&mut ibd_b_mother, &mut ibd_b_father, &b_inbreed)
    } else {
        let mut ret = Vec::with_capacity(ibd_b_mother.len() + ibd_b_father.len());
        ret.extend_from_slice(&ibd_b_mother);
        ret.extend_from_slice(&ibd_b_father);
        ret
    }

}

fn merge_overlaps(mut ibd_a: Vec<(u32, u32)>, mut ibd_b: Vec<(u32, u32)>) -> Vec<(u32, u32)> {
    let total_length = ibd_a.len() + ibd_b.len();
    if total_length == 0 {
        return Vec::new();
    }
    let mut all_ibd = Vec::with_capacity(ibd_a.len() + ibd_b.len());
    all_ibd.append(&mut ibd_a);
    all_ibd.append(&mut ibd_b);
    all_ibd.sort_unstable_by(|a, b| compare_tuples(a, b).reverse());
    let mut ret = Vec::with_capacity(all_ibd.len());
    ret.push(all_ibd.pop().unwrap());
    while let Some(current) = all_ibd.pop() {
        let top = *ret.last().unwrap();
        if current.0 <= top.1 {
            let end = max(current.1, top.1);
            *ret.last_mut().unwrap() = (top.0, end);
        } else {
            ret.push(current);
        }
    }
    ret
}

fn remove_inbreeding(ibd_a: &mut Vec<(u32, u32)>, ibd_b: &mut Vec<(u32, u32)>, inbreeding: &[(u32, u32)]) -> Vec<(u32, u32)> {
    for region in inbreeding {
        let a_overlap = size_of_overlap(ibd_a, *region);
        let b_overlap = size_of_overlap(ibd_b, *region);
        assert_eq!(a_overlap, b_overlap);
        if a_overlap == 0 || b_overlap == 0 {
            continue;
        }
        if a_overlap < b_overlap {
            subtract_region(ibd_a, region);
        } else {
            subtract_region(ibd_b, region);
        }
    }
    let mut ret = Vec::with_capacity(ibd_a.len() + ibd_b.len());
    ret.append(ibd_a);
    ret.append(ibd_b);
    ret

}

fn compare_tuples(a: &(u32, u32), b: &(u32, u32)) -> Ordering {
    let first = a.0.cmp(&b.0);
    if first == Ordering::Equal {
        a.1.cmp(&b.1)
    } else {
        first
    }
}

///Subtract the regions in b from the regions in a
fn subtract_regions(a: &mut Vec<(u32, u32)>, b: &Vec<(u32, u32)>) {
    for region in b {
        subtract_region(a, region);
    }
}

fn subtract_region(a: &mut Vec<(u32, u32)>, region: &(u32, u32)) {
    let mut start_i: i32 = 0;
    if a.len() == 0 {
        return;
    }
    let mut stop_i = a.len() as i32 - 1;
    while (start_i as usize) < a.len() && a[start_i as usize].1 <= region.0 {
        start_i += 1;
    }

    while 0 <= stop_i && region.1 <= a[stop_i as usize].0 {
        stop_i -= 1;
    }

    if stop_i < start_i {
        return;
    }

    let mut to_insert = Vec::with_capacity(2);
    if a[start_i as usize].0 < region.0 {
        to_insert.push((a[start_i as usize].0, region.0));
    }
    if region.1 < a[stop_i as usize].1 {
        to_insert.push((region.1, a[stop_i as usize].1));
    }
    a.drain(start_i as usize ..stop_i as usize + 1);
    while let Some(insert) = to_insert.pop() {
        a.insert(start_i as usize, insert);
    }
}

pub fn size_of_overlap(regions: &[(u32, u32)], query_region: (u32, u32)) -> u32 {
    let mut overlap = 0;
    for region in regions {
        let start = max(region.0, query_region.0);
        let end = min(region.1, query_region.1);
        if start < end {
            overlap += end - start;
        }
    }
    overlap
}

pub fn shared_segment_length_genomes(a: &Genome, b: &Genome, cm_converter: &CmConverter) -> f64 {
    let lengths = common_segments_inbreeding(a, b);
    cm_converter.cm_lengths(&lengths).iter().sum()
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

#[test]
fn test_size_of_overlap() {
    assert_eq!(size_of_overlap(&vec![(1, 2)], (1, 2)), 1);
    assert_eq!(size_of_overlap(&vec![(1, 10)], (5, 10)), 5);
    assert_eq!(size_of_overlap(&vec![(1, 10)], (1, 5)), 4);
    assert_eq!(size_of_overlap(&vec![(1, 10)], (2, 8)), 6);
    assert_eq!(size_of_overlap(&vec![(1, 10)], (15, 20)), 0);
    assert_eq!(size_of_overlap(&vec![(10, 20)], (1, 10)), 0);

    assert_eq!(size_of_overlap(&vec![(1, 10), (20, 30), (40, 50)], (1, 10)), 9);
    assert_eq!(size_of_overlap(&vec![(1, 10), (20, 30), (40, 50)], (20, 30)), 10);
    assert_eq!(size_of_overlap(&vec![(1, 10), (20, 30), (40, 50)], (40, 50)), 10);

    assert_eq!(size_of_overlap(&vec![(1, 10), (20, 30), (40, 50)], (60, 70)), 0);
    assert_eq!(size_of_overlap(&vec![(1, 10), (20, 30), (40, 50)], (0, 1)), 0);
    assert_eq!(size_of_overlap(&vec![(1, 10), (20, 30), (40, 50)], (12, 18)), 0);


    assert_eq!(size_of_overlap(&vec![(1, 10), (20, 30), (40, 50)], (5, 10)), 5);

    assert_eq!(size_of_overlap(&vec![(1, 10), (20, 30), (40, 50)], (5, 12)), 5);

    assert_eq!(size_of_overlap(&vec![(1, 10), (20, 30), (40, 50)], (5, 25)), 10);
    assert_eq!(size_of_overlap(&vec![(1, 10), (20, 30), (40, 50)], (5, 45)), 20);
    assert_eq!(size_of_overlap(&vec![(1, 10), (20, 30), (40, 50)], (5, 50)), 25);
}

#[test]
fn test_subtract_region() {
    let mut input = vec![(1, 2)];
    subtract_region(&mut input, &(1, 2));
    assert_eq!(input, vec![]);

    input = vec![(1, 10)];
    subtract_region(&mut input, &(5, 10));
    assert_eq!(input, vec![(1, 5)]);

    input = vec![(1, 10)];
    subtract_region(&mut input, &(1, 5));
    assert_eq!(input, vec![(5, 10)]);

    input = vec![(1, 10)];
    subtract_region(&mut input, &(2, 8));
    assert_eq!(input, vec![(1, 2), (8, 10)]);

    input = vec![(1, 10)];
    subtract_region(&mut input, &(15, 20));
    assert_eq!(input, vec![(1, 10)]);

    input = vec![(10, 20)];
    subtract_region(&mut input, &(1, 10));
    assert_eq!(input, vec![(10, 20)]);

    input = vec![(1, 10), (20, 30), (40, 50)];
    subtract_region(&mut input, &(1, 10));
    assert_eq!(input, vec![(20, 30), (40, 50)]);

    input = vec![(1, 10), (20, 30), (40, 50)];
    subtract_region(&mut input, &(20, 30));
    assert_eq!(input, vec![(1, 10), (40, 50)]);

    input = vec![(1, 10), (20, 30), (40, 50)];
    subtract_region(&mut input, &(40, 50));
    assert_eq!(input, vec![(1, 10), (20, 30)]);

    input = vec![(1, 10), (20, 30), (40, 50)];
    subtract_region(&mut input, &(5, 10));
    assert_eq!(input, vec![(1, 5), (20, 30), (40, 50)]);

    input = vec![(1, 10), (20, 30), (40, 50)];
    subtract_region(&mut input, &(5, 12));
    assert_eq!(input, vec![(1, 5), (20, 30), (40, 50)]);

    input = vec![(1, 10), (20, 30), (40, 50)];
    subtract_region(&mut input, &(5, 25));
    assert_eq!(input, vec![(1, 5), (25, 30), (40, 50)]);

    input = vec![(1, 10), (20, 30), (40, 50)];
    subtract_region(&mut input, &(5, 45));
    assert_eq!(input, vec![(1, 5), (45, 50)]);

    input = vec![(1, 10), (20, 30), (40, 50)];
    subtract_region(&mut input, &(5, 50));
    assert_eq!(input, vec![(1, 5)]);

}

#[test]
fn test_merge_overlaps() {
    assert_eq!(merge_overlaps(vec![], vec![]), vec![]);
    assert_eq!(merge_overlaps(vec![], vec![(1, 2)]), vec![(1, 2)]);
    assert_eq!(merge_overlaps(vec![(1, 2)], vec![]), vec![(1, 2)]);

    assert_eq!(merge_overlaps(vec![(1, 2)], vec![(1, 2)]), vec![(1, 2)]);
    assert_eq!(merge_overlaps(vec![(1, 2)], vec![(4, 5)]), vec![(1, 2), (4, 5)]);
    assert_eq!(merge_overlaps(vec![(4, 5)], vec![(1, 2)]), vec![(1, 2), (4, 5)]);
    assert_eq!(merge_overlaps(vec![(1, 2)], vec![(1, 4)]), vec![(1, 4)]);
    assert_eq!(merge_overlaps(vec![(1, 4)], vec![(1, 2)]), vec![(1, 4)]);

    assert_eq!(merge_overlaps(vec![(1, 2), (5, 8), (10, 15), (20, 25)], vec![(9, 18)]),
               vec![(1, 2), (5, 8), (9, 18), (20, 25)]);

}

#[test]
fn test_remove_inbreeding() {
    assert_eq!(remove_inbreeding(&mut vec![(10, 20)], &mut vec![(10, 20)], &vec![(10, 20)]),
               vec![(10, 20)]);
    assert_eq!(remove_inbreeding(&mut vec![(10, 20)], &mut vec![(10, 20)], &vec![]),
               vec![(10, 20), (10, 20)]);
    assert_eq!(remove_inbreeding(&mut vec![(10, 20)], &mut vec![(10, 20)], &vec![(30, 40)]),
               vec![(10, 20), (10, 20)]);
    assert_eq!(remove_inbreeding(&mut vec![(10, 20), (20, 30)], &mut vec![(10, 20)], &vec![(10, 20)]),
               vec![(10, 20), (20, 30)]);
    assert_eq!(remove_inbreeding(&mut vec![(10, 20)], &mut vec![(10, 20), (20, 30)], &vec![(10, 20)]),
               vec![(10, 20), (20, 30)]);
}
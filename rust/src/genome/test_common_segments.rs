#[allow(unused_imports)]
use genome::diploid::Diploid;
// use genome::genome::Genome;
#[allow(unused_imports)]
use genome::common_segments::{consolidate_sequence, common_homolog_segments};

// common_homolog_segments tests

#[test]
fn test_single_same_segment() {
    let a = Diploid {
        starts: vec![0],
        end: 10,
        founder: vec![1]
    };
    let shared = common_homolog_segments(&a, &a);
    assert_eq!(shared, [(0, 10)]);
}

#[test]
fn test_single_different_segment() {
    let a = Diploid {
        starts: vec![0],
        end: 10,
        founder: vec![0]
    };
    let b = Diploid {
        starts: vec![0],
        end: 10,
        founder: vec![1]
    };
    let shared = common_homolog_segments(&a, &b);
    assert_eq!(shared, []);
}

#[test]
fn test_multiple_same_segment() {
    let a = Diploid {
        starts: vec![0],
        end: 10,
        founder: vec![0]
    };
    let b = Diploid {
        starts: vec![0, 5],
        end: 10,
        founder: vec![0, 0]
    };
    let shared = common_homolog_segments(&a, &b);
    assert_eq!(shared, [(0, 10)]);
}

#[test]
fn test_two_same_segments() {
    let a = Diploid {
        starts: vec![0, 5],
        end: 10,
        founder: vec![1, 2]
    };
    let shared = common_homolog_segments(&a, &a);
    assert_eq!(shared, [(0, 10)]);
}

fn test_two_segments_different_boundary() {
    let a = Diploid {
        starts: vec![0, 5],
        end: 10,
        founder: vec![1, 1]
    };
    let b = Diploid {
        starts: vec![0, 6],
        end: 10,
        founder: vec![1, 1]
    };
    let shared = common_homolog_segments(&a, &b);
    assert_eq!(shared, [(0, 10)]);
}

fn test_single_element_vs_many_match_in_back() {
    let a = Diploid {
        starts: vec![0],
        end: 10,
        founder: vec![0]
    };
    let b = Diploid {
        starts: vec![0, 2, 4, 8],
        end: 10,
        founder: vec![1, 2, 3, 0]
    };
    let shared = common_homolog_segments(&a, &b);
    assert_eq!(shared, [(8, 10)]);
}

fn test_single_element_vs_many_match_in_front() {
    let a = Diploid {
        starts: vec![0],
        end: 10,
        founder: vec![0]
    };
    let b = Diploid {
        starts: vec![0, 2, 4, 8],
        end: 10,
        founder: vec![0, 2, 3, 4]
    };
    let shared = common_homolog_segments(&a, &b);
    assert_eq!(shared, [(0, 2)]);
}


// consoloidate tests
#[test]
fn test_two_elements_merge() {
    let seq = [(0, 5), (5, 10)];
    let con = consolidate_sequence(&seq);
    assert_eq!(con, [(0, 10)]);
}

#[test]
fn test_two_elements_disjoint() {
    let seq = [(0, 5), (6, 10)];
    let con = consolidate_sequence(&seq);
    assert_eq!(con, seq);
}

#[test]
fn test_first_two_merge() {
    let seq = [(0, 4), (4, 8), (9, 10)];
    let con = consolidate_sequence(&seq);
    assert_eq!(con, [(0, 8), (9, 10)]);
}

#[test]
fn test_last_two_merge() {
    let seq = [(0, 4), (5, 8), (8, 10)];
    let con = consolidate_sequence(&seq);
    assert_eq!(con, [(0, 4), (5, 10)]);
}

#[test]
fn test_middle_two_merge() {
    let seq = [(0, 3), (4, 6), (6, 8), (9, 10)];
    let con = consolidate_sequence(&seq);
    assert_eq!(con, [(0, 3), (4, 8), (9, 10)]);
}

#[test]
fn test_many_elements() {
    let seq = [(0, 2), (2, 4), (4, 8), (8, 10)];
    let con = consolidate_sequence(&seq);
    assert_eq!(con, [(0, 10)]);
}

#[test]
fn test_many_emements_single_point() {
    let seq = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5),
               (5, 6), (6, 7), (7, 8), (8, 9), (9, 10)];
    let con = consolidate_sequence(&seq);
    assert_eq!(con, [(0, 10)]);
}

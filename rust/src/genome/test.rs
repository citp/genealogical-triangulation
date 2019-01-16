#[cfg(test)]
#[allow(unused_imports)]
use genome::diploid::Diploid;
#[allow(unused_imports)]
use genome::genome::Genome;
#[allow(unused_imports)]
use genome::recombinator::{swap_at_locations, new_sequence, cum_sum};

#[test]
fn test_single_element_middle() {
    let diploid = Diploid {
        starts: vec![0],
        end: 10,
        founder: vec![1]
    };
    let new = new_sequence(&diploid, &vec![5]);
    assert_eq!(new.starts, vec![0, 5]);
    assert_eq!(new.founder, vec![1, 1]);
}

#[test]
fn test_single_element_start() {
    let diploid = Diploid {
        starts: vec![0],
        end: 10,
        founder: vec![1]
    };
    let new = new_sequence(&diploid, &vec![0]);
    assert_eq!(new.starts, vec![0]);
    assert_eq!(new.founder, vec![1]);
}

#[test]
fn test_single_element_multiple() {
    let diploid = Diploid {
        starts: vec![0],
        end: 10,
        founder: vec![1]
    };
    let new = new_sequence(&diploid, &vec![0, 4, 6, 8]);
    assert_eq!(new.starts, vec![0, 4, 6, 8]);
    assert_eq!(new.founder, vec![1, 1, 1, 1]);
}

#[test]
fn test_end_boundary_two_element() {
    let diploid = Diploid {
        starts: vec![0, 10],
        end: 20,
        founder: vec![1, 2]
    };
    let new = new_sequence(&diploid, &vec![10]);
    assert_eq!(new.starts, vec![0, 10]);
    assert_eq!(new.founder, vec![1, 2]);
}

#[test]
fn test_start_boundary_two_element() {
    let diploid = Diploid {
        starts: vec![0, 10],
        end: 20,
        founder: vec![1, 2]
    };
    let new = new_sequence(&diploid, &vec![0]);
    assert_eq!(new.starts, vec![0, 10]);
    assert_eq!(new.founder, vec![1, 2]);
}

#[test]
fn test_middle_boundary_two_element() {
    let diploid = Diploid {
        starts: vec![0, 10],
        end: 20,
        founder: vec![1, 2]
    };
    let new = new_sequence(&diploid, &vec![5]);
    assert_eq!(new.starts, vec![0, 5, 10]);
    assert_eq!(new.founder, vec![1, 1, 2]);
}

#[test]
fn test_middle_boundary_two_element_multiple_breaks() {
    let diploid = Diploid {
        starts: vec![0, 10],
        end: 20,
        founder: vec![1, 2]
    };
    let new = new_sequence(&diploid, &vec![5, 15]);
    assert_eq!(new.starts, vec![0, 5, 10, 15]);
    assert_eq!(new.founder, vec![1, 1, 2, 2]);
}

#[test]
fn test_single_element_end() {
    let diploid = Diploid {
        starts: vec![0],
        end: 10,
        founder: vec![1]
    };
    let new = new_sequence(&diploid, &vec![10]);
    assert_eq!(new.starts, vec![0]);
    assert_eq!(new.founder, vec![1]);
}

#[test]
fn test_empty_locations() {
    let diploid = Diploid {
        starts: vec![0],
        end: 10,
        founder: vec![1]
    };
    let new = new_sequence(&diploid, &vec![]);
    assert_eq!(new.starts, vec![0]);
    assert_eq!(new.founder, vec![1]);
}

    
#[test]
fn test_single_location_left_boundary() {
    let mother = Diploid {
        starts: vec![0],
        founder: vec![1],
        end: 10
    };
    let father = Diploid {
        starts: vec![0],
        founder: vec![2],
        end: 10
    };
    let genome = Genome {mother: mother, father: father};
    let locations = vec![0, 5];
    let res = swap_at_locations(&genome, &locations);
    let new_mother = res.mother;
    let new_father = res.father;
    assert_eq!(new_mother.starts, vec![0, 5]);
    assert_eq!(new_mother.founder, vec![2, 1]);
    assert_eq!(new_father.starts, vec![0, 5]);
    assert_eq!(new_father.founder, vec![1, 2]);
}

#[test]
fn test_single_location_right_boundary() {
    let mother = Diploid {
        starts: vec![0],
        founder: vec![1],
        end: 10
    };
    let father = Diploid {
        starts: vec![0],
        founder: vec![2],
        end: 10
    };
    let genome = Genome {mother: mother, father: father};
    let locations = vec![5, 10];
    let res = swap_at_locations(&genome, &locations);
    let new_mother = res.mother;
    let new_father = res.father;
    assert_eq!(new_mother.starts, vec![0, 5]);
    assert_eq!(new_mother.founder, vec![1, 2]);
    assert_eq!(new_father.starts, vec![0, 5]);
    assert_eq!(new_father.founder, vec![2, 1]);
}

#[test]
fn test_single_location_middle_boundary() {
    let mother = Diploid {
        starts: vec![0],
        founder: vec![1],
        end: 10
    };
    let father = Diploid {
        starts: vec![0],
        founder: vec![2],
        end: 10
    };
    let genome = Genome {mother: mother, father: father};
    let locations = vec![2, 8];
    let res = swap_at_locations(&genome, &locations);
    let new_mother = res.mother;
    let new_father = res.father;
    assert_eq!(new_mother.starts, vec![0, 2, 8]);
    assert_eq!(new_mother.founder, vec![1, 2, 1]);
    assert_eq!(new_father.starts, vec![0, 2, 8]);
    assert_eq!(new_father.founder, vec![2, 1, 2]);
}

#[test]
fn test_multiple_locations_single_segment() {
    let mother = Diploid {
        starts: vec![0],
        founder: vec![1],
        end: 10
    };
    let father = Diploid {
        starts: vec![0],
        founder: vec![2],
        end: 10
    };
    let genome = Genome {mother: mother, father: father};
    let locations = vec![0, 4, 6, 10];
    let res = swap_at_locations(&genome, &locations);
    let new_mother = res.mother;
    let new_father = res.father;
    assert_eq!(new_mother.starts, vec![0, 4, 6]);
    assert_eq!(new_mother.founder, vec![2, 1, 2]);
    assert_eq!(new_father.starts, vec![0, 4, 6]);
    assert_eq!(new_father.founder, vec![1, 2, 1]);
}

#[test]
fn test_single_location_two_segments_first_segment() {
    let mother = Diploid {
        starts: vec![0, 10],
        founder: vec![1, 2],
        end: 20
    };
    let father = Diploid {
        starts: vec![0, 10],
        founder: vec![3, 4],
        end: 20
    };
    let genome = Genome {mother: mother, father: father};
    let locations = vec![0, 10];
    let res = swap_at_locations(&genome, &locations);
    let new_mother = res.mother;
    let new_father = res.father;
    assert_eq!(new_mother.starts, vec![0, 10]);
    assert_eq!(new_mother.founder, vec![3, 2]);
    assert_eq!(new_father.starts, vec![0, 10]);
    assert_eq!(new_father.founder, vec![1, 4]);
}

#[test]
fn test_single_location_two_segments_last_segment() {
    let mother = Diploid {
        starts: vec![0, 10],
        founder: vec![1, 2],
        end: 20
    };
    let father = Diploid {
        starts: vec![0, 10],
        founder: vec![3, 4],
        end: 20
    };
    let genome = Genome {mother: mother, father: father};
    let locations = vec![10, 20];
    let res = swap_at_locations(&genome, &locations);
    let new_mother = res.mother;
    let new_father = res.father;
    assert_eq!(new_mother.starts, vec![0, 10]);
    assert_eq!(new_mother.founder, vec![1, 4]);
    assert_eq!(new_father.starts, vec![0, 10]);
    assert_eq!(new_father.founder, vec![3, 2]);
}

#[test]
fn test_single_location_overlapping_two_segments() {
    let mother = Diploid {
        starts: vec![0, 10],
        founder: vec![1, 2],
        end: 20
    };
    let father = Diploid {
        starts: vec![0, 10],
        founder: vec![3, 4],
        end: 20
    };
    let genome = Genome {mother: mother, father: father};
    let locations = vec![5, 15];
    let res = swap_at_locations(&genome, &locations);
    let new_mother = res.mother;
    let new_father = res.father;
    assert_eq!(new_mother.starts, vec![0, 5, 10, 15]);
    assert_eq!(new_mother.founder, vec![1, 3, 4, 2]);
    assert_eq!(new_father.starts, vec![0, 5, 10, 15]);
    assert_eq!(new_father.founder, vec![3, 1, 2, 4]);
}

#[test]
fn test_two_locations_at_boundary_two_segments() {
    let mother = Diploid {
        starts: vec![0, 10],
        founder: vec![1, 2],
        end: 20
    };
    let father = Diploid {
        starts: vec![0, 10],
        founder: vec![3, 4],
        end: 20
    };
    let genome = Genome {mother: mother, father: father};
    let locations = vec![0, 5, 15, 20];
    let res = swap_at_locations(&genome, &locations);
    let new_mother = res.mother;
    let new_father = res.father;
    assert_eq!(new_mother.starts, vec![0, 5, 10, 15]);
    assert_eq!(new_mother.founder, vec![3, 1, 2, 4]);
    assert_eq!(new_father.starts, vec![0, 5, 10, 15]);
    assert_eq!(new_father.founder, vec![1, 3, 4, 2]);
}

#[test]
fn test_two_locations_middle_two_segments() {
    let mother = Diploid {
        starts: vec![0, 10],
        founder: vec![1, 2],
        end: 20
    };
    let father = Diploid {
        starts: vec![0, 10],
        founder: vec![3, 4],
        end: 20
    };
    let genome = Genome {mother: mother, father: father};
    let locations = vec![2, 5, 15, 18];
    let res = swap_at_locations(&genome, &locations);
    let new_mother = res.mother;
    let new_father = res.father;
    assert_eq!(new_mother.starts, vec![0, 2, 5, 10, 15, 18]);
    assert_eq!(new_mother.founder, vec![1, 3, 1, 2, 4, 2]);
    assert_eq!(new_father.starts, vec![0, 2, 5, 10, 15, 18]);
    assert_eq!(new_father.founder, vec![3, 1, 3, 4, 2, 4]);
}

#[test]
fn test_single_location_two_segments_uneven() {
    let mother = Diploid {
        starts: vec![0, 5],
        founder: vec![1, 2],
        end: 20
    };
    let father = Diploid {
        starts: vec![0, 10],
        founder: vec![3, 4],
        end: 20
    };
    let genome = Genome {mother: mother, father: father};
    let locations = vec![2, 11];
    let res = swap_at_locations(&genome, &locations);
    let new_mother = res.mother;
    let new_father = res.father;
    assert_eq!(new_mother.starts, vec![0, 2, 10, 11]);
    assert_eq!(new_mother.founder, vec![1, 3, 4, 2]);
    assert_eq!(new_father.starts, vec![0, 2, 5, 11]);
    assert_eq!(new_father.founder, vec![3, 1, 2, 4]);
}

#[test]
fn test_empty_cum_sum() {
    let in_seq: Vec<u32> = vec![];
    assert_eq!(cum_sum(&in_seq), in_seq);
}

#[test]
fn test_single_element_cum_sum() {
    assert_eq!(cum_sum(&vec![1]), vec![1]);
}

#[test]
fn test_two_element_cum_sum() {
    assert_eq!(cum_sum(&vec![1, 2]), vec![1, 3]);
}

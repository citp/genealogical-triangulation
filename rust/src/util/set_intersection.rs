use std::arch::x86_64::{_mm_loadu_si128, __m128i, _mm_cmpeq_epi32, _mm_shuffle_epi32, _mm_or_si128, _mm_movemask_epi8};
use std::mem::transmute;

const CYCLIC_SHIFT: i32 = 57;

// adapted from https://highlyscalable.wordpress.com/2012/06/05/fast-intersection-sorted-lists-sse/
pub fn nonzero_intersection_simd(founders_a: &[u32], founders_b: &[u32]) -> bool {
    let mut i = 0;
    let mut j = 0;
    let a_newlen = (founders_a.len() / 4) * 4;
    let b_newlen = (founders_b.len() / 4) * 4;

    while i < a_newlen && j < b_newlen {
        unsafe {
            let a_ptr = founders_a[i..i+4].as_ptr();
            let v_a = _mm_loadu_si128(transmute::<*const u32, *const __m128i>(a_ptr));
            let b_ptr = founders_b[j..j+4].as_ptr();
            let mut v_b = _mm_loadu_si128(transmute::<*const u32, *const __m128i>(b_ptr));

            let cmp_mask1 = _mm_cmpeq_epi32(v_a, v_b);
            v_b = _mm_shuffle_epi32(v_b, CYCLIC_SHIFT);       // shuffling
            let cmp_mask2 = _mm_cmpeq_epi32(v_a, v_b);    // again...
            v_b = _mm_shuffle_epi32(v_b, CYCLIC_SHIFT);
            let cmp_mask3 = _mm_cmpeq_epi32(v_a, v_b);    // and again...
            v_b = _mm_shuffle_epi32(v_b, CYCLIC_SHIFT);
            let cmp_mask4 = _mm_cmpeq_epi32(v_a, v_b);
            let cmp_mask = _mm_or_si128(
                _mm_or_si128(cmp_mask1, cmp_mask2),
                _mm_or_si128(cmp_mask3, cmp_mask4)
            );
            let mask = _mm_movemask_epi8(cmp_mask);
            if mask != 0 {
                return true;
            }
        }

        if founders_a[i + 3] < founders_b[j + 3] {
            i += 4;
        } else {
            j += 4;
        }
    }

    return nonzero_intersection(&founders_a[i..], &founders_b[j..]);
}

pub fn nonzero_intersection(founders_a: &[u32], founders_b: &[u32]) -> bool {
    let mut i = 0;
    let mut j = 0;
    while i < founders_a.len() && j < founders_b.len() {
        if founders_a[i] < founders_b[j] {
            i += 1;
        } else if founders_b[j] < founders_a[i] {
            j += 1;
        } else {
            return true;
        }
    }

    false
}

#[cfg(test)]
mod test {
    use rand::prelude::*;
    use fnv::FnvHashSet;
    use std::iter::FromIterator;

    use util::set_intersection::{nonzero_intersection, nonzero_intersection_simd};

    // from https://users.rust-lang.org/t/how-to-serialize-a-u32-into-byte-array/986/5?u=paul-e
    fn transform_u32_to_array_of_u8(x:u32) -> [u8;4] {
        let b1 : u8 = ((x >> 24) & 0xff) as u8;
        let b2 : u8 = ((x >> 16) & 0xff) as u8;
        let b3 : u8 = ((x >> 8) & 0xff) as u8;
        let b4 : u8 = (x & 0xff) as u8;
        return [b1, b2, b3, b4]
    }

    fn get_founders(seed: u32) -> Vec<u32> {
        let mut seed_bytes = [0; 32];
        seed_bytes[0..4].copy_from_slice(&transform_u32_to_array_of_u8(seed));
        let mut rng: StdRng = SeedableRng::from_seed(seed_bytes);

        let mut founders = Vec::new();
        for _i in 0..255 {
            let num = rng.gen_range(0, 100000) as u32;
            founders.push(num);
        }
        founders
    }

    fn get_founder_pair_disjoint(seed: u32) -> (Vec<u32>, Vec<u32>) {
        //let mut rng = thread_rng();
        let mut seed_bytes = [0; 32];
        seed_bytes[0..4].copy_from_slice(&transform_u32_to_array_of_u8(seed));
        let mut rng: StdRng = SeedableRng::from_seed(seed_bytes);

        let mut founders = FnvHashSet::default();
        while founders.len() < 600 {
            let num = rng.gen_range(0, 100000) as u32;
            founders.insert(num);
        }
        let founder_vec: Vec<u32> = founders.iter().map(|x| *x).collect();
        (founder_vec[0..300].to_vec(), founder_vec[300..].to_vec())
    }

    fn get_founder_pair(seed: Option<u32>) -> (Vec<u32>, Vec<u32>) {
        let unwraped_seed = seed.unwrap_or(0);
        get_founder_pair_disjoint(unwraped_seed)
        //(get_founders(unwraped_seed), get_founders(unwraped_seed + 1))
        //((0..255).collect(), (255..500).collect())
    }


    #[test]
    fn test_same_output() {
        println!("Running same output.");
        for i in 0..1000 {
            println!("{}", i);
            let (mut founders_a, mut founders_b) = get_founder_pair(Some(i));
            founders_a.sort();
            founders_b.sort();

            let founders_a_set = FnvHashSet::from_iter(&founders_a);
            let founders_b_set = FnvHashSet::from_iter(&founders_b);
            assert!(nonzero_intersection(&founders_a, &founders_b) == founders_a_set.intersection(&founders_b_set).next().is_some());
        }
    }

    #[test]
    fn test_same_output_simd() {
        for i in 0..1000 {
            //println!("{}", i);
            let (mut founders_a, mut founders_b) = get_founder_pair(Some(i));
            founders_a.sort();
            founders_b.sort();

            let founders_a_set = FnvHashSet::from_iter(&founders_a);
            let founders_b_set = FnvHashSet::from_iter(&founders_b);
            let sorted_intersection = nonzero_intersection_simd(&founders_a, &founders_b);
            let set_intersection = founders_a_set.intersection(&founders_b_set).next().is_some();
            assert!(sorted_intersection == set_intersection, "Wrong intersection results on\n{:?}\n vs\n {:?}\n. Sorted intersection returned: {}\n", founders_a, founders_b, sorted_intersection);

            let a = vec![10, 20, 30, 40, 50, 60, 70, 80, 90];
            let b = vec![15, 22, 35, 45, 55, 65, 75, 85, 90];
            assert!(nonzero_intersection_simd(&a, &b));
        }
    }

    #[test]
    fn test_single_element_intersection_simd() {
        let a = vec![90];
        let b = vec![15, 25, 35, 45, 55, 65, 75, 85, 90];
        assert!(nonzero_intersection_simd(&a, &b));
        for i in 0..9 {
            assert!(nonzero_intersection_simd(&a, &b[1..]), "b slice: {:?}\n", &b[i..]);
        }

    }

    #[test]
    fn test_single_element_intersection() {
        let a = vec![90];
        let b = vec![15, 25, 35, 45, 55, 65, 75, 85, 90];
        assert!(nonzero_intersection(&a, &b));
        for i in 0..9 {
            assert!(nonzero_intersection(&a, &b[1..]), "b slice: {:?}\n", &b[i..]);
        }
    }
}
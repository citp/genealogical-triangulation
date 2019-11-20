extern crate rgsl;
extern crate csv;

use std::cell::RefCell;
use std::collections::HashMap;
use std::ops::{Add, AddAssign, DerefMut};
use std::path::Path;

use self::rgsl::randist::binomial;
use self::rgsl::RngType;
use rand::thread_rng;
// TODO: Fix this
#[allow(deprecated)]
use rand::distributions::{IndependentSample, Range};


use genome::ordered_float::OrderedFloat;

use genome::genome::{Genome, CHROMOSOMES};
use genome::diploid::Diploid;

pub static RECOMB_FILES: &'static [(u32, &'static str)] = &[(1, "genetic_map_chr1_b36.txt"), (2, "genetic_map_chr2_b36.txt"), (3, "genetic_map_chr3_b36.txt"), (4, "genetic_map_chr4_b36.txt"), (5, "genetic_map_chr5_b36.txt"), (6, "genetic_map_chr6_b36.txt"), (7, "genetic_map_chr7_b36.txt"), (8, "genetic_map_chr8_b36.txt"), (9, "genetic_map_chr9_b36.txt"), (10, "genetic_map_chr10_b36.txt"), (11, "genetic_map_chr11_b36.txt"), (12, "genetic_map_chr12_b36.txt"), (13, "genetic_map_chr13_b36.txt"), (14, "genetic_map_chr14_b36.txt"), (15, "genetic_map_chr15_b36.txt"), (16, "genetic_map_chr16_b36.txt"), (17, "genetic_map_chr17_b36.txt"), (18, "genetic_map_chr18_b36.txt"), (19, "genetic_map_chr19_b36.txt"), (20, "genetic_map_chr20_b36.txt"), (21, "genetic_map_chr21_b36.txt"), (22, "genetic_map_chr22_b36.txt")];

pub struct RecombinatorPair {
    pub male: Recombinator,
    pub female: Recombinator
}

pub struct Recombinator {
    pub chromosomes: HashMap<u32, ChromRecombData>,
    rng: RefCell<rgsl::Rng>
}

pub struct ChromRecombData {
    pub num_bases: u32,
    pub num_centimorgans: f64,
    pub end_points: Vec<OrderedFloat<f64>>,
    pub end_point_range: HashMap<OrderedFloat<f64>, (u32, u32)>,
    pub chrom_start_offset: u32
}

impl Recombinator {
    pub fn new(recombination_data: &HashMap<u32, Vec<(u32, f64, f64)>>)
               -> Recombinator {
        let mut chromosomes: HashMap<u32, ChromRecombData> = HashMap::new();
        let mut cum_bases = 0;
        for chromosome in CHROMOSOMES {
            let data = &recombination_data[chromosome];
            let last_line = data.last().unwrap();
            // println!("Chrom: {}, last_line: {:?}", chromosome, last_line);
            let num_bases = last_line.0;
            let num_centimorgans = last_line.2;
            let chrom_end_points: Vec<OrderedFloat<f64>> = data[1..].iter()
                .map(|row| OrderedFloat::from(row.2)).collect();
            let mut end_point_range: HashMap<OrderedFloat<f64>,
                                             (u32, u32)> = HashMap::new();
            for i in 0..data.len() - 1 {
                let pos_1 = data[i].0;
                let pos_2 = data[i + 1].0;
                let end_point: OrderedFloat<f64>
                    = OrderedFloat::from(data[i + 1].2);
                assert!(end_point == chrom_end_points[i]);
                end_point_range.insert(end_point, (pos_1, pos_2));
            }
            assert!(end_point_range.len() == chrom_end_points.len());
            for (i, end_point) in chrom_end_points.iter().enumerate() {
                if !end_point_range.contains_key(&end_point) {
                    println!("Missing end point: {:.32} for chrom {} at line {} of {}",
                             end_point, chromosome, i, chrom_end_points.len() - 1);
                }
                assert!(end_point_range.contains_key(&end_point));
            }
            let chrom_data = ChromRecombData{ num_bases: num_bases,
                                              num_centimorgans: num_centimorgans,
                                              end_points: chrom_end_points,
                                              end_point_range: end_point_range,
                                              chrom_start_offset: cum_bases };
            cum_bases += num_bases;
            chromosomes.insert(*chromosome, chrom_data);
        }
        rgsl::RngType::env_setup();
        let t: RngType = rgsl::rng::default();
        let r = rgsl::Rng::new(&t).unwrap();
        Recombinator { chromosomes: chromosomes, rng: RefCell::new(r) }
    }

    ///
    /// Returns recombination locations for the given chromosome.
    fn _recombination_locations(&self, chromosome_num: u32) -> Vec<u32> {
        let chromosome = self.chromosomes.get(&chromosome_num).unwrap();
        let p = (chromosome.num_centimorgans * 0.01)
            / (chromosome.num_bases as f64);
        let num_events = binomial::binomial(self.rng.borrow_mut().deref_mut(),
                                            p, chromosome.num_bases);
        let mut loci: Vec<u32> = Vec::with_capacity(num_events as usize);
        if num_events == 0 {
            return loci;
        }
        let mut rng = thread_rng();
        let loc_range: Range<f64> = Range::new(0.0,
                                               chromosome.num_centimorgans);
        let mut locations = Vec::with_capacity(num_events as usize);
        for _ in 0..num_events {
            // TODO: Fix this
            #[allow(deprecated)]
            locations.push(OrderedFloat::from(loc_range.ind_sample(&mut rng)));
        }
        locations.sort();
        // println!("Chrom locations: {:?}", locations);
        for location in locations {
            let index = match chromosome.end_points.binary_search(&location) {
                Ok(i) => i,
                Err(i) => i
            };
            let end_point = chromosome.end_points[index];
            let range = chromosome.end_point_range.get(&end_point).unwrap();
            let (start, stop) = (range.0, range.1);
            let start_point = match index {
                0 => OrderedFloat::from(0.0),
                _ => chromosome.end_points[index - 1]
            };
            assert!(location > start_point);
            let f64_location: f64 = location.into();
            let f64_start_point: f64 = start_point.into();
            let f64_end_point: f64 = end_point.into();
            let fraction_in: f64 = (f64_location - f64_start_point) /
                (f64_end_point - f64_start_point);
            let recomb_spot = ((stop - start) as f64 * fraction_in
                               + start as f64) as u32;
            if loci.len() == 0 || recomb_spot != loci[loci.len() - 1] {
                loci.push(recomb_spot);
            }
        }
        // println!("Chrom locations: {:?}", loci);
        loci
    }

    pub fn recombination(&self, genome: &Genome) -> Genome {
        let mut global_locations: Vec<u32> = Vec::with_capacity(30);
        for chrom_num in CHROMOSOMES {
            let chromosome = self.chromosomes.get(chrom_num).unwrap();
            let mut locations = self._recombination_locations(*chrom_num);
            if locations.len() % 2 == 1 {
                assert!(locations[locations.len() - 1] <= chromosome.num_bases);
                locations.push(chromosome.num_bases);
            }
            let offset = chromosome.chrom_start_offset;
            // println!("Chrom num: {}, Chrom offset: {}", chrom_num, offset);
            for location in locations {
                global_locations.push(location + offset);
            }
        }
        // println!("Swapping at locations {:?}", global_locations);
        swap_at_locations(genome, &global_locations)
    }
                     
}

pub fn recombinators_from_directory(directory: &str) -> RecombinatorPair {
    /* Given a directory of files downloaded from, returns a Recombinator
    object for those files.
    Genome maps are from
    http://hapmap.ncbi.nlm.nih.gov/downloads/recombination/
    Sex based centimorgan lengths are from decode doi:10.1038/ng917. */
    let mut files: HashMap<u32, String> = HashMap::new();
    for &(chrom, filename) in RECOMB_FILES {
        let full_path = Path::new(directory).join(filename);
        files.insert(chrom, String::from(full_path.to_str().unwrap()));
    }
    let sex_lengths = get_sex_lengths();
    recombinators_from_hapmap_files(&files, sex_lengths)
}

fn recombinators_from_hapmap_files(hapmap_files: &HashMap<u32, String>,
                                   sex_lengths: SexCmLengths)
                                   -> RecombinatorPair {
    let mut chrom_data = HashMap::new();
    for (chromosome, filename) in hapmap_files {
        chrom_data.insert(*chromosome, read_recombination_file(filename));
    }

    let male_adjusted_data = adjust_chromosomes(&chrom_data,
                                                    &sex_lengths.male);
    let female_adjusted_data = adjust_chromosomes(&chrom_data,
                                                  &sex_lengths.female);
    let male_recombinator = Recombinator::new(&male_adjusted_data);
    let female_recombinator = Recombinator::new(&female_adjusted_data);
    RecombinatorPair {
        male: male_recombinator,
        female: female_recombinator
    }
}

fn adjust_chromosomes(chrom_data: &HashMap<u32, Vec<(u32, f64, f64)>>,
                      sex_lengths: &HashMap<u32, f64>)
                      -> HashMap<u32, Vec<(u32, f64, f64)>> {
    let mut adjusted_data = HashMap::new();
    for (chrom, length) in sex_lengths {
        let data = &chrom_data[chrom];
        let original_length = data[data.len() - 1].2;
        let ratio: f64 = length / original_length;
        adjusted_data.insert(*chrom, adjust_centimorgans(data, ratio));
    }
    adjusted_data
}

fn adjust_centimorgans(rows: &Vec<(u32, f64, f64)>, multiplier: f64) ->
    Vec<(u32, f64, f64)> {
    /* Adjusts the number of centimorgans in rows by some multiplier.
    eg if we have the rows

    a 0.05 0
    b 0.1  0.05
    c 0.05 0.15

    and we apply a multiplier of 2, then we are doubling the number of
    centimorgans, giving us.

    a 0.1 0
    b 0.2 0.1
    c 0.1 0.3 */
    rows.iter().map(|&(pos, rate, distance)|
                    (pos, rate * multiplier, distance * multiplier))
               .collect()
}

pub fn read_recombination_file(filename: &str) -> Vec<(u32, f64, f64)> {
    /*
    Reads a recombination file and returns the rows, with columns
    converted to the appropriate numeric types.
     */
    let mut ret = Vec::new();
    let mut rdr = csv::Reader::from_file(filename).unwrap()
        .has_headers(true)
        .delimiter(b' ');
    for row in rdr.decode() {
        let (bp, rate, distance): (u32, f64, f64) = row.unwrap();
        ret.push((bp, rate, distance));
    }
    ret
}

// Stores the lengths of the chromosomes (in CM) for each sex
struct SexCmLengths {
    pub male: HashMap<u32, f64>,
    pub female: HashMap<u32, f64>
}

fn list_to_map(list: &[(u32, f64)]) -> HashMap<u32, f64> {
    let mut map = HashMap::with_capacity(list.len());
    for &(chrom, length) in list {
        map.insert(chrom, length);
    }
    map
}
fn get_sex_lengths() -> SexCmLengths {
    // Hard coded hack.
    // TODO: read this from the file in the future.
    let male_list = [(1, 195.12), (2, 189.55), (3, 160.71), (4,
    146.54), (5, 151.2), (6, 137.62), (7, 128.35), (8, 107.94), (9,
    117.25), (10, 133.89), (11, 109.36), (12, 135.54), (13, 101.31),
    (14, 94.62), (15, 102.57), (16, 108.1), (17, 108.56), (18, 98.62),
    (19, 92.64), (20, 74.72), (21, 47.31), (22, 48.96)];

    let female_list = [(1, 345.41), (2, 325.41), (3, 275.64), (4,
    259.06), (5, 260.19), (6, 241.59), (7, 230.33), (8, 209.94), (9,
    198.2), (10, 218.13), (11, 195.53), (12, 206.64), (13, 155.88),
    (14, 142.36), (15, 154.96), (16, 149.62), (17, 161.53), (18,
    142.57), (19, 126.82), (20, 121.97), (21, 76.4), (22, 82.76)];

    SexCmLengths {male: list_to_map(&male_list),
                  female: list_to_map(&female_list)}
}


fn search_left(to_search: &Vec<u32>, find: u32) -> usize {
    let search = to_search.binary_search(&find);
    match search {
        Ok(mut i) => {
            while i > 0 && to_search[i - 1] == find {
                i -= 1;
            }
            i
        },
        Err(i) => i
    }
}

pub fn cum_sum<T: Add + Copy + AddAssign>(in_seq: &[T]) -> Vec<T> {
    let mut ret = Vec::with_capacity(in_seq.len());
    if in_seq.len() == 0 {
        return ret;
    }
    ret.push(in_seq[0]);
    let mut accum = in_seq[0];
    for item in &in_seq[1..] {
        accum += *item;
        ret.push(accum);
    }
    ret
}

pub fn swap_at_locations(genome: &Genome, locations: &Vec<u32>) -> Genome {
    // println!("swap_at_locations");
    let temp_mother = new_sequence(&genome.mother, locations);
    let temp_father = new_sequence(&genome.father, locations);
    let mut new_father = Vec::with_capacity(temp_father.starts.len());
    let mut new_mother = Vec::with_capacity(temp_mother.starts.len());
    let mut new_father_founder = Vec::with_capacity(new_father.len());
    let mut new_mother_founder = Vec::with_capacity(new_mother.len());
    let mut mother_prev = 0;
    let mut father_prev = 0;
    for pair in locations.chunks(2) {
        // println!("{:?}", pair);
        let start = pair[0];
        let stop = pair[1];
        // println!("temp_mother.starts: {:?}", &temp_mother.starts);
        let mother_start_i = search_left(&temp_mother.starts, start);
        // println!("mother_start_i: {}", mother_start_i);
        let mother_stop_i = search_left(&temp_mother.starts, stop);
        // println!("mother_stop_i: {}", mother_stop_i);
        let father_start_i = search_left(&temp_father.starts, start);
        let father_stop_i = search_left(&temp_father.starts, stop);
        // println!("To extend: {:?}", &temp_mother.starts[mother_prev..mother_start_i]);
        new_mother.extend_from_slice(&temp_mother.starts[mother_prev..mother_start_i]);
        new_mother.extend_from_slice(&temp_father.starts[father_start_i..father_stop_i]);
        new_mother_founder.extend_from_slice(&temp_mother.founder[mother_prev..mother_start_i]);
        new_mother_founder.extend_from_slice(&temp_father.founder[father_start_i..father_stop_i]);
        mother_prev = mother_stop_i;


        new_father.extend_from_slice(&temp_father.starts[father_prev..father_start_i]);
        new_father.extend_from_slice(&temp_mother.starts[mother_start_i..mother_stop_i]);
        new_father_founder.extend_from_slice(&temp_father.founder[father_prev..father_start_i]);
        new_father_founder.extend_from_slice(&temp_mother.founder[mother_start_i..mother_stop_i]);
        father_prev = father_stop_i;
    }
    new_mother.extend_from_slice(&temp_mother.starts[mother_prev..]);
    new_mother_founder.extend_from_slice(&temp_mother.founder[mother_prev..]);
    new_father.extend_from_slice(&temp_father.starts[father_prev..]);
    new_father_founder.extend_from_slice(&temp_father.founder[father_prev..]);
    let mother = Diploid {starts: new_mother, founder: new_mother_founder,
                          end: temp_mother.end};
    let father = Diploid {starts: new_father, founder: new_father_founder,
                          end: temp_father.end};
    Genome {mother: mother, father: father}
}

pub fn new_sequence(diploid: &Diploid, locations: &Vec<u32>) -> Diploid {
    if locations.len() == 0 {
        return diploid.clone();
    }
    let max_locations = locations.len() + diploid.starts.len();
    let mut new_starts: Vec<u32> = Vec::with_capacity(max_locations);
    let mut new_founder: Vec<u32> = Vec::with_capacity(max_locations);
    let trunc_locations;
    if locations[locations.len() - 1] == diploid.end {
        trunc_locations = &locations[0..locations.len() - 1];
    } else {
        trunc_locations = locations;
    }
    let mut locations_i = 0;
    let mut starts_i = 0;
    // let mut founder_i = 0;
    assert!(diploid.starts[0] <= locations[0]);
    for _ in 0..(max_locations) {
        if trunc_locations.len() <= locations_i {
            new_starts.extend_from_slice(&diploid.starts[starts_i..]);
            new_founder.extend_from_slice(&diploid.founder[starts_i..]);
            break;
        }
        if diploid.starts.len() <= starts_i {
            new_starts.extend_from_slice(&trunc_locations[locations_i..]);
            let founder = new_founder[new_founder.len() - 1];
            for _ in 0..(trunc_locations[locations_i..].len() as u32) {
                new_founder.push(founder);
            }
            break;
        }
        let diploid_loci = diploid.starts[starts_i];
        let break_loci = trunc_locations[locations_i];
        if diploid_loci < break_loci {
            new_starts.push(diploid_loci);
            new_founder.push(diploid.founder[starts_i]);
            starts_i += 1;
        } else if diploid_loci > break_loci {
            new_starts.push(break_loci);
            new_founder.push(diploid.founder[starts_i - 1]);
            locations_i += 1;
        } else { // diploid_loci == break_loci
            new_starts.push(diploid_loci);
            new_founder.push(diploid.founder[starts_i]);
            starts_i += 1;
            locations_i += 1;
        }
    }
    Diploid { starts: new_starts, founder: new_founder, end: diploid.end  }
}

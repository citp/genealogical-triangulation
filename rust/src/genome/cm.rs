use std::collections::HashMap;
use std::path::Path;

use genome::recombinator::{RECOMB_FILES, read_recombination_file};


pub struct CmConverter{
    bases: Vec<u32>,
    cm: Vec<f64>,
    rates: Vec<f64>
}

impl CmConverter{
    pub fn new(recomb_dir: &str) -> Self {
        let mut recomb_data = HashMap::new();
        for &(chrom, filename) in RECOMB_FILES {
            let full_path = Path::new(recomb_dir).join(filename);
            let path_str = full_path.to_str().unwrap();
            let data = read_recombination_file(&path_str);
            recomb_data.insert(chrom, data);
        }
        let mut chroms: Vec<_> = RECOMB_FILES.iter().map(|x| x.0).collect();
        chroms.sort_unstable();
        let mut bases = Vec::new();
        let mut cm = Vec::new();
        let mut rates = Vec::new();
        let mut bases_accum = 0;
        let mut cm_accum = 0.0;
        for chrom in chroms {
            let data = recomb_data.get(&chrom).unwrap();
            for (loci, rate, cumulative_cm) in data {
                bases.push(loci + bases_accum);
                cm.push(cumulative_cm + cm_accum);
                rates.push(rate / 1000000.0);
            }
            let last_entry = data.last().unwrap();
            bases_accum += last_entry.0;
            cm_accum += last_entry.2;
        }
        CmConverter {bases, cm, rates}
    }

    // Locations must be sorted.
    fn cumulative_cm(&self, locations: &[u32]) -> Vec<f64> {
        //println!("Cumulative cm called");
        let last_base = self.bases.last().unwrap();
        let last_cm = *self.cm.last().unwrap();
        let mut ret = Vec::with_capacity(locations.len());
        for location in locations {
            assert!(location <= last_base);
            let current_index = match self.bases.binary_search(location) {
                Ok(i) => i,
                Err(i) => i
            };
            //println!("current_index: {}, location: {}", current_index, location);
            let cm_distance = self.cm[current_index];
            let bp_difference = self.bases[current_index] - location;
            let cm_difference = if current_index == 0 {
                 0.0
            } else {
                bp_difference as f64 * self.rates[current_index - 1]
            };
            let adjusted_cm_distance = cm_distance - cm_difference;
            assert!(0.0 <= adjusted_cm_distance,
                    "current_index:{}, \nlocation: {}, recomb_location: {}, recomb_location - 1: {}\nadjusted distance: {}, cm_distance: {}, cm_difference: {}\nbp_difference: {},",
                    current_index, location, self.bases[current_index], self.bases[current_index - 1], adjusted_cm_distance, cm_distance, cm_difference, bp_difference);
            assert!(adjusted_cm_distance <= last_cm);
            ret.push(adjusted_cm_distance);
        }
        ret
    }

    pub fn cm_lengths(&self, regions: &[(u32, u32)]) -> Vec<f64> {
        let mut starts = Vec::with_capacity(regions.len());
        let mut stops = Vec::with_capacity(regions.len());
        for region in regions {
            starts.push(region.0);
            stops.push(region.1);
        }
        let cm_starts = self.cumulative_cm(&starts);
        let cm_stops = self.cumulative_cm(&stops);
        cm_starts.iter().zip(cm_stops.iter()).map(|(start, stop)| stop - start).collect()
    }
}
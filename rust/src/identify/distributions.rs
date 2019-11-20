// https://docs.rs/GSL/0.4.27/rgsl/randist/gamma/index.html
// https://www.gnu.org/software/gsl/manual/html_node/The-Gamma-Distribution.html
use rgsl::randist::gamma::gamma_pdf;
use fnv::FnvHashMap;

use population::node_pair::NodeIdPair;
use statistics::gamma::HurdleGammaParams;

#[derive(Copy, Clone, Serialize, Deserialize, Debug)]
struct HurdleGammaParamsWithLog {
    gamma: HurdleGammaParams,
    log_zero_prob: f64
}

impl From<HurdleGammaParams> for HurdleGammaParamsWithLog {
    fn from(other: HurdleGammaParams) -> Self {
        Self {
            gamma: other,
            log_zero_prob: other.zero_prob.ln()
        }
    }
}

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct Distribution {
    distributions: FnvHashMap<NodeIdPair, HurdleGammaParamsWithLog>,
    cryptic_distribution: HurdleGammaParamsWithLog,
    pub labeled_nodes: Vec<u32>
}


impl Distribution {

    pub fn new(mut distributions: FnvHashMap<NodeIdPair, HurdleGammaParams>, cryptic_distribution: HurdleGammaParams, labeled_nodes: &[u32]) -> Distribution {
        let mut new_distributions = FnvHashMap::default();
        for (key, value) in distributions.drain() {
            new_distributions.insert(key, value.into());
        }
        Distribution {
            distributions: new_distributions,
            cryptic_distribution: cryptic_distribution.into(),
            labeled_nodes: labeled_nodes.to_vec(),
        }
    }
    pub fn get_probability(&self, length: u64, pair: &NodeIdPair) -> f64 {
        match self.distributions.get(pair) {
            Some(params) => gamma_prob(length, params, false),
            None => gamma_prob(length, &self.cryptic_distribution, false)
        }
    }

    pub fn get_log_probability(&self, length: f64, pair: &NodeIdPair) -> f64 {
        match self.distributions.get(pair) {
            Some(params) => gamma_log_prob(length, params),
            None => gamma_log_prob(length, &self.cryptic_distribution)
        }
    }

}

pub struct DistributionAlt {
    distributions: FnvHashMap<u32, Vec<(u32, HurdleGammaParamsWithLog)>>,
    cryptic_distribution: HurdleGammaParamsWithLog,
    pub labeled_nodes: Vec<u32>
}

impl From<&Distribution> for DistributionAlt {
    fn from(other: &Distribution) -> Self {
        let mut dist = FnvHashMap::default();
        for (key, value) in other.distributions.iter() {
            let mut list = dist.entry(key.unlabeled).or_insert_with(Vec::new);
            list.push((key.labeled, *value));
        }

        for value in dist.values_mut() {
            value.sort_unstable_by_key(|x| x.0);
        }

        Self {
            distributions: dist,
            cryptic_distribution: other.cryptic_distribution.clone(),
            labeled_nodes: other.labeled_nodes.clone()
        }
    }
}

impl DistributionAlt {
    pub fn get_log_probabilities(&self, unlabeled: u32, lengths: &[(u32, f64)]) -> Vec<f64> {
        let empty = Vec::new();
        let distributions = self.distributions.get(&unlabeled).unwrap_or(&empty);
        let mut dist_i = 0;
        let mut ret = Vec::with_capacity(lengths.len());
        for (labeled, length) in lengths.iter() {
            if dist_i < distributions.len() && *labeled == distributions[dist_i].0 {
                ret.push(gamma_log_prob(*length, &distributions[dist_i].1));
                dist_i += 1;
            } else {
                ret.push(gamma_log_prob(*length, &self.cryptic_distribution));
            }
        }
        assert_eq!(dist_i, distributions.len());
//        match self.distributions.get(pair) {
//            Some(params) => gamma_log_prob(length, params),
//            None => gamma_log_prob(length, &self.cryptic_distribution)
//        }
        ret
    }
}

static mut LENGTH: f64 = 0.0;
static mut SHAPE: f64 = 0.0;
static mut SCALE: f64 = 0.0;


fn gamma_prob(length: u64, params: &HurdleGammaParamsWithLog, log_prob: bool) -> f64 {
    if length == 0 {
        if log_prob {
            params.log_zero_prob
        } else {
            params.gamma.zero_prob
        }
    } else {
//        unsafe {
//            LENGTH = length as f64;
//            SHAPE = params.shape;
//            SCALE = params.scale;
//        }
        let gamma_params = params.gamma;
        let gamma_pdf = gamma_pdf(length as f64, gamma_params.shape, gamma_params.scale);
        let mut gamma_prob = gamma_pdf * (1.0 - gamma_params.zero_prob);
        if gamma_prob <= 0.0 {
            gamma_prob = 1e-12;
        }
        if log_prob {
            gamma_prob.ln()
        } else {
            gamma_prob
        }

    }
}

fn gamma_log_prob(length: f64, params: &HurdleGammaParamsWithLog) -> f64 {
    if length == 0.0 {
        params.log_zero_prob
    } else {
        let gamma_params = params.gamma;
        let gamma_pdf = gamma_pdf(length, gamma_params.shape, gamma_params.scale);
        let mut gamma_prob = gamma_pdf * (1.0 - gamma_params.zero_prob);
        if gamma_prob <= 0.0 {
            gamma_prob = 1e-12;
        }
        gamma_prob.ln()
    }
}


pub fn get_length() -> f64 {
    unsafe {
        LENGTH
    }
}


pub fn get_shape() -> f64 {
    unsafe {
        SHAPE
    }
}

pub fn get_scale() -> f64 {
    unsafe {
        SCALE
    }
}

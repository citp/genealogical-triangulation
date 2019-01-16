// https://docs.rs/GSL/0.4.27/rgsl/randist/gamma/index.html
// https://www.gnu.org/software/gsl/manual/html_node/The-Gamma-Distribution.html
use rgsl::randist::gamma::gamma_pdf;
use fnv::FnvHashMap;

use population::node_pair::NodeIdPair;
use statistics::gamma::HurdleGammaParams;

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct Distribution {
    pub distributions: FnvHashMap<NodeIdPair, HurdleGammaParams>,
    pub cryptic_distribution: HurdleGammaParams,
    pub labeled_nodes: Vec<u32>,
}

impl Distribution {
    pub fn get_probability(&self, length: u64, pair: &NodeIdPair) -> f64 {
        match self.distributions.get(pair) {
            Some(params) => gamma_prob(length, params),
            None => gamma_prob(length, &self.cryptic_distribution),
        }
    }
}

fn gamma_prob(length: u64, params: &HurdleGammaParams) -> f64 {
    if length == 0 {
        params.zero_prob
    } else {
        let gamma_pdf = gamma_pdf(length as f64, params.shape, params.scale);
        let mut gamma_prob = gamma_pdf * (1.0 - params.zero_prob);
        if gamma_prob <= 0.0 {
            gamma_prob = 1e-12;
        }
        gamma_prob
    }
}

fn smoothing_prob(length: u64) -> f64 {
    if length == 0 {
        0.98162761
    } else if length < 30000000 {
        0.01963307
    } else {
        0.00117543
    }
}

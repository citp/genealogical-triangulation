use rand;
use rgsl::psi::diagamma;
use rgsl::psi::trigamma;

use rand::distributions::{Distribution, Uniform};

const SUFFICIENT_DATA_POINTS: usize = 5;

#[derive(PartialEq, Debug, Clone, Copy)]
pub struct GammaParams {
    shape: f64,
    scale: f64
}

#[derive(PartialEq, Clone, Copy, Serialize, Deserialize, Debug)]
pub struct HurdleGammaParams {
    pub zero_prob: f64,
    pub shape: f64,
    pub scale: f64
}

impl HurdleGammaParams {
    pub fn has_nan_parameters(&self) -> bool {
        self.zero_prob.is_nan() || self.shape.is_nan() || self.scale.is_nan()
    }
}


pub fn fit_gamma(data: &[f64]) -> GammaParams {
    let mut sum: f64 = 0.0;
    let mut sum_of_log: f64 = 0.0;
    let noise = Uniform::new(1e-8, 10000.0);
    let mut rng = rand::thread_rng();
    for val in data {
        let noisy_data = *val + noise.sample(&mut rng);
        sum += noisy_data;
        sum_of_log += noisy_data.ln();
    }
    let data_mean = sum / (data.len() as f64);
    let mean_of_logs = sum_of_log / (data.len() as f64);
    let log_of_mean = (data_mean).ln();
    let log_diff = mean_of_logs - log_of_mean;
    let mut shape = 0.5 / (log_of_mean - mean_of_logs);
    let mut shape_reciprocal = 1.0 / shape;
    let mut difference = 1.0;
    while difference > 0.000005 {
        let numerator = log_diff + shape.ln() - diagamma::psi(shape);
        let denominator = shape.powi(2) *
            (shape_reciprocal - trigamma::psi_1(shape));
        let tmp_shape_reciprocal = shape_reciprocal + (numerator / denominator);
        let tmp_shape = 1.0 / tmp_shape_reciprocal;
        difference = (tmp_shape - shape).abs();
        shape = tmp_shape;
        shape_reciprocal = tmp_shape_reciprocal;
    }
    GammaParams { shape: shape, scale: data_mean / shape }
}

pub fn fit_hurdle_gamma(data: &[u64]) -> Option<HurdleGammaParams> {
    let nonzero: Vec<_> = data.iter().filter(|&x| *x != 0)
        .map(|x| *x as f64)
        .collect();
    if nonzero.len() < SUFFICIENT_DATA_POINTS {
        return None;
    }
    let zero_prob = (data.len() - nonzero.len()) as f64 / data.len() as f64;
    let gamma = fit_gamma(&nonzero);
    if gamma.scale.is_nan() || gamma.shape.is_nan() {
        return None;
    }
    Some(HurdleGammaParams {zero_prob: zero_prob,
                       shape: gamma.shape,
                       scale: gamma.scale})
}

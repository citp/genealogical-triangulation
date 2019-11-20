# Python's built in log is faster for non vector operations
from math import log, isnan
from warnings import warn

import numpy as np
import numpy.ma as ma
from scipy.special import digamma, polygamma

SUFFICIENT_DATA_POINTS = 5


def fit_gamma(data):
    """
    Gamma parameter estimation using the algorithm described in
    http://research.microsoft.com/en-us/um/people/minka/papers/minka-gamma.pdf
    Returns a tuple with (shape, scale) parameters.
    """
    # data += 1e-8 # Add small number to avoid 0s in the data causing issues.
    # Add small amount of noise to avoid 0s in the data causing issues
    # or all values being identical causing issues.
    data += np.random.uniform(1e-8, 10000, len(data))
    data_mean = np.mean(data)
    log_of_mean = log(data_mean)
    mean_of_logs = np.mean(np.log(data))
    log_diff = mean_of_logs - log_of_mean
    shape = 0.5 / (log_of_mean - mean_of_logs)
    shape_reciprocal = 1  / shape
    difference = 1
    while difference > 0.000005:
        numerator = log_diff + log(shape) - digamma(shape)
        denominator = (shape ** 2) * (shape_reciprocal - polygamma(1, shape))
        tmp_shape_reciprocal = shape_reciprocal + numerator / denominator
        tmp_shape = 1 / tmp_shape_reciprocal
        difference = abs(tmp_shape - shape)
        shape = tmp_shape
        shape_reciprocal = tmp_shape_reciprocal
    return (shape, data_mean / shape)

def fit_hurdle_gamma(data):
    data = np.asarray(data, dtype = np.float64)
    nonzero = data != 0.0
    num_nonzero = np.sum(nonzero)
    num_zero = len(data) - num_nonzero
    if num_nonzero <= SUFFICIENT_DATA_POINTS:
        # Handle this case better
        return (None, None, None)
    prob_zero = num_zero / len(data)
    data = data[nonzero]
    # data += 1e-8 # Add small number to avoid 0s in the data causing issues.
    # Add small amount of noise to avoid 0s in the data causing issues
    # or all values being identical causing issues.
    data += np.random.uniform(1e-8, 0.2, len(data))
    data_mean = np.mean(data)
    log_of_mean = log(data_mean)
    mean_of_logs = np.mean(np.log(data))
    log_diff = mean_of_logs - log_of_mean
    shape = 0.5 / (log_of_mean - mean_of_logs)
    shape_reciprocal = 1  / shape
    difference = 1
    while difference > 0.000005:
        numerator = log_diff + log(shape) - digamma(shape)
        denominator = (shape ** 2) * (shape_reciprocal - polygamma(1, shape))
        tmp_shape_reciprocal = shape_reciprocal + numerator / denominator
        tmp_shape = 1 / tmp_shape_reciprocal
        difference = abs(tmp_shape - shape)
        shape = tmp_shape
        shape_reciprocal = tmp_shape_reciprocal
    scale = data_mean / shape
    if isnan(shape) or isnan(scale):
        warn("NaN shape or scale value")
    return (shape, scale, prob_zero)

def fit_hurdle_gamma_vector(data):
    nonzero = data != 0.0
    num_nonzero = np.sum(nonzero, axis = 1, keepdims = True)
    num_zero = data.shape[1] - num_nonzero
    insufficient_data = (num_nonzero <= SUFFICIENT_DATA_POINTS).flatten()
    prob_zero = num_zero / data.shape[1]
    data = ma.array(data = data, mask = ~nonzero)
    # Add small number to avoid 0s in the data causing issues.
    # Add small amount of noise to avoid 0s in the data causing issues
    # or all values being identical causing issues.
    data += np.random.uniform(1e-8, 0.2, size = data.shape)
    data_mean = np.mean(data, axis = 1)
    log_of_mean = np.log(data_mean)
    mean_of_logs = np.mean(ma.log(data), axis = 1)
    log_diff = mean_of_logs - log_of_mean
    shape = 0.5 / (log_of_mean - mean_of_logs)
    shape_reciprocal = 1  / shape
    difference = 1
    while difference > 0.000005:
        numerator = log_diff + np.log(shape) - digamma(shape)
        denominator = (shape ** 2) * (shape_reciprocal - polygamma(1, shape))
        tmp_shape_reciprocal = shape_reciprocal + numerator / denominator
        tmp_shape = 1 / tmp_shape_reciprocal
        difference = np.max(np.abs(tmp_shape - shape))
        shape = tmp_shape
        shape_reciprocal = tmp_shape_reciprocal
    scale = data_mean / shape
    if np.any(np.isnan(shape)) or np.any(np.isnan(scale)):
        warn("NaN shape or scale value")
    return (shape.data, scale.data, prob_zero.flatten(), insufficient_data)

# Python's built in log is faster for non vector operations
from math import log, isnan
from warnings import warn

import numpy as np
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
    nonzero = data != 0
    num_nonzero = np.sum(nonzero)
    num_zero = len(data) - num_nonzero
    if num_nonzero <= SUFFICIENT_DATA_POINTS:
        # Handle this case better
        return (None, None, None)
    prob_zero = num_zero / len(data)
    data = np.array(data[nonzero], dtype = np.float64)
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
    scale = data_mean / shape
    if isnan(shape) or isnan(scale):
        warn("NaN shape or scale value")
    return (shape, scale, prob_zero)

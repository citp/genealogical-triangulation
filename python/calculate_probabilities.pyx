# cython: language_level=3

import numpy as np
from array import array

cimport numpy as np
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef calculate_probabilities(np.ndarray[np.uint32_t, ndim=1] dist_labeled,
                              np.ndarray[np.float64_t, ndim=1] shape,
                              np.ndarray[np.float64_t, ndim=1] scale,
                              np.ndarray[np.float64_t, ndim=1] zero_prob,
                              np.ndarray[np.float64_t, ndim=1] lens,
                              np.ndarray[np.uint32_t, ndim=1] labeled,
                              np.ndarray[np.float64_t, ndim=1] cryptic_prob,
                              unsigned long cur_node):
  assert len(shape) == len(scale) == len(zero_prob) == len(dist_labeled)
  assert len(cryptic_prob) == len(lens) == len(labeled)
  cpdef unsigned long i = 0
  cpdef unsigned long dist_i = 0
  cdef unsigned long labeled_len = len(dist_labeled)
  cpdef unsigned long labeled_node
  cdef list ret_shapes = []
  cdef list ret_scales = []
  cdef list ret_zero_prob = []
  cdef list ret_lengths = []
  cdef double ret_cryptic_prob = 0.0
  for i in range(len(labeled)):
    labeled_node = labeled[i]
    if labeled_node == cur_node:
        continue
    while dist_i < labeled_len and dist_labeled[dist_i] < labeled_node:
        dist_i += 1
    if dist_i < labeled_len and dist_labeled[dist_i] == labeled_node:
        ret_shapes.append(shape[dist_i])
        ret_scales.append(scale[dist_i])
        ret_zero_prob.append(zero_prob[dist_i])
        ret_lengths.append(lens[i])
        dist_i += 1
    else:
        ret_cryptic_prob += cryptic_prob[i]
    

  return (ret_lengths, ret_shapes, ret_scales, ret_zero_prob, ret_cryptic_prob)

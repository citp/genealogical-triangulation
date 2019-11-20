# cython: language_level=3

import numpy as np
from array import array

cimport numpy as np
cimport cython

cpdef unpack_params(params):
  cdef unsigned long params_len = len(params)
  cdef np.ndarray[np.float64_t, ndim=1] shapes = np.empty(params_len,
                                                          dtype = np.float64)
  cdef np.ndarray[np.float64_t, ndim=1] scales = np.empty(params_len,
                                                          dtype = np.float64)
  cdef np.ndarray[np.float64_t, ndim=1] zero_probs = np.empty(params_len,
                                                              dtype = np.float64)
  cdef unsigned long i = 0
  for i in range(params_len):
      shape, scale, zero_prob = params[i]
      shapes[i] = shape
      scales[i] = scale
      zero_probs[i] = zero_prob
  return (shapes, scales, zero_probs)


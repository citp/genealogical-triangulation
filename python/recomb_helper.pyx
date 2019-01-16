# cython: profile=True

from diploid import Diploid

import numpy as np

cimport numpy as np
cimport cython

def new_sequence(diploid, np.ndarray[np.uint32_t, ndim=1] locations):
    """
    Return a new sequence, broken up at the given start, stop locations.
    Eg the sequence starts:  0  10 20
                    founder: 1  2  3
    passed with the locations [15, 25] should produce the sequence
                    starts:  0  10 15 20 25
                    founder: 1  2  2  3  3
    Assumes that locations is sorted. new_sequence works like a merge
    sort which only keeps unique values in output starts array, but
    with modifications to handle the founder array.
    """
    if locations.shape[0] == 0:
        return Diploid(diploid.starts.tolist(), diploid.end,
                       diploid.founder.tolist())
    cdef np.ndarray[np.uint32_t, ndim=1] starts = diploid.starts
    cdef np.ndarray[np.uint32_t, ndim=1] founder = diploid.founder
    cdef unsigned long break_index, break_location
    cdef list new_starts = []
    cdef list new_founder = []
    cdef unsigned long starts_i, locations_i
    cdef unsigned long total_len, locations_len, starts_len
    cdef unsigned long start_loci, break_loci
    if locations.shape[0] > 0 and locations[-1] == diploid.end:
        locations = locations[:-1]
    starts_i = 0
    locations_i = 0
    founder_i = 0
    total_i = 0
    locations_len = locations.shape[0]
    starts_len = starts.shape[0]
    total_len = starts_len + locations_len
    while (starts_i + locations_i) < total_len:
        if (locations_len == locations_i or
            (starts_len != starts_i and
             starts[starts_i] < locations[locations_i])):
            new_starts.append(starts[starts_i])
            new_founder.append(founder[starts_i])
            starts_i += 1
        elif (starts_len == starts_i or
              starts[starts_i] > locations[locations_i]):
            new_starts.append(locations[locations_i])
            new_founder.append(founder[starts_i - 1])
            locations_i += 1
        else: # starts[starts_i] == locations[locations_i])
            new_starts.append(starts[starts_i])
            new_founder.append(founder[starts_i])
            starts_i += 1
            locations_i += 1
    return Diploid(new_starts, diploid.end, new_founder)

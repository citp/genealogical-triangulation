"""# cython: profile=True"""
import numpy as np
from array import array

cimport numpy as np
cimport cython

cpdef common_segment_ibd(genome_a, genome_b):
    """
    Given two genomes returns a list of integers for each autosome,
    corresponding to the length of segments that are shared between
    the two autosomes.
    """
    cdef list ibd_segments = []
    ibd_segments.extend(common_homolog_segments(genome_a.mother,
                                            genome_b.mother))
    ibd_segments.extend(common_homolog_segments(genome_a.father,
                                           genome_b.mother))
    ibd_segments.extend(common_homolog_segments(genome_a.mother,
                                           genome_b.father))
    ibd_segments.extend(common_homolog_segments(genome_a.father,
                                           genome_b.father))
    return ibd_segments

def common_segment_lengths(genome_a, genome_b):
    """
    Given two genomes returns a list of integers for each autosome,
    corresponding to the length of segments that are shared between
    the two autosomes.
    """
    return lengths(common_segment_ibd(genome_a, genome_b))

# We don't need bounds checks, because the condition of the while loop
# ensures array access is within bounds.
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef list common_homolog_segments(homolog_a, homolog_b):
    """
    Given two autosome homologs, returns a list of ranges (a, b), (b, c), ...
    where the two autosomes have the same underlying sequence.
    """
    cdef unsigned long len_a, len_b, index_a, index_b, start, stop
    cdef unsigned long a_start, a_stop, a_id, b_start, b_stop, b_id
    cdef np.ndarray[np.uint32_t, ndim=1] starts_a = homolog_a.starts
    cdef np.ndarray[np.uint32_t, ndim=1] founder_a = homolog_a.founder
    cdef np.ndarray[np.uint32_t, ndim=1] starts_b = homolog_b.starts
    cdef np.ndarray[np.uint32_t, ndim=1] founder_b = homolog_b.founder
    cdef unsigned long end = homolog_a.end
    len_a = len(starts_a)
    len_b = len(starts_b)
    index_a = 0
    index_b = 0
    cdef list shared_segments = []
    while index_a < len_a and index_b < len_b:
        a_start = starts_a[index_a]
        if index_a + 1 < len_a:
            a_stop = starts_a[index_a + 1]
        else:
            a_stop = end
        a_id = founder_a[index_a]

        b_start = starts_b[index_b]
        if index_b + 1 < len_b:
            b_stop = starts_b[index_b + 1]
        else:
            b_stop = end
        # assert a_start < a_stop
        # assert b_start < b_stop
        b_id = founder_b[index_b]
        if a_id == b_id:
            if a_start > b_start:
                start = a_start
            else:
                start = b_start
            if a_stop < b_stop:
                stop = a_stop
            else:
                stop = b_stop
            shared_segments.append((start, stop))
        if a_stop == b_stop:
            index_a += 1
            index_b += 1
        elif a_stop > b_stop:
            index_b += 1
        else:
            index_a += 1
    if len(shared_segments) <= 1:
        return shared_segments
    # consolidate contiguous segments eg if we have shared segments
    # (0, 5) and (5, 10), then we should merge them into (0, 10).
    return _consolidate_sequence(shared_segments)



cpdef list lengths(list segments):
    """
    Takes a list of segments and returns a list of lengths.
    """
    cdef unsigned long a, b
    return [b - a for a, b in segments]

cpdef list _consolidate_sequence(sequence):
    """
    Takes a list of elements of the form (a, b), (c, d), ...  and
    merges elements where b = c such that (a, b), (c, d) becomes (a, d)
    """
    cdef int i, j
    assert len(sequence) > 1
    i = 0
    j = 1
    cdef list consolidated = []
    while j < len(sequence):
        if sequence[j - 1][1] != sequence[j][0]:
            consolidated.append((sequence[i][0], sequence[j - 1][1]))
            i = j
        j += 1
    consolidated.append((sequence[i][0], sequence[j - 1][1]))
    return consolidated

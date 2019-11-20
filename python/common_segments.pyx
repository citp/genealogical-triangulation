# cython: language_level=3
"""# cython: profile=True"""

import numpy as np
from array import array
from bisect import bisect_left
from collections import defaultdict

cimport numpy as np
cimport cython

cpdef common_segment_ibd(genome_a, genome_b, segment_filter):
    """
    Given two genomes returns a list of integers for each autosome,
    corresponding to the length of segments that are shared between
    the two autosomes.
    """
    # TODO: Move segment detector call inside this function
    cdef list ibd_segments = []
    cdef list b_father_ibd = []
    cdef list b_mother_ibd = []
    cdef list b_ibd = []
    cdef list a_mother_ibd
    cdef list a_father_ibd
    cdef list b_inbreed
    
    a_mother_ibd = segment_filter(common_homolog_segments(genome_a.mother,
                                                          genome_b.mother))
    a_father_ibd = segment_filter(common_homolog_segments(genome_a.father,
                                                          genome_b.mother))

    # Any repeat regions in the previously calculated IBD comes from
    # inbreeding in genome_a. To prevent overcounting inbreeding IBD,
    # we merge together the IBD segments.
    # The merging step is a part of handling the S1 case in:
    # https://web.archive.org/web/20190319202929/https://openi.nlm.nih.gov/imgs/512/114/3579841/PMC3579841_pone.0057003.g001.png
    b_mother_ibd = merge_overlaps(a_mother_ibd, a_father_ibd)
    
    a_mother_ibd = segment_filter(common_homolog_segments(genome_a.mother,
                                                          genome_b.father))
    a_father_ibd = segment_filter(common_homolog_segments(genome_a.father,
                                                          genome_b.father))

    b_father_ibd = merge_overlaps(a_mother_ibd, a_father_ibd)

    b_inbreed = common_homolog_segments(genome_b.mother, genome_b.father)
    if len(b_inbreed) > 0:
        a_inbreed = common_homolog_segments(genome_a.mother, genome_a.father)
        # Subtracting is so we are sure we aren't under counting in the S1 condition
        b_inbreed_only = subtract_regions(b_inbreed, a_inbreed)
        return remove_inbreeding(b_mother_ibd, b_father_ibd, b_inbreed_only)
    else:
        return b_mother_ibd + b_father_ibd

cpdef list remove_inbreeding(list ibd_a, list ibd_b, list inbreeding):
    cdef unsigned long a_overlap, b_overlap
    cdef tuple inbred_region
    # This function essentially handles the S3/S5 case of the
    # previously linked image
    for inbred_region in inbreeding:
        a_overlap = size_of_overlap(ibd_a, inbred_region)
        b_overlap = size_of_overlap(ibd_b, inbred_region)
        # Without segment filtering, if one overlap were nonzero, both
        # would be nonzero, but the filter may drop ibd segments in
        # inbred regions.
        if a_overlap == 0 or b_overlap == 0:
            continue
        if a_overlap < b_overlap:
            subtract_region(ibd_a, inbred_region)
        else:
            subtract_region(ibd_b, inbred_region)
    return ibd_a + ibd_b

cpdef unsigned long size_of_overlap(list a, tuple region):
    if len(a) == 0:
        return 0

    overlap = 0
    for a_region in a:
        start = max(a_region[0], region[0])
        end = min(a_region[1], region[1])
        size = end - start
        if 0 < size:
            overlap += size
    return overlap

cpdef list subtract_region(list a, tuple region):
    cdef long start_i
    cdef long stop_i
    cdef unsigned long a_len = len(a)
    if a_len == 0:
        return a
    start_i = 0
    stop_i = a_len - 1
    while start_i < a_len and a[start_i][1] <= region[0]:
        start_i += 1

    while 0 <= stop_i and region[1] <= a[stop_i][0]:
        stop_i -= 1

    # implicit in this check is a bounds check that the indicies are
    # between 0 and len(a) - 1
    if stop_i < start_i:
        return a

    to_insert = []
    if a[start_i][0] < region[0]:
        to_insert.append((a[start_i][0], region[0]))
    if region[1] < a[stop_i][1]:
        to_insert.append((region[1], a[stop_i][1]))
    del a[start_i:stop_i+1]
    while len(to_insert) > 0:
        a.insert(start_i, to_insert.pop())
    return a
    
cpdef list subtract_regions(list a, list b):
    """
    Subtract the regions in b from the regions in a
    """
    ret = sorted(a)
    for region in b:
        subtract_region(ret, region)
    return ret
    

cpdef list merge_overlaps(list ibd_a, list ibd_b):
    cdef list all_ibd = ibd_a + ibd_b
    cdef list ret = []
    if len(all_ibd) <= 1:
        return all_ibd
    all_ibd.sort(reverse = True)
    ret.append(all_ibd.pop())
    while len(all_ibd) > 0:
        top = ret[-1]
        current = all_ibd.pop()
        if current[0] <= top[1]:
            if current[1] < top[1]:
                end = top[1]
            else:
                end = current[1]
            ret.pop()
            ret.append((top[0], end))
        else:
            ret.append(current)
    return ret

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

cpdef update_segments(to_update, segments_dict):
    for founder, segments in segments_dict.items():
        to_update[founder].extend(segments)

cpdef common_segment_ibd_by_founders(genome_a, genome_b):
    """
    Given two genomes returns a list of integers for each autosome,
    corresponding to the length of segments that are shared between
    the two autosomes.
    """
    ibd_segments = defaultdict(list)
    d = common_homolog_segments_by_founder(genome_a.mother,
                                           genome_b.mother)
    update_segments(ibd_segments, d)

    d = common_homolog_segments_by_founder(genome_a.father,
                                           genome_b.mother)
    update_segments(ibd_segments, d)

    d = common_homolog_segments_by_founder(genome_a.mother,
                                           genome_b.father)
    update_segments(ibd_segments, d)

    d = common_homolog_segments_by_founder(genome_a.father,
                                           genome_b.father)
    update_segments(ibd_segments, d)

    return ibd_segments

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef common_homolog_segments_by_founder(homolog_a, homolog_b):
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
    shared_segments = defaultdict(list)
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
            shared_segments[a_id].append((start, stop))
        if a_stop == b_stop:
            index_a += 1
            index_b += 1
        elif a_stop > b_stop:
            index_b += 1
        else:
            index_a += 1
    # consolidate contiguous segments eg if we have shared segments
    # (0, 5) and (5, 10), then we should merge them into (0, 10).
    for key, shared in shared_segments.items():
        if len(shared) > 1:
            shared_segments[key] = _consolidate_sequence(shared)
    return shared_segments


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

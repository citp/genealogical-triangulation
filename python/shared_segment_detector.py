import numpy as np

from common_segments import common_segment_ibd
from cm import cm_lengths


class SharedSegmentDetector:
    """
    This class compares two genomes and returns the # of basepairs
    they share IBD. It handles detection cutoffs.
    """
    def __init__(self, minimum_base_length = 5000000, minimum_cm_length = 0,
                 recomb_data = None):
        assert minimum_base_length >= 0
        assert minimum_cm_length >= 0
        if minimum_cm_length > 0:
            assert recomb_data is not None
        self.minimum_base_length = minimum_base_length
        self.minimum_cm_length = minimum_cm_length
        self.recomb_data = recomb_data

    def shared_segment_length(self, genome_a, genome_b):
        segments = common_segment_ibd(genome_a, genome_b)
        if len(segments) == 0:
            return 0

        starts, stops = zip(*segments)
        np_starts = np.array(starts, dtype = np.uint32)
        np_stops = np.array(stops, dtype = np.uint32)
        lengths = np_stops - np_starts
        assert np.all(lengths >= 0)
        above_base_cutoff = lengths >= self.minimum_base_length

        if self.minimum_cm_length > 0:
            cm_l = cm_lengths(starts, stops, self.recomb_data)
            cm_cutoff = cm_l >= self.minimum_cm_length
        else:
            cm_cutoff = np.full(len(lengths), True)

        return np.sum(lengths[above_base_cutoff & cm_cutoff])

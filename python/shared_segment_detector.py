import numpy as np

from common_segments import common_segment_ibd
from cm import cm_lengths


class SharedSegmentDetector:
    """
    This class compares two genomes and returns the # of basepairs
    they share IBD. It handles detection cutoffs.
    """
    def __init__(self, recomb_data,
                 minimum_base_length = 5000000, minimum_cm_length = 0.0):
        assert minimum_base_length >= 0
        assert minimum_cm_length >= 0
        self.minimum_base_length = minimum_base_length
        self.minimum_cm_length = float(minimum_cm_length)
        self.recomb_data = recomb_data

    def _segment_filter(self, segments):
        if len(segments) == 0:
            return segments

        starts, stops = zip(*segments)
        np_starts = np.array(starts, dtype = np.uint32)
        np_stops = np.array(stops, dtype = np.uint32)
        lengths = np_stops - np_starts
        assert np.all(lengths >= 0)
        above_base_cutoff = lengths >= self.minimum_base_length

        if self.minimum_cm_length > 0.0:
            cm_l = cm_lengths(starts, stops, self.recomb_data)
            cm_cutoff = cm_l >= self.minimum_cm_length
        else:
            cm_cutoff = np.full(len(lengths), True)

        include_segments = above_base_cutoff & cm_cutoff
        return [segment for segment, boolean in zip(segments, include_segments)
                if boolean]

    def shared_segment_length(self, genome_a, genome_b):
        segments = common_segment_ibd(genome_a, genome_b, self._segment_filter)
        if len(segments) == 0:
            return 0.0
        starts, stops = list(zip(*segments))
        np_starts = np.array(starts, dtype = np.uint32)
        np_stops = np.array(stops, dtype = np.uint32)
        return float(np.sum(cm_lengths(starts, stops, self.recomb_data)))

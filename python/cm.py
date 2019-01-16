from collections import namedtuple

import numpy as np
import pdb

from recomb_genome import hapmap_filenames, CHROMOSOME_ORDER, _read_recombination_file

CentimorganData = namedtuple("CentimorganData", ["bases", "cm", "rates"])

def centimorgan_data_from_directory(directory):
    """
    Given a directory with hapmap data, returns CentimorganData for data.
    TODO: Some of this functionality is redundant with functionality in
    recomb_genome. A refactor at some point should clean this up.
    """
    file_names = hapmap_filenames(directory)
    bp = []
    cm = []
    rates = []
    bp_accum = 0
    cm_accum = 0
    for chrom in CHROMOSOME_ORDER:
        filename = file_names[chrom]
        data = _read_recombination_file(filename)
        for loci, rate, cumulative_cm in data:
            bp.append(loci + bp_accum)
            cm.append(cumulative_cm + cm_accum)
            rates.append(rate / 1000000)
        bp_accum += data[-1][0]
        cm_accum += data[-1][2]

    np_bases = np.array(bp, dtype = np.uint32)
    np_cm = np.array(cm)
    np_rates = np.array(rates)
    return CentimorganData(np_bases, np_cm, np_rates)
        
        
def cumulative_cm(locations, recombination_data):
    """
    Calculates the cumulative centimorgans for the given genome
    locations based on the given recombination data.
    """
    np_locations = np.array(locations, dtype = np.uint32, copy = False)
    base_ends = recombination_data.bases
    assert np.all(np_locations <= base_ends[-1])
    end_index = np.searchsorted(base_ends, np_locations, side = "left")
    
    cm_ends = recombination_data.cm
    cm_distance = cm_ends[end_index]
    
    bp_difference = base_ends[end_index] - locations
    rates = recombination_data.rates
    cm_difference = bp_difference * rates[end_index - 1]
    cm_difference[end_index == 0] = 0

    adjusted_cm_distance = cm_distance - cm_difference
    assert np.all(0 <= adjusted_cm_distance) , "Expected positive cM distances, got {}".format(adjusted_cm_distance[adjusted_cm_distance < 0])
    assert np.all(adjusted_cm_distance <= cm_ends[-1])
    return adjusted_cm_distance

def cm_lengths(starts, stops, recombination_data):
    """
    Computes the centimorgan length for regions of the genome defined
    by starts and stops, where stop[i] corresponds to start[i]
    """
    np_starts = np.array(starts, dtype = np.uint32, copy = False)
    np_stops = np.array(stops, dtype = np.uint32, copy = False)
    cm_starts = cumulative_cm(np_starts, recombination_data)
    cm_stops = cumulative_cm(np_stops, recombination_data)
    return cm_stops - cm_starts

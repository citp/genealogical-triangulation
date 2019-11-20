#!/usr/bin/env python3

import unittest
from unittest.mock import MagicMock

import numpy as np

# import pyximport; pyximport.install()
from common_segments import common_homolog_segments, _consolidate_sequence, merge_overlaps, subtract_region, size_of_overlap, remove_inbreeding

from cm import cm_lengths, cumulative_cm

uint32 = np.uint32

class TestCommonHomologSegments(unittest.TestCase):
    def test_single_same_segment(self):
        a = MagicMock()
        a.starts = np.array([0], dtype = uint32)
        a.founder = np.array([0], dtype = uint32)
        a.end = 10
        shared = common_homolog_segments(a, a)
        self.assertEqual(shared, [(0, 10)])

    def test_single_different_segment(self):
        a = MagicMock()
        a.starts = np.array([0], dtype = uint32)
        a.founder = np.array([0], dtype = uint32)
        a.end = 10
        b = MagicMock()
        b.starts = np.array([0], dtype = uint32)
        b.founder = np.array([1], dtype = uint32)
        b.end = 10
        shared = common_homolog_segments(a, b)
        self.assertEqual(shared, [])

    def test_multiple_same_segment(self):
        a = MagicMock()
        a.starts = np.array([0], dtype = uint32)
        a.founder = np.array([0], dtype = uint32)
        a.end = 10
        b = MagicMock()
        b.starts = np.array([0, 5], dtype = uint32)
        b.founder = np.array([0, 0], dtype = uint32)
        b.end = 10
        shared = common_homolog_segments(a, b)
        self.assertEqual(shared, [(0, 10)])

    def test_two_same_segments(self):
        a = MagicMock()
        a.starts = np.array([0, 5], dtype = uint32)
        a.founder = np.array([1, 2], dtype = uint32)
        a.end = 10
        shared = common_homolog_segments(a, a)
        self.assertEqual(shared, [(0, 10)])

    def test_two_segments_different_boundary(self):
        a = MagicMock()
        a.starts = np.array([0, 5], dtype = uint32)
        a.founder = np.array([1, 1], dtype = uint32)
        a.end = 10
        b = MagicMock()
        b.starts = np.array([0, 6], dtype = uint32)
        b.founder = np.array([1, 1], dtype = uint32)
        b.end = 10
        shared = common_homolog_segments(a, b)
        self.assertEqual(shared, [(0, 10)])

    def test_single_element_vs_many_match_in_back(self):
        a = MagicMock()
        a.starts = np.array([0], dtype = uint32)
        a.founder = np.array([0], dtype = uint32)
        a.end = 10
        b = MagicMock()
        b.starts = np.array([0, 2, 4, 8], dtype = uint32)
        b.founder = np.array([1, 2, 3, 0], dtype = uint32)
        b.end = 10
        shared = common_homolog_segments(a, b)
        self.assertEqual(shared, [(8, 10)])

    def test_single_element_vs_many_match_in_front(self):
        a = MagicMock()
        a.starts = np.array([0], dtype = uint32)
        a.founder = np.array([0], dtype = uint32)
        a.end = 10
        b = MagicMock()
        b.starts = np.array([0, 2, 4, 8], dtype = uint32)
        b.founder = np.array([0, 2, 3, 4], dtype = uint32)
        b.end = 10
        shared = common_homolog_segments(a, b)
        self.assertEqual(shared, [(0, 2)])

class TestConsolidateSequence(unittest.TestCase):
    def test_two_elements_merge(self):
        seq = [(0, 5), (5, 10)]
        con = _consolidate_sequence(seq)
        self.assertEqual(con, [(0, 10)])

    def test_two_elements_disjoint(self):
        seq = [(0, 5), (6, 10)]
        con = _consolidate_sequence(seq)
        self.assertEqual(con, seq)

    def test_first_two_merge(self):
        seq = [(0, 4), (4, 8), (9, 10)]
        con = _consolidate_sequence(seq)
        self.assertEqual(con, [(0, 8), (9, 10)])

    def test_last_two_merge(self):
        seq = [(0, 4), (5, 8), (8, 10)]
        con = _consolidate_sequence(seq)
        self.assertEqual(con, [(0, 4), (5, 10)])
        
    def test_middle_two_merge(self):
        seq = [(0, 3), (4, 6), (6, 8), (9, 10)]
        con = _consolidate_sequence(seq)
        self.assertEqual(con, [(0, 3), (4, 8), (9, 10)])

    def test_many_elements(self):
        seq = [(0, 2), (2, 4), (4, 8), (8, 10)]
        con = _consolidate_sequence(seq)
        self.assertEqual(con, [(0, 10)])

    def test_many_emements_single_point(self):
        seq = list(zip(range(0, 10), range(1, 11)))
        con = _consolidate_sequence(seq)
        self.assertEqual(con, [(0, 10)])

class TestCumulativeCm(unittest.TestCase):
    def test_single_length(self):
        """
        data
        5	0.0002000000 0.00000
        10	0.0000500000 0.00100
        15	0.0001800000 0.00125
        20	0.0002000000 0.00215
        """

        recombination = MagicMock()
        recombination.bases = np.array([5, 10, 15, 20], dtype = np.uint32)
        recombination.cm = np.array([0.00000, 0.001, 0.00125, 0.00215])
        recombination.rates = np.array([0.0002000000, 0.0000500000,
                                        0.0001800000, 0.0002000000])

        locations = [5]
        lengths = cumulative_cm(locations, recombination)
        expect = np.array([0.0], dtype = np.float64)
        np.testing.assert_almost_equal(lengths, expect)

        locations = [10]
        lengths = cumulative_cm(locations, recombination)
        expect = np.array([0.001], dtype = np.float64)
        np.testing.assert_almost_equal(lengths, expect)

        locations = [15]
        lengths = cumulative_cm(locations, recombination)
        expect = np.array([0.00125], dtype = np.float64)
        np.testing.assert_almost_equal(lengths, expect)

        locations = [20]
        lengths = cumulative_cm(locations, recombination)
        expect = np.array([0.00215], dtype = np.float64)
        np.testing.assert_almost_equal(lengths, expect)

        locations = [3]
        lengths = cumulative_cm(locations, recombination)
        expect = np.array([0.00], dtype = np.float64)
        np.testing.assert_almost_equal(lengths, expect)

        locations = [4]
        lengths = cumulative_cm(locations, recombination)
        expect = np.array([0.0], dtype = np.float64)
        np.testing.assert_almost_equal(lengths, expect)

        locations = [9]
        lengths = cumulative_cm(locations, recombination)
        expect = np.array([0.0008], dtype = np.float64)
        np.testing.assert_almost_equal(lengths, expect)

        locations = [11]
        lengths = cumulative_cm(locations, recombination)
        expect = np.array([0.00105], dtype = np.float64)
        np.testing.assert_almost_equal(lengths, expect)

        locations = [14]
        lengths = cumulative_cm(locations, recombination)
        expect = np.array([0.0012], dtype = np.float64)
        np.testing.assert_almost_equal(lengths, expect)

        locations = [19]
        lengths = cumulative_cm(locations, recombination)
        expect = np.array([0.00197], dtype = np.float64)
        np.testing.assert_almost_equal(lengths, expect)

    def test_multiple_length(self):
        """
        data
        5	0.0002000000 0.00000
        10	0.0000500000 0.00100
        15	0.0001800000 0.00125
        20	0.0002000000 0.00215
        """

        recombination = MagicMock()
        recombination.bases = np.array([5, 10, 15, 20], dtype = np.uint32)
        recombination.cm = np.array([0.00000, 0.001, 0.00125, 0.00215])
        recombination.rates = np.array([0.0002000000, 0.0000500000,
                                        0.0001800000, 0.0002000000])

        locations = [5, 10, 15, 20, 3, 4, 9, 11, 14, 19]
        expect = np.array([0.0, 0.001, 0.00125, 0.00215, 0.00, 0.0,
                           0.0008, 0.00105, 0.0012, 0.00197])
        lengths = cumulative_cm(locations, recombination)
        np.testing.assert_almost_equal(lengths, expect)
        

class TestCmLengths(unittest.TestCase):
    def test_single_segment(self):
        """
        data
        5	0.0002000000 0.00000
        10	0.0000500000 0.00100
        15	0.0001800000 0.00125
        20	0.0002000000 0.00215
        """

        recombination = MagicMock()
        recombination.bases = np.array([5, 10, 15, 20], dtype = np.uint32)
        recombination.cm = np.array([0.00000, 0.001, 0.00125, 0.00215])
        recombination.rates = np.array([0.0002000000, 0.0000500000,
                                        0.0001800000, 0.0002000000])

        starts = [0]
        stops = [5]
        lengths = cm_lengths(starts, stops, recombination)
        expect = np.array([0.000], dtype = np.float64)
        np.testing.assert_almost_equal(lengths, expect)

        starts = [5]
        stops = [10]
        lengths = cm_lengths(starts, stops, recombination)
        expect = np.array([0.001], dtype = np.float64)
        np.testing.assert_almost_equal(lengths, expect)

        starts = [10]
        stops = [15]
        lengths = cm_lengths(starts, stops, recombination)
        expect = np.array([0.00025], dtype = np.float64)
        np.testing.assert_almost_equal(lengths, expect)

        starts = [15]
        stops = [20]
        lengths = cm_lengths(starts, stops, recombination)
        expect = np.array([0.0009], dtype = np.float64)
        np.testing.assert_almost_equal(lengths, expect)

        starts = [11]
        stops = [14]
        lengths = cm_lengths(starts, stops, recombination)
        expect = np.array([0.00015], dtype = np.float64)
        np.testing.assert_almost_equal(lengths, expect)

        starts = [11]
        stops = [19]
        lengths = cm_lengths(starts, stops, recombination)
        expect = np.array([0.0009199999999999999], dtype = np.float64)
        np.testing.assert_almost_equal(lengths, expect)

        starts = [11]
        stops = [18]
        lengths = cm_lengths(starts, stops, recombination)
        expect = np.array([0.0007399999999999999], dtype = np.float64)
        np.testing.assert_almost_equal(lengths, expect)

        starts = [0]
        stops = [20]
        lengths = cm_lengths(starts, stops, recombination)
        expect = np.array([0.00215], dtype = np.float64)
        np.testing.assert_almost_equal(lengths, expect)
        
class TestMergeOverlaps(unittest.TestCase):
    def test_empty(self):
        self.assertEqual(merge_overlaps([], []), [])
        self.assertEqual(merge_overlaps([], [(1, 2)]), [(1, 2)])
        self.assertEqual(merge_overlaps([(1, 2)], []), [(1, 2)])

    def test_single(self):
        self.assertEqual(merge_overlaps([(1, 2)], [(1, 2)]), [(1, 2)])
        self.assertEqual(merge_overlaps([(1, 2)], [(4, 5)]), [(1, 2), (4, 5)])
        self.assertEqual(merge_overlaps([(4, 5)], [(1, 2)]), [(1, 2), (4, 5)])
        self.assertEqual(merge_overlaps([(1, 2)], [(1, 4)]), [(1, 4)])
        self.assertEqual(merge_overlaps([(1, 4)], [(1, 2)]), [(1, 4)])

    def test_many(self):
        self.assertEqual(merge_overlaps([(1, 2), (5, 8), (10, 15), (20, 25)],
                                        [(9, 18)]),
                         [(1, 2), (5, 8), (9, 18), (20, 25)])

class TestSubtractRegion(unittest.TestCase):

    def test_single(self):
        self.assertEqual(subtract_region([(1, 2)], (1, 2)), [])
        self.assertEqual(subtract_region([(1, 10)], (5, 10)), [(1, 5)])
        self.assertEqual(subtract_region([(1, 10)], (1, 5)), [(5, 10)])
        self.assertEqual(subtract_region([(1, 10)], (2, 8)), [(1, 2), (8, 10)])
        self.assertEqual(subtract_region([(1, 10)], (15, 20)), [(1, 10)])
        self.assertEqual(subtract_region([(10, 20)], (1, 10)), [(10, 20)])

    def test_multiple(self):
        # Exact matches
        self.assertEqual(subtract_region([(1, 10), (20, 30), (40, 50)],
                                         (1, 10)), [(20, 30), (40, 50)])
        self.assertEqual(subtract_region([(1, 10), (20, 30), (40, 50)],
                                         (20, 30)), [(1, 10), (40, 50)])
        self.assertEqual(subtract_region([(1, 10), (20, 30), (40, 50)],
                                         (40, 50)), [(1, 10), (20, 30)])

        # no matches
        self.assertEqual(subtract_region([(1, 10), (20, 30), (40, 50)],
                                         (60, 70)),
                         [(1, 10), (20, 30), (40, 50)])
        self.assertEqual(subtract_region([(1, 10), (20, 30), (40, 50)],
                                         (0, 1)),
                         [(1, 10), (20, 30), (40, 50)])
        self.assertEqual(subtract_region([(1, 10), (20, 30), (40, 50)],
                                         (12, 18)),
                         [(1, 10), (20, 30), (40, 50)])
        

        # sub region matches
        self.assertEqual(subtract_region([(1, 10), (20, 30), (40, 50)],
                                         (5, 10)),
                         [(1, 5), (20, 30), (40, 50)])

        self.assertEqual(subtract_region([(1, 10), (20, 30), (40, 50)],
                                         (5, 12)),
                         [(1, 5), (20, 30), (40, 50)])

    def test_remove_multiple(self):
        self.assertEqual(subtract_region([(1, 10), (20, 30), (40, 50)],
                                         (5, 25)),
                         [(1, 5), (25, 30), (40, 50)])
        self.assertEqual(subtract_region([(1, 10), (20, 30), (40, 50)],
                                         (5, 45)),
                         [(1, 5), (45, 50)])
        self.assertEqual(subtract_region([(1, 10), (20, 30), (40, 50)],
                                         (5, 50)),
                         [(1, 5)])


class TestSizeOfOverlap(unittest.TestCase):

    def test_single(self):
        self.assertEqual(size_of_overlap([(1, 2)], (1, 2)), 1)
        self.assertEqual(size_of_overlap([(1, 10)], (5, 10)), 5)
        self.assertEqual(size_of_overlap([(1, 10)], (1, 5)), 4)
        self.assertEqual(size_of_overlap([(1, 10)], (2, 8)), 6)
        self.assertEqual(size_of_overlap([(1, 10)], (15, 20)), 0)
        self.assertEqual(size_of_overlap([(10, 20)], (1, 10)), 0)

    def test_multiple(self):
        # Exact matches
        self.assertEqual(size_of_overlap([(1, 10), (20, 30), (40, 50)],
                                         (1, 10)), 9)
        self.assertEqual(size_of_overlap([(1, 10), (20, 30), (40, 50)],
                                         (20, 30)), 10)
        self.assertEqual(size_of_overlap([(1, 10), (20, 30), (40, 50)],
                                         (40, 50)), 10)

        # no matches
        self.assertEqual(size_of_overlap([(1, 10), (20, 30), (40, 50)],
                                         (60, 70)), 0)
        self.assertEqual(size_of_overlap([(1, 10), (20, 30), (40, 50)],
                                         (0, 1)), 0)
        self.assertEqual(size_of_overlap([(1, 10), (20, 30), (40, 50)],
                                         (12, 18)), 0)
        

        # sub region matches
        self.assertEqual(size_of_overlap([(1, 10), (20, 30), (40, 50)],
                                         (5, 10)), 5)

        self.assertEqual(size_of_overlap([(1, 10), (20, 30), (40, 50)],
                                         (5, 12)), 5)

    def test_remove_multiple(self):
        self.assertEqual(size_of_overlap([(1, 10), (20, 30), (40, 50)],
                                         (5, 25)), 10)
        self.assertEqual(size_of_overlap([(1, 10), (20, 30), (40, 50)],
                                         (5, 45)), 20)
        self.assertEqual(size_of_overlap([(1, 10), (20, 30), (40, 50)],
                                         (5, 50)), 25)

class TestRemoveInbreeding(unittest.TestCase):
    def test_single(self):
        self.assertEqual(remove_inbreeding([(10, 20)], [(10, 20)], [(10, 20)]),
                         [(10, 20)])
        self.assertEqual(remove_inbreeding([(10, 20)], [(10, 20)], []),
                         [(10, 20), (10, 20)])
        self.assertEqual(remove_inbreeding([(10, 20)], [(10, 20)], [(30, 40)]),
                         [(10, 20), (10, 20)])
        self.assertEqual(remove_inbreeding([(10, 20), (20, 30)],
                                           [(10, 20)], [(10, 20)]),
                         [(10, 20), (20, 30)])
        self.assertEqual(remove_inbreeding([(10, 20)],
                                           [(10, 20), (20, 30)], [(10, 20)]),
                         [(10, 20), (20, 30)])

    
if __name__ == '__main__':
    unittest.main()

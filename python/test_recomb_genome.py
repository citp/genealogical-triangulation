#!/usr/bin/env python3

from bisect import bisect_left
from unittest.mock import MagicMock
import unittest

import numpy as np
# import pyximport; pyximport.install()

import recomb_genome
from recomb_helper import new_sequence

def ar(locs):
    return np.array(locs, dtype = np.uint32)

def break_sequence_wrapper(sequence, location):
    return recomb_genome._break_sequence(sequence, location,
                                         bisect_left(sequence, (location,)))


class TestNewSequence(unittest.TestCase):
    def test_single_element_middle(self):
        diploid = MagicMock()
        diploid.starts = np.array([0], dtype = np.uint32)
        diploid.end = 10
        diploid.founder = np.array([1], dtype = np.uint32)
        ret_diploid = new_sequence(diploid, ar([5]))
        self.assertEqual(ret_diploid.starts, [0, 5])
        self.assertEqual(ret_diploid.founder, [1, 1])


    def test_single_element_start(self):
        diploid = MagicMock()
        diploid.starts = np.array([0], dtype = np.uint32)
        diploid.end = 10
        diploid.founder = np.array([1], dtype = np.uint32)
        ret_diploid = new_sequence(diploid, ar([0]))
        self.assertEqual(ret_diploid.starts, [0])
        self.assertEqual(ret_diploid.founder, [1])

    def test_single_element_multiple(self):
        diploid = MagicMock()
        diploid.starts = np.array([0], dtype = np.uint32)
        diploid.end = 10
        diploid.founder = np.array([1], dtype = np.uint32)
        ret_diploid = new_sequence(diploid, ar([0, 4, 6, 8]))
        self.assertEqual(ret_diploid.starts, [0, 4, 6, 8])
        self.assertEqual(ret_diploid.founder, [1, 1, 1, 1])

    def test_end_boundary_two_element(self):
        diploid = MagicMock()
        diploid.starts = np.array([0, 10], dtype = np.uint32)
        diploid.end = 20
        diploid.founder = np.array([1, 2], dtype = np.uint32)
        ret_diploid = new_sequence(diploid, ar([10]))
        self.assertEqual(ret_diploid.starts, [0, 10])
        self.assertEqual(ret_diploid.founder, [1, 2])

    def test_start_boundary_two_element(self):
        diploid = MagicMock()
        diploid.starts = np.array([0, 10], dtype = np.uint32)
        diploid.end = 20
        diploid.founder = np.array([1, 2], dtype = np.uint32)
        ret_diploid = new_sequence(diploid, ar([0]))
        self.assertEqual(ret_diploid.starts, [0, 10])
        self.assertEqual(ret_diploid.founder, [1, 2])

    def test_middle_boundary_two_element(self):
        diploid = MagicMock()
        diploid.starts = np.array([0, 10], dtype = np.uint32)
        diploid.end = 20
        diploid.founder = np.array([1, 2], dtype = np.uint32)
        ret_diploid = new_sequence(diploid, ar([5]))
        self.assertEqual(ret_diploid.starts, [0, 5, 10])
        self.assertEqual(ret_diploid.founder, [1, 1, 2])

    def test_middle_boundary_two_element_multiple_breaks(self):
        diploid = MagicMock()
        diploid.starts = np.array([0, 10], dtype = np.uint32)
        diploid.end = 20
        diploid.founder = np.array([1, 2], dtype = np.uint32)
        ret_diploid = new_sequence(diploid, ar([5, 15]))
        self.assertEqual(ret_diploid.starts, [0, 5, 10, 15])
        self.assertEqual(ret_diploid.founder, [1, 1, 2, 2])

    def test_single_element_end(self):
        diploid = MagicMock()
        diploid.starts = np.array([0], dtype = np.uint32)
        diploid.end = 10
        diploid.founder = np.array([1], dtype = np.uint32)
        ret_diploid = new_sequence(diploid, ar([10]))
        self.assertEqual(ret_diploid.starts, [0])
        self.assertEqual(ret_diploid.founder, [1])

    def test_empty_locations(self):
        diploid = MagicMock()
        diploid.starts = np.array([0], dtype = np.uint32)
        diploid.end = 10
        diploid.founder = np.array([1], dtype = np.uint32)
        ret_diploid = new_sequence(diploid, ar([]))
        self.assertEqual(ret_diploid.starts, [0])
        self.assertEqual(ret_diploid.founder, [1])

class TestSwapAtLocations(unittest.TestCase):
    def test_single_location_left_boundary(self):
        mother = MagicMock()
        mother.starts = np.array([0], dtype = np.uint32)
        mother.end = 10
        mother.founder = np.array([1], dtype = np.uint32)
        father = MagicMock()
        father.starts = np.array([0], dtype = np.uint32)
        father.end = 10
        father.founder = np.array([2], dtype = np.uint32)
        locations = [(0, 5)]
        new_mother, new_father = recomb_genome._swap_at_locations(mother,
                                                                  father,
                                                                  locations)
        np.testing.assert_array_equal(new_mother.starts,
                                      np.array([0, 5], dtype = np.uint32))
        np.testing.assert_array_equal(new_mother.founder,
                                      np.array([2, 1], dtype = np.uint32))
        np.testing.assert_array_equal(new_father.starts,
                                      np.array([0, 5], dtype = np.uint32))
        np.testing.assert_array_equal(new_father.founder,
                                      np.array([1, 2], dtype = np.uint32))
        
    def test_single_location_right_boundary(self):
        mother = MagicMock()
        mother.starts = np.array([0], dtype = np.uint32)
        mother.end = 10
        mother.founder = np.array([1], dtype = np.uint32)
        father = MagicMock()
        father.starts = np.array([0], dtype = np.uint32)
        father.end = 10
        father.founder = np.array([2], dtype = np.uint32)
        locations = [(5, 10)]
        new_mother, new_father = recomb_genome._swap_at_locations(mother,
                                                                  father,
                                                                  locations)
        np.testing.assert_array_equal(new_mother.starts,
                                      np.array([0, 5], dtype = np.uint32))
        np.testing.assert_array_equal(new_mother.founder,
                         np.array([1, 2], dtype = np.uint32))
        np.testing.assert_array_equal(new_father.starts,
                                      np.array([0, 5], dtype = np.uint32))
        np.testing.assert_array_equal(new_father.founder,
                                      np.array([2, 1], dtype = np.uint32))

    def test_single_location_middle_boundary(self):
        mother = MagicMock()
        mother.starts = np.array([0], dtype = np.uint32)
        mother.end = 10
        mother.founder = np.array([1], dtype = np.uint32)
        father = MagicMock()
        father.starts = np.array([0], dtype = np.uint32)
        father.end = 10
        father.founder = np.array([2], dtype = np.uint32)
        locations = [(2, 8)]
        new_mother, new_father = recomb_genome._swap_at_locations(mother,
                                                                  father,
                                                                  locations)
        np.testing.assert_array_equal(new_mother.starts,
                         np.array([0, 2, 8], dtype = np.uint32))
        np.testing.assert_array_equal(new_mother.founder,
                         np.array([1, 2, 1], dtype = np.uint32))
        np.testing.assert_array_equal(new_father.starts,
                         np.array([0, 2, 8], dtype = np.uint32))
        np.testing.assert_array_equal(new_father.founder,
                         np.array([2, 1, 2], dtype = np.uint32))

    def test_multiple_locations_single_segment(self):
        mother = MagicMock()
        mother.starts = np.array([0], dtype = np.uint32)
        mother.end = 10
        mother.founder = np.array([1], dtype = np.uint32)
        father = MagicMock()
        father.starts = np.array([0], dtype = np.uint32)
        father.end = 10
        father.founder = np.array([2], dtype = np.uint32)
        locations = [(0, 4), (6, 10)]
        new_mother, new_father = recomb_genome._swap_at_locations(mother,
                                                                  father,
                                                                  locations)
        np.testing.assert_array_equal(new_mother.starts,
                                      np.array([0, 4, 6], dtype = np.uint32))
        np.testing.assert_array_equal(new_mother.founder,
                                      np.array([2, 1, 2], dtype = np.uint32))
        np.testing.assert_array_equal(new_father.starts,
                                      np.array([0, 4, 6], dtype = np.uint32))
        np.testing.assert_array_equal(new_father.founder,
                                      np.array([1, 2, 1], dtype = np.uint32))

    def test_single_location_two_segments_first_segment(self):
        mother = MagicMock()
        mother.starts = np.array([0, 10], dtype = np.uint32)
        mother.end = 20
        mother.founder = np.array([1, 2], dtype = np.uint32)
        father = MagicMock()
        father.starts = np.array([0, 10], dtype = np.uint32)
        father.end = 20
        father.founder = np.array([3, 4], dtype = np.uint32)
        locations = [(0, 10)]
        new_mother, new_father = recomb_genome._swap_at_locations(mother,
                                                                  father,
                                                                  locations)
        np.testing.assert_array_equal(new_mother.starts,
                         np.array([0, 10], dtype = np.uint32))
        np.testing.assert_array_equal(new_mother.founder,
                         np.array([3, 2], dtype = np.uint32))
        np.testing.assert_array_equal(new_father.starts,
                         np.array([0, 10], dtype = np.uint32))
        np.testing.assert_array_equal(new_father.founder,
                         np.array([1, 4], dtype = np.uint32))

    def test_single_location_two_segments_last_segment(self):
        mother = MagicMock()
        mother.starts = np.array([0, 10], dtype = np.uint32)
        mother.end = 20
        mother.founder = np.array([1, 2], dtype = np.uint32)
        father = MagicMock()
        father.starts = np.array([0, 10], dtype = np.uint32)
        father.end = 20
        father.founder = np.array([3, 4], dtype = np.uint32)
        locations = [(10, 20)]
        new_mother, new_father = recomb_genome._swap_at_locations(mother,
                                                                  father,
                                                                  locations)
        np.testing.assert_array_equal(new_mother.starts,
                         np.array([0, 10], dtype = np.uint32))
        np.testing.assert_array_equal(new_mother.founder,
                         np.array([1, 4], dtype = np.uint32))
        np.testing.assert_array_equal(new_father.starts,
                         np.array([0, 10], dtype = np.uint32))
        np.testing.assert_array_equal(new_father.founder,
                         np.array([3, 2], dtype = np.uint32))


    def test_single_location_overlapping_two_segments(self):
        mother = MagicMock()
        mother.starts = np.array([0, 10], dtype = np.uint32)
        mother.end = 20
        mother.founder = np.array([1, 2], dtype = np.uint32)
        father = MagicMock()
        father.starts = np.array([0, 10], dtype = np.uint32)
        father.end = 20
        father.founder = np.array([3, 4], dtype = np.uint32)
        locations = [(5, 15)]
        new_mother, new_father = recomb_genome._swap_at_locations(mother,
                                                                  father,
                                                                  locations)
        np.testing.assert_array_equal(new_mother.starts,
                                      np.array([0, 5, 10, 15],
                                               dtype = np.uint32))
        np.testing.assert_array_equal(new_mother.founder,
                                      np.array([1, 3, 4, 2],
                                               dtype = np.uint32))
        np.testing.assert_array_equal(new_father.starts,
                                      np.array([0, 5, 10, 15],
                                               dtype = np.uint32))
        np.testing.assert_array_equal(new_father.founder,
                                      np.array([3, 1, 2, 4],
                                               dtype = np.uint32))

    def test_two_locations_at_boundary_two_segments (self):
        mother = MagicMock()
        mother.starts = np.array([0, 10], dtype = np.uint32)
        mother.end = 20
        mother.founder = np.array([1, 2], dtype = np.uint32)
        father = MagicMock()
        father.starts = np.array([0, 10], dtype = np.uint32)
        father.end = 20
        father.founder = np.array([3, 4], dtype = np.uint32)
        locations = [(0, 5), (15, 20)]
        new_mother, new_father = recomb_genome._swap_at_locations(mother,
                                                                  father,
                                                                  locations)
        np.testing.assert_array_equal(new_mother.starts,
                                      np.array([0, 5, 10, 15],
                                               dtype = np.uint32))
        np.testing.assert_array_equal(new_mother.founder,
                                      np.array([3, 1, 2, 4],
                                               dtype = np.uint32))
        np.testing.assert_array_equal(new_father.starts,
                                      np.array([0, 5, 10, 15],
                                               dtype = np.uint32))
        np.testing.assert_array_equal(new_father.founder,
                                      np.array([1, 3, 4, 2],
                                               dtype = np.uint32))

    def test_two_locations_middle_two_segments(self):
        mother = MagicMock()
        mother.starts = np.array([0, 10], dtype = np.uint32)
        mother.end = 20
        mother.founder = np.array([1, 2], dtype = np.uint32)
        father = MagicMock()
        father.starts = np.array([0, 10], dtype = np.uint32)
        father.end = 20
        father.founder = np.array([3, 4], dtype = np.uint32)
        locations = [(2, 5), (15, 18)]
        new_mother, new_father = recomb_genome._swap_at_locations(mother,
                                                                  father,
                                                                  locations)
        np.testing.assert_array_equal(new_mother.starts,
                                      np.array([0, 2, 5, 10, 15, 18],
                                               dtype = np.uint32))
        np.testing.assert_array_equal(new_mother.founder,
                                      np.array([1, 3, 1, 2, 4, 2],
                                               dtype = np.uint32))
        np.testing.assert_array_equal(new_father.starts,
                                      np.array([0, 2, 5, 10, 15, 18],
                                               dtype = np.uint32))
        np.testing.assert_array_equal(new_father.founder,
                                      np.array([3, 1, 3, 4, 2, 4],
                                               dtype = np.uint32))

    
    def test_single_location_two_segments_uneven(self):
        mother = MagicMock()
        mother.starts = np.array([0, 5], dtype = np.uint32)
        mother.end = 20
        mother.founder = np.array([1, 2], dtype = np.uint32)
        father = MagicMock()
        father.starts = np.array([0, 10], dtype = np.uint32)
        father.end = 20
        father.founder = np.array([3, 4], dtype = np.uint32)
        locations = [(2, 11)]
        new_mother, new_father = recomb_genome._swap_at_locations(mother,
                                                                  father,
                                                                  locations)
        np.testing.assert_array_equal(new_mother.starts,
                                      np.array([0, 2, 10, 11],
                                               dtype = np.uint32))
        np.testing.assert_array_equal(new_mother.founder,
                                      np.array([1, 3, 4, 2],
                                               dtype = np.uint32))
        np.testing.assert_array_equal(new_father.starts,
                                      np.array([0, 2, 5, 11],
                                               dtype = np.uint32))
        np.testing.assert_array_equal(new_father.founder,
                                      np.array([3, 1, 2, 4],
                                               dtype = np.uint32))



    
if __name__ == '__main__':
    unittest.main()

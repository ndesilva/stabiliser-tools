import numpy as np
import unittest
import benchmarking.brute_force.stab_state_check as ssc
import math

class Test_Brute_Force_Stabiliser_Check(unittest.TestCase):
    
    def test_is_stabiliser_state_accepts(self):
        state_vector = np.array([1] * 8)
        self.assertTrue(ssc.is_stab_from_paulis(state_vector))

    def test_is_stabiliser_state_rejects(self):
        state_vector = np.array(
        [1, 1, 1, 1, 1, 1, 1, (1/math.sqrt(2))*(1+1j)])/math.sqrt(8)
        self.assertFalse(ssc.is_stab_from_paulis(state_vector))
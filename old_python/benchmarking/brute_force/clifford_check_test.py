import unittest
import numpy as np
import benchmarking.brute_force.clifford_check as cc

class Test_Brute_Force_Clifford_Check(unittest.TestCase):

    @staticmethod
    def get_three_qubit_clifford() -> np.ndarray:
        return np.array([[ 0. +0.5j,  0. +0.5j,  0. +0.j ,  0. +0.j ,  0. -0.5j,  0. -0.5j,
         0. +0.j ,  0. +0.j ],
       [ 0. +0.j ,  0. +0.j ,  0. +0.5j,  0. -0.5j,  0. +0.j ,  0. +0.j ,
         0. -0.5j,  0. +0.5j],
       [ 0. +0.j ,  0. +0.j ,  0. +0.5j,  0. +0.5j,  0. +0.j ,  0. +0.j ,
         0. -0.5j,  0. -0.5j],
       [ 0. +0.5j,  0. -0.5j,  0. +0.j ,  0. +0.j ,  0. -0.5j,  0. +0.5j,
         0. +0.j ,  0. +0.j ],
       [ 0. +0.j ,  0. +0.j ,  0.5+0.j ,  0.5+0.j ,  0. +0.j ,  0. +0.j ,
         0.5+0.j ,  0.5+0.j ],
       [ 0.5+0.j , -0.5+0.j ,  0. +0.j ,  0. +0.j ,  0.5+0.j , -0.5+0.j ,
         0. +0.j ,  0. +0.j ],
       [ 0.5+0.j ,  0.5+0.j ,  0. +0.j ,  0. +0.j ,  0.5+0.j ,  0.5+0.j ,
         0. +0.j ,  0. +0.j ],
       [ 0. +0.j ,  0. +0.j ,  0.5+0.j , -0.5+0.j ,  0. +0.j ,  0. +0.j ,
         0.5+0.j , -0.5+0.j ]])
    
    def test_is_clifford_accepts(self):
        matrix = self.get_three_qubit_clifford()

        self.assertTrue(cc.is_clifford(matrix))

    def test_is_clifford_rejects_when_columns_have_wrong_relative_phase(self):
        matrix = self.get_three_qubit_clifford()

        matrix[:, 7] *= 1j

        self.assertFalse(cc.is_clifford(matrix))
import unittest
import numpy as np
import clifford.Clifford as cl

class Test_Clifford_Check(unittest.TestCase):
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
    
    def test_get_matrix(self):
        matrix = self.get_three_qubit_clifford()

        clifford = cl.Clifford.from_matrix(matrix)
        new_matrix = clifford.get_matrix()

        self.assertTrue(np.array_equal(matrix, new_matrix))

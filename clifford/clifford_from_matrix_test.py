import unittest
import numpy as np
import math
import clifford.clifford_from_matrix as cc
import pauli.Pauli as p
import pauli.pauli_check as pc

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

    def test_is_clifford_accepts(self):
        matrix = self.get_three_qubit_clifford()
        
        clifford = cc.Clifford_From_Matrix(matrix, only_testing = True)

        self.assertTrue(clifford.is_clifford)

    def test_is_clifford_rejects_when_inconsistent(self):
        matrix = self.get_three_qubit_clifford()
        matrix[:, 7] = matrix [:, 0]

        clifford = cc.Clifford_From_Matrix(matrix, only_testing = True)

        self.assertFalse(clifford.is_clifford)

    def test_is_clifford_rejects_when_remaining_column_not_stabilised(self):
        matrix = self.get_three_qubit_clifford()
        matrix[0,7] = 1

        clifford = cc.Clifford_From_Matrix(matrix, only_testing = True)

        self.assertFalse(clifford.is_clifford)

    def test_is_clifford_rejects_when_initial_column_not_stabilised(self):
        matrix = self.get_three_qubit_clifford()
        matrix[0,1] = 0

        clifford = cc.Clifford_From_Matrix(matrix, only_testing = True)

        self.assertFalse(clifford.is_clifford)

    def test_is_clifford_rejects_when_first_column_not_stabiliser_state(self):
        matrix = self.get_three_qubit_clifford()
        matrix[0,0] = 0

        clifford = cc.Clifford_From_Matrix(matrix, only_testing = True)

        self.assertFalse(clifford.is_clifford)

    def test_is_clifford_rejects_with_global_factor_when_flagged(self):
        matrix = self.get_three_qubit_clifford()
        matrix *= (1+1j)/math.sqrt(2)
        
        clifford = cc.Clifford_From_Matrix(matrix, only_testing = True, allow_global_factor = False)

        self.assertFalse(clifford.is_clifford)

    def test_is_clifford_accepts_with_global_factor(self):
        matrix = self.get_three_qubit_clifford()
        matrix *= (1+1j)/math.sqrt(2)
        
        clifford = cc.Clifford_From_Matrix(matrix, only_testing = True)

        self.assertTrue(clifford.is_clifford)

    def test_is_clifford_rejects_when_reminaing_column_has_wrong_relative_phase(self):
        matrix = self.get_three_qubit_clifford()

        matrix[:, 7] *= 1j

        clifford = cc.Clifford_From_Matrix(matrix, only_testing = True)

        self.assertFalse(clifford.is_clifford)

    def test_is_clifford_rejects_when_intial_column_has_nonstabiliser_relative_phase(self):
        matrix = self.get_three_qubit_clifford()

        matrix[:, 2] *= (1+ 1j)/math.sqrt(2)

        clifford = cc.Clifford_From_Matrix(matrix, only_testing = True)

        self.assertFalse(clifford.is_clifford)

    def test_get_clifford_raises_error_on_invalid_matrix(self):
        matrix = self.get_three_qubit_clifford()
        matrix[:, 7] *= 1j

        clifford = cc.Clifford_From_Matrix(matrix, allow_global_factor = True)

        self.assertRaises(ValueError, clifford.get_clifford)

    def test_get_clifford_without_assuming_clifford(self):
        n = 3
        matrix = self.get_three_qubit_clifford()

        self.assertEqual(matrix.shape[0], 1 << n)

        clifford = cc.Clifford_From_Matrix(matrix).get_clifford()

        for i in range(n):
            z_i = p.Pauli(n, 0, 1 << i, 0, 0).generate_matrix()
            expected_u_i = matrix @ z_i @ matrix.conj().T
            
            self.assertTrue(pc.is_pauli(expected_u_i))

            u_i = clifford.z_conjugates[i].generate_matrix()

            self.assertTrue(np.array_equal(u_i, expected_u_i))
      
        for i in range(n):
          x_i = p.Pauli(n, 1<<i , 0, 0, 0).generate_matrix()
          expected_v_i = matrix @ x_i @ matrix.conj().T
            
          self.assertTrue(pc.is_pauli(expected_v_i))
          
          v_i = clifford.x_conjugates[i].generate_matrix()

          self.assertTrue(np.array_equal(v_i, expected_v_i))

    def test_get_clifford_with_assuming_clifford(self):
        n = 3
        matrix = self.get_three_qubit_clifford()

        self.assertEqual(matrix.shape[0], 1 << n)

        clifford = cc.Clifford_From_Matrix(matrix, assume_clifford = True).get_clifford()

        for i in range(n):
            z_i = p.Pauli(n, 0, 1 << i, 0, 0).generate_matrix()
            expected_u_i = matrix @ z_i @ matrix.conj().T
            
            self.assertTrue(pc.is_pauli(expected_u_i))

            u_i = clifford.z_conjugates[i].generate_matrix()

            self.assertTrue(np.array_equal(u_i, expected_u_i))
      
        for i in range(n):
          x_i = p.Pauli(n, 1<<i , 0, 0, 0).generate_matrix()
          expected_v_i = matrix @ x_i @ matrix.conj().T
            
          self.assertTrue(pc.is_pauli(expected_v_i))

          v_i = clifford.x_conjugates[i].generate_matrix()

          self.assertTrue(np.array_equal(v_i, expected_v_i))

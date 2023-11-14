import unittest
import numpy as np
import math
import clifford.clifford_from_matrix as cc
import pauli.Pauli as p
import pauli.pauli_check as pc

class Test_Clifford_Check(unittest.TestCase):
    unscaled_hadmard = np.array([[1,1],[1,-1]])
    hadamard = (1/math.sqrt(2)) * np.array([[1,1],[1,-1]])

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

    def test_unscaled_multiply_by_hadmard_at_index_0(self):
        test_matrix = np.random.rand(8,8)
        single_hadamard = np.kron(np.eye(4), self.unscaled_hadmard)

        expected_product = test_matrix @ single_hadamard
        cc.unscaled_multiply_by_hadmard_at(test_matrix, 0, 3)

        self.assertTrue(np.allclose(test_matrix, expected_product))

    def test_unscaled_multiply_by_hadmard_at_largest_index(self):
        test_matrix = np.random.rand(8,8)
        single_hadamard = np.kron(self.unscaled_hadmard, np.eye(4))

        expected_product = test_matrix @ single_hadamard
        cc.unscaled_multiply_by_hadmard_at(test_matrix, 2, 3)

        self.assertTrue(np.allclose(test_matrix, expected_product))

    def test_unscaled_multiply_by_hadmard_at_middle_index(self):
        test_matrix = np.random.rand(8,8)
        single_hadamard = np.kron(np.eye(2), np.kron(self.unscaled_hadmard, np.eye(2)) )

        expected_product = test_matrix @ single_hadamard
        cc.unscaled_multiply_by_hadmard_at(test_matrix, 1, 3)

        self.assertTrue(np.allclose(test_matrix, expected_product))

    def test_multiply_by_hadmard_product(self):
        test_matrix = np.random.rand(8,8)
        hadamard_product = np.kron(self.hadamard, np.kron(self.hadamard, self.hadamard) )

        expected_product = test_matrix @ hadamard_product
        matrix_product = cc.multiply_by_hadamard_product(test_matrix, 3)

        self.assertTrue(np.allclose(matrix_product, expected_product))

    # def test_columns_consistent_accepts(self): TODO look at better testing for this
    #     matrix = self.get_three_qubit_clifford()
    #     clifford = cc.
        
    #     self.assertTrue(cc.columns_consistent(matrix, 3, False))

    # def test_columns_consistent_accepts_with_incorrect_relative_phase(self):
    #     matrix = self.get_three_qubit_clifford()
    #     matrix[:,7] *= (1+1j)/math.sqrt(2)
        
    #     self.assertTrue(cc.columns_consistent(matrix, 3, False))

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

    def test_is_clifford_rejects_when_columns_have_wrong_relative_phase(self):
        matrix = self.get_three_qubit_clifford()

        matrix[:, 7] *= 1j

        clifford = cc.Clifford_From_Matrix(matrix, only_testing = True)

        self.assertFalse(clifford.is_clifford)

    def test_get_clifford_raises_error_on_invalid_matrix(self):
        matrix = self.get_three_qubit_clifford()
        matrix[:, 7] *= 1j

        clifford = cc.Clifford_From_Matrix(matrix, allow_global_factor = True)

        self.assertRaises(ValueError, clifford.get_clifford)

    def test_get_clifford(self):
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
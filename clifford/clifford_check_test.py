import unittest
import numpy as np
import math
import clifford.clifford_check as cc

class Test_Clifford_Check(unittest.TestCase):
    unscaled_hadmard = np.array([[1,1],[1,-1]])
    hadamard = (1/math.sqrt(2)) * np.array([[1,1],[1,-1]])
    
    def test_is_full_rank_accepts_when_full_rank(self):
        vectors = [18, 21, 30] # 10010, 10101, 11110

        self.assertTrue(cc.is_full_rank(vectors, 3))

    def test_is_full_rank_rejects_when_lin_dep(self):
        vectors = [45, 26, 10, 61] 

        self.assertTrue(vectors[0] ^ vectors[1] == vectors[2] ^ vectors[3])
        self.assertFalse(cc.is_full_rank(vectors, 4))

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
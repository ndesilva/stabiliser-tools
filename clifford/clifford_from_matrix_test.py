import unittest
import numpy as np
import math
import clifford.clifford_from_matrix as cc
import pauli.Pauli as p
import pauli.pauli_check as pc
import clifford.Clifford as c
import benchmarking.generators as gs

NUM_REPETITIONS = 5

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

    def test_is_clifford_on_random_cliffords(self):
        number_qubits = 6

        for _ in range(NUM_REPETITIONS):
            matrix = gs.random_almost_clifford(number_qubits)
            
            clifford = cc.Clifford_From_Matrix(matrix, only_testing = True)
            clifford2 = cc.Clifford_From_Matrix(matrix)

            self.assertFalse(clifford.is_clifford)
            self.assertFalse(clifford2.is_clifford)

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

        self.assertTrue(extracted_clifford_is_correct(matrix, clifford, n))

    def test_get_clifford_with_assuming_clifford(self):
        n = 3
        matrix = self.get_three_qubit_clifford()

        self.assertEqual(matrix.shape[0], 1 << n)

        clifford = cc.Clifford_From_Matrix(matrix, assume_clifford = True).get_clifford()

        self.assertTrue(extracted_clifford_is_correct(matrix, clifford, n))

    def test_get_clifford_on_random_cliffords(self):
        number_qubits = 6

        for _ in range(NUM_REPETITIONS):
            random_clifford = gs.random_clifford(number_qubits)
            
            clifford = cc.Clifford_From_Matrix(random_clifford, assume_clifford = True).get_clifford()
            clifford2 = cc.Clifford_From_Matrix(random_clifford).get_clifford()

            is_clifford = cc.Clifford_From_Matrix(random_clifford, only_testing = True).is_clifford

            self.assertTrue(extracted_clifford_is_correct(random_clifford, clifford, number_qubits))
            self.assertTrue(extracted_clifford_is_correct(random_clifford, clifford2, number_qubits))
            self.assertTrue(is_clifford)

def extracted_clifford_is_correct(matrix : np.ndarray, clifford : c.Clifford, num_qubits : int) -> bool:
    for i in range(num_qubits):
        z_i = p.Pauli(num_qubits, 0, 1 << i, 0, 0).generate_matrix()
        expected_u_i = matrix @ z_i @ matrix.conj().T

        u_i = clifford.z_conjugates[i].generate_matrix()

        if not np.allclose(u_i, expected_u_i):
            return False
      
    for i in range(num_qubits):
        x_i = p.Pauli(num_qubits, 1<<i , 0, 0, 0).generate_matrix()
        expected_v_i = matrix @ x_i @ matrix.conj().T

        v_i = clifford.x_conjugates[i].generate_matrix()

        if not np.allclose(v_i, expected_v_i):
            return False

    return True

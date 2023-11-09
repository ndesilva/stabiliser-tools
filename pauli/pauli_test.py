import unittest
import numpy as np
from .Pauli import Pauli

class Test_Pauli_Class(unittest.TestCase):
    
    def test_generate_pauli_single_qubit(self):
        pauli = Pauli(1, 1, 1, 1, 1)

        Y = np.array([[0, -1j],[1j, 0]]) # type: ignore

        self.assertTrue(np.array_equal(pauli.generate_matrix(), Y))

    def test_generate_pauli_three_qubits(self):
        pauli = Pauli(3, 6, 5, 1, 0)

        expected_pauli = np.array([
            [0  ,0  ,0  ,0  ,0  ,0  ,1  ,0],
            [0  ,0  ,0  ,0  ,0  ,0  ,0  ,-1],
            [0  ,0  ,0  ,0  ,1  ,0  ,0  ,0],
            [0  ,0  ,0  ,0  ,0  ,-1 ,0  ,0],
            [0  ,0  ,-1 ,0  ,0  ,0  ,0  ,0],
            [0  ,0  ,0  ,1  ,0  ,0  ,0  ,0],
            [-1 ,0  ,0  ,0  ,0  ,0  ,0  ,0],
            [0  ,1  ,0  ,0  ,0  ,0  ,0  ,0]
        ])

        self.assertTrue(np.array_equal(pauli.generate_matrix(), expected_pauli))  

    def test_get_sign_eigenvalue_raises_error_with_different_dimensions(self):
        state_vector = np.array([1,2])

        pauli = Pauli(2, 0, 0, 0, 0)

        self.assertRaises(ValueError, pauli.get_sign_eigenvalue, state_vector)

    def test_get_sign_eigenvalue_gives_0(self):
        pauli = Pauli(4, 12, 5, 1, 1) # X Y 1 Z
        state_vector = np.array([0, 1, 0, 0, 0, -1j, 0, 0, 0, 1, 0, 0, 0, -1j, 0, 0])

        self.assertEqual(pauli.get_sign_eigenvalue(state_vector), 0)

    def test_get_sign_eigenvalue_gives_1(self):
        pauli = Pauli(4, 14, 5, 1, 1) # X Y X Z
        state_vector = np.array([1, 0, -1, 0, 1j, 0, -1j, 0, 1, 0, -1, 0, 1j, 0, -1j, 0])

        self.assertEqual(pauli.get_sign_eigenvalue(state_vector), 1)

    def test_get_sign_eigenvalue_returns_none_when_not_eigenstate(self):
        pauli = Pauli(1, 1, 0, 0, 0) # X
        state_vector = np.array([1,0])

        self.assertIsNone(pauli.get_sign_eigenvalue(state_vector))

    def test_get_sign_eigenvalue_returns_none_when_non_sign_eigenstate(self):
        pauli = Pauli(2, 2, 2, 0, 0) # -iY 1
        state_vector = np.array([1, 0, 1j, 0])

        self.assertTrue(np.array_equal(-1j*state_vector, pauli.generate_matrix()@state_vector))
        self.assertIsNone(pauli.get_sign_eigenvalue(state_vector))

if __name__ == '__main__':
    unittest.main()
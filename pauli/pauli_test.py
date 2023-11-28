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

    def test_get_sign_eigenvalue_returns_even_when_not_eigenstate_when_flagged(self):
        pauli = Pauli(2, 0, 1, 0, 0) # 1 Z
        state_vector = np.array([1, 1, 1, 1])

        self.assertEqual(pauli.get_sign_eigenvalue(state_vector, assume_equation_holds = True), 0)

    def test_multiply_on_the_right_no_overlap(self):
        pauli_one = Pauli(2, 3, 2, 0, 0) # XZ X
        pauli_two = Pauli(2, 1, 2, 1, 1) #(-1)(-i) Z X

        expected_pauli = Pauli(2, 2, 0, 1, 1) # (-1)(-i) X 1

        pauli_one.multiply_by_pauli_on_right(pauli_two)

        self.assertEqual(pauli_one, expected_pauli)

    def test_multiply_on_the_right_overlap(self):
        pauli_one = Pauli(3, 6, 3, 1, 1) #(-1)(-i) X XZ Z
        pauli_two = Pauli(3, 2, 5, 1, 1) #(-1)(-i) Z X  Z

        expected_pauli = Pauli(3, 4, 6, 0, 0) # XZ Z 1

        pauli_one.multiply_by_pauli_on_right(pauli_two)

        self.assertEqual(pauli_one, expected_pauli)

    def test_anticommutes_with_with_anticommuting_paulis(self):
        pauli_one = Pauli(3, 6, 1, 1, 1) #(-1)(-i) X X Z
        pauli_two = Pauli(3, 2, 5, 1, 1) #(-1)(-i) Z X Z

        self.assertEqual(pauli_one.anticommutes_with(pauli_two), 1)

    def test_anticommutes_with_with_commuting_paulis(self):
        pauli_one = Pauli(3, 6, 3, 1, 1) #(-1)(-i) X XZ Z
        pauli_two = Pauli(3, 2, 5, 1, 1) #(-1)(-i) Z X  Z

        self.assertEqual(pauli_one.anticommutes_with(pauli_two), 0)

    def test_commutes_with_with_anticommuting_paulis(self):
        pauli_one = Pauli(3, 6, 1, 1, 1) #(-1)(-i) X X Z
        pauli_two = Pauli(3, 2, 5, 1, 1) #(-1)(-i) Z X Z

        self.assertEqual(pauli_one.commutes_with(pauli_two), 0)

    def test_commutes_with_with_commuting_paulis(self):
        pauli_one = Pauli(3, 6, 3, 1, 1) #(-1)(-i) X XZ Z
        pauli_two = Pauli(3, 2, 5, 1, 1) #(-1)(-i) Z X  Z

        self.assertEqual(pauli_one.commutes_with(pauli_two), 1)

    def test_is_hermitian_returns_true(self):
        pauli = Pauli(3, 6, 3, 1, 1) # X Y Z

        self.assertTrue(pauli.is_hermitian())

    def test_is_hermitian_returns_false(self):
        pauli = Pauli(3, 5, 6, 1, 0) # i Y Z X

        self.assertFalse(pauli.is_hermitian())

if __name__ == '__main__':
    unittest.main()
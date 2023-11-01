import unittest
import numpy as np
from .pauli_class import Pauli

class Test_Pauli_Class(unittest.TestCase):
    def test_generate_pauli_single_qubit(self):
        pauli = Pauli(1, 1, 1, 0, 1)

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
    
    def test_phase_mod2product_returns_1(self):
        x = 26 #binary 11010
        y = 19 #binary 10011

        self.assertEqual(Pauli.phase_mod2product(x,y), 1)

    def test_phase_mod2product_returns_minus_1(self):
        x = 51 #binary 110011
        y = 49 #binary 110101

        self.assertEqual(Pauli.phase_mod2product(x,y), -1)

if __name__ == '__main__':
    unittest.main()
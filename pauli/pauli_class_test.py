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

if __name__ == '__main__':
    unittest.main()
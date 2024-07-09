print("### LAUNCHING PYTHON TESTS ###")

import sys, unittest
from math import sqrt
import numpy as np

PATH_TO_LIBRARY = './build/ninja-multi-vcpkg/cpp/src/Release'

sys.path.append(PATH_TO_LIBRARY)

import fast as fst

print("imported library")

class TestStabiliserStateMethods(unittest.TestCase):
    def test_stabiliser_state_consistency(self):
        XXX = fst.Pauli(3, 7, 0, 0, 0)
        ZZI = fst.Pauli(3, 0, 6, 0, 0)
        ZIZ = fst.Pauli(3, 0, 5, 0, 0)

        pauli_group = [XXX, ZZI, ZIZ]
        check_matrix = fst.Check_Matrix(pauli_group)
        stabiliser_state = fst.Stabiliser_State(check_matrix)

        state_vector = np.array(stabiliser_state.get_state_vector())

        for pauli in pauli_group:
            matrix = np.array(pauli.get_matrix())
            self.assertTrue(np.array_equal(state_vector, matrix@state_vector))

    def test_is_stabiliser_state(self):
        stabiliser_statevector = self.get_uniform_stabiliser_state(3)
        almost_stabiliser_statevector = self.get_non_stabiliser_statevector(3)

        self.assertTrue(fst.is_stabiliser_state(stabiliser_statevector))
        self.assertFalse(fst.is_stabiliser_state(almost_stabiliser_statevector))

    def test_stabiliser_state_from_statevector(self):
        expected_statevector = self.get_uniform_stabiliser_state(3)
        stabiliser_state = fst.stabiliser_state_from_statevector(expected_statevector)

        statevector = np.array(stabiliser_state.get_state_vector())

        self.assertTrue(np.allclose(expected_statevector, statevector))

    def test_almost_stabiliser_state(self):
        almost_stabiliser_statevector = self.get_non_stabiliser_statevector(3)

        with self.assertRaises(ValueError):
            fst.stabiliser_state_from_statevector(almost_stabiliser_statevector)

        # Check doesn't Raise and exception
        fst.stabiliser_state_from_statevector(almost_stabiliser_statevector, assume_valid = True)

        
    def get_uniform_stabiliser_state(self, number_qubits : int):
        support_size = 1 << number_qubits
        return np.ones(support_size, dtype = complex)/sqrt(support_size)

    def get_non_stabiliser_statevector(self, number_qubits : int):
        non_stabiliser_statevector = self.get_uniform_stabiliser_state(number_qubits)
        non_stabiliser_statevector[-1] *= -1
        return non_stabiliser_statevector

unittest.main()
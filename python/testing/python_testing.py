print("### LAUNCHING PYTHON TESTS ###")

import sys, unittest
from math import sqrt
import numpy as np

PATH_TO_LIBRARY = './build/ninja-multi-vcpkg/cpp/src/Release'

sys.path.append(PATH_TO_LIBRARY)

import fast as fst

print("imported library")

# TODO : add the python tests to the workflow?

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

        # Check doesn't raise an exception
        fst.stabiliser_state_from_statevector(almost_stabiliser_statevector, assume_valid = True)

    def test_consistency_again(self):
        stabiliser_statevector = np.array([0, 1, 0, 0, 0, 0, 1, 0]) / np.sqrt(2)
        self.assertTrue(fst.is_stabiliser_state(stabiliser_statevector))
        
        check_matrix = fst.Check_Matrix( fst.stabiliser_state_from_statevector(stabiliser_statevector) )

        intermediate_stab_state = fst.Stabiliser_State( check_matrix )
        output_statevector = np.array( check_matrix.get_state_vector() )
        
        self.assertTrue( np.linalg.norm(stabiliser_statevector - output_statevector) <= 1e-7 )
        
    def get_uniform_stabiliser_state(self, number_qubits : int):
        support_size = 1 << number_qubits
        return np.ones(support_size, dtype = complex)/sqrt(support_size)

    def get_non_stabiliser_statevector(self, number_qubits : int):
        non_stabiliser_statevector = self.get_uniform_stabiliser_state(number_qubits)
        non_stabiliser_statevector[-1] *= -1
        return non_stabiliser_statevector

class TestCliffordMethods(unittest.TestCase):
    def test_clifford_consistency(self):
        X = fst.Pauli(1,1,0,0,0)
        Z = fst.Pauli(1,0,1,0,0)

        clifford = fst.Clifford([X], [Z], -1)
        matrix = np.array(clifford.get_matrix())

        X_matrix = np.array(X.get_matrix())
        Z_matrix = np.array(Z.get_matrix())

        self.assertTrue(np.allclose(X_matrix, matrix@Z_matrix@matrix.conj().T))
        self.assertTrue(np.allclose(Z_matrix, matrix@X_matrix@matrix.conj().T))

    def test_is_clifford_matrix(self):
        matrix = self.get_hadamard_tensor_hadamard()
        almost_clifford = self.get_almost_clifford_matrix()

        self.assertTrue(fst.is_clifford_matrix(matrix))
        self.assertFalse(fst.is_clifford_matrix(almost_clifford))

    def test_clifford_from_matrix(self):
        expected_matrix = np.array(self.get_hadamard_tensor_hadamard())
        clifford = fst.clifford_from_matrix(expected_matrix)

        matrix = np.array(clifford.get_matrix())

        self.assertTrue(np.allclose(expected_matrix, matrix))

    def test_almost_clifford(self):
        almost_hadamard = self.get_almost_clifford_matrix()

        with self.assertRaises(ValueError):
            fst.clifford_from_matrix(almost_hadamard)

        # Check doesn't Raise an exception
        fst.clifford_from_matrix(almost_hadamard, assume_valid = True)

    def get_hadamard_tensor_hadamard(self):
        return [[.5, .5, .5, .5], [.5, -.5, .5, -.5], [.5, .5, -.5, -.5], [.5, -.5, -.5, .5]]
    
    def get_almost_clifford_matrix(self):
        matrix = self.get_hadamard_tensor_hadamard()
        matrix[3][3] *= -1

        return matrix
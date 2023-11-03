import unittest
import numpy as np
from stabiliser_state.Stabiliser_State import Stabiliser_State

class Test_Stabiliser_State_Class(unittest.TestCase):
    
    def test_from_state_vector(self):
        state_vector = np.array([0,-1j, 1, 0, 1j, 0, 0, 1])

        state = Stabiliser_State.from_statevector(state_vector)

        self.assertEqual(state.shift, 1)
        self.assertEqual(state.vector_basis, [3, 5])
        self.assertEqual(state.real_linear_part, 2)
        self.assertEqual(state.imaginary_part, 1)
        self.assertEqual(state.quadratic_form, [3])
        self.assertEqual(state.global_factor, -2j)

    def test_from_state_vector_raises_error_when_not_stab_state(self):
        non_stab_state = (1/np.sqrt(2))*np.array([1, 1/np.sqrt(2)*(1 + 1j)]) # T state is not stabiliser!

        self.assertRaises(ValueError, Stabiliser_State.from_statevector, non_stab_state)

    def test_generate_state_vector_case_one(self):
        number_qubits = 1
        quadratic_form = []
        real_linear_part = 0
        imag_part = 1
        vector_basis = [1]
        shift = 0

        state = Stabiliser_State(number_qubits, quadratic_form, real_linear_part, imag_part, vector_basis, shift)
        expected_vector = 1/(np.sqrt(2))*np.array([1, 1j])

        self.assertTrue(np.array_equal(state.generate_state_vector(), expected_vector))

    def test_generate_state_vector_case_two(self):
        number_qubits = 1
        quadratic_form = []
        real_linear_part = 1
        imag_part = 0
        vector_basis = [0]
        shift = 1

        state = Stabiliser_State(number_qubits, quadratic_form, real_linear_part, imag_part, vector_basis, shift)
        expected_vector = 1/(np.sqrt(2))*np.array([0, -1])

        self.assertTrue(np.array_equal(state.generate_state_vector(), expected_vector))

    def test_generate_state_vector_with_global_factor(self):
        number_qubits = 3
        quadratic_form = [3]
        real_linear_part = 1
        imag_part = 2
        vector_basis = [6,1]
        shift = 4
        global_factor = 1j

        state = Stabiliser_State(number_qubits, quadratic_form, real_linear_part, imag_part, vector_basis, shift, global_factor = global_factor)
        expected_vector = .5*np.array([0,0,-1j,-1,1j,-1,0,0])

        self.assertTrue(np.array_equal(state.generate_state_vector(), expected_vector))

    def test_row_reduce_basis(self):
        number_qubits = 3
        quadratic_form = []
        real_linear_part = 0
        imag_part = 0
        vector_basis = [13, 22, 26]
        shift = 2

        state = Stabiliser_State(number_qubits, quadratic_form, real_linear_part, imag_part, vector_basis, shift)
        state.row_reduce_basis()

        self.assertEqual(state.vector_basis, [22, 12, 1])

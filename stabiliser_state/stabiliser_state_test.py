import math
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
        non_stab_state = (1/math.sqrt(2))*np.array([1, 1/math.sqrt(2)*(1 + 1j)]) # T state is not stabiliser!

        self.assertRaises(ValueError, Stabiliser_State.from_statevector, non_stab_state)

    def test_generate_state_vector_case_one(self):
        number_qubits = 1
        quadratic_form = []
        real_linear_part = 0
        imag_part = 1
        vector_basis = [1]
        shift = 0

        state = Stabiliser_State(number_qubits, quadratic_form, real_linear_part, imag_part, vector_basis, shift)
        expected_vector = 1/(math.sqrt(2))*np.array([1, 1j])

        self.assertTrue(np.array_equal(state.get_state_vector(), expected_vector))

    def test_generate_state_vector_case_two(self):
        number_qubits = 1
        quadratic_form = []
        real_linear_part = 1
        imag_part = 0
        vector_basis = [0]
        shift = 1

        state = Stabiliser_State(number_qubits, quadratic_form, real_linear_part, imag_part, vector_basis, shift)
        expected_vector = 1/(math.sqrt(2))*np.array([0, -1])

        self.assertTrue(np.array_equal(state.get_state_vector(), expected_vector))

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

        self.assertTrue(np.array_equal(state.get_state_vector(), expected_vector))

    def test_get_quadratic_form_as_dictionary(self):
        number_qubits = 3
        quadratic_form = [3]
        real_linear_part = 1
        imag_part = 2
        vector_basis = [4,2,1]
        shift = 0
        global_factor = 1j

        state = Stabiliser_State(number_qubits, quadratic_form, real_linear_part, imag_part, vector_basis, shift, global_factor = global_factor)

        expected_dict = {0:0, 1:1, 2:0, 4:0, 3:1, 5:0, 6:0}

        self.assertEqual(expected_dict, state._Stabiliser_State__get_quadratic_form_as_dictionary())

    def test_update_class_quadratic_form(self):
        number_qubits = 3
        quadratic_form = [3]
        real_linear_part = 1
        imag_part = 2
        vector_basis = [4,2,1]
        shift = 0
        global_factor = 1j

        state = Stabiliser_State(number_qubits, quadratic_form, real_linear_part, imag_part, vector_basis, shift, global_factor = global_factor)

        quadratic_dictionary = {0:0, 1:1, 2:0, 4:1, 3:0, 5:1, 6:1}
        state._Stabiliser_State__set_quadratic_form_from_dict(quadratic_dictionary)

        self.assertEqual(5, state.real_linear_part)
        self.assertEqual([5,6], state.quadratic_form)

    def test_row_reduce_basis(self):
        number_qubits = 5
        quadratic_form = [3]
        real_linear_part = 1
        imag_part = 1
        vector_basis = [13, 22, 26]
        shift = 2
        gloal_factor = math.sqrt(8)

        state = Stabiliser_State(number_qubits, quadratic_form, real_linear_part, imag_part, vector_basis, shift, global_factor = gloal_factor)
        initial_state_vector = state.get_state_vector()
        state._Stabiliser_State__row_reduce_basis()

        self.assertEqual(state.vector_basis, [12, 22, 1])
        self.assertTrue(np.array_equal(initial_state_vector, state.get_state_vector()))

    def test_get_stabiliser_group_generators(self):
        number_qubits = 5
        vector_basis = [16, 9, 2] # already row reduced
        shift = 1
        linear_part = 5
        imag_part = 5
        quadratic_form = [6, 3]
        global_phase = math.sqrt(8)

        state = Stabiliser_State(number_qubits, quadratic_form, linear_part, imag_part, vector_basis, shift, global_factor = global_phase)
        state_vector = state.get_state_vector()
        check_matrix = state.get_check_matrix()

        for pauli in check_matrix.paulis:
            self.assertTrue(np.array_equal(state_vector, pauli.generate_matrix()@state_vector))

        self.assertEqual(len(check_matrix.paulis), 5)

    def test_get_stabiliser_group_generators_basis_not_row_reduced(self):
        number_qubits = 5
        vector_basis = [17, 20 , 30] # 10001, 10100, 11110
        shift = 1
        linear_part = 5
        imag_part = 5
        quadratic_form = [6, 3]
        global_phase = math.sqrt(8)

        state = Stabiliser_State(number_qubits, quadratic_form, linear_part, imag_part, vector_basis, shift, global_factor = global_phase)
        
        state_vector = state.get_state_vector()
        check_matrix = state.get_check_matrix()

        state.vector_basis.sort()
        self.assertEqual(state.vector_basis, [5, 10, 17])

        for pauli in check_matrix.paulis:
            self.assertTrue(np.array_equal(state_vector, pauli.generate_matrix()@state_vector))

        self.assertEqual(len(check_matrix.paulis), 5)
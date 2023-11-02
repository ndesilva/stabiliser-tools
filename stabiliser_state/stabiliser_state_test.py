import unittest
import numpy as np
from .Stabiliser_State import Stabiliser_State

class Test_Stabiliser_State_Class(unittest.TestCase):
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

    def test_generate_state_vector_case_three(self):
        number_qubits = 3
        quadratic_form = [3]
        real_linear_part = 1
        imag_part = 2
        vector_basis = [6,1]
        shift = 4

        state = Stabiliser_State(number_qubits, quadratic_form, real_linear_part, imag_part, vector_basis, shift)
        expected_vector = .5*np.array([0,0,-1,1j,1,1j,0,0])

        self.assertTrue(np.array_equal(state.generate_state_vector(), expected_vector))

import math
import unittest
import numpy as np
import stabiliser_state.stabiliser_from_state_vector as ssv
from stabiliser_state.Stabiliser_State import Stabiliser_State
import benchmarking.generators as gs

NUM_REPETITIONS = 10

class Test_Stabiliser_From_State_Vector(unittest.TestCase):

    three_stab_state = Stabiliser_State(3, [3], 1, 2, [6,1], 4)
    five_stab_state = Stabiliser_State(5, [3,5], 4, 2, [9, 22, 16], 1)
    
    def test_is_stab_state_accepts_basis_state(self):
        state = np.array([0,1])

        self.assertTrue(ssv.Stabiliser_From_State_Vector(state).is_stab_state)

    def test_is_stab_state_accepts_three_qubit_state(self):
        state = self.three_stab_state.get_state_vector()

        self.assertTrue(ssv.Stabiliser_From_State_Vector(state).is_stab_state)

    def test_is_stab_state_rejects_on_all_zero_state(self):
        state = np.zeros(4)

        self.assertFalse(ssv.Stabiliser_From_State_Vector(state).is_stab_state)

    def test_is_stab_state_rejects_with_incorrect_support_size(self):
        state = self.three_stab_state.get_state_vector()
        
        state[0] = 1 # add element to the support

        self.assertFalse(ssv.Stabiliser_From_State_Vector(state).is_stab_state)

    def test_is_stab_state_rejects_when_support_not_affine_space(self):
        state = self.five_stab_state.get_state_vector()

        state[1] = 0 # remove element of affine space
        state[0] = 1 # add non-affine space element to support

        self.assertFalse(ssv.Stabiliser_From_State_Vector(state).is_stab_state)

    def test_is_stab_state_rejects_with_global_factor_when_flagged(self):
        state = np.exp(2j)*self.three_stab_state.get_state_vector()

        self.assertFalse(ssv.Stabiliser_From_State_Vector(state, allow_global_factor = False).is_stab_state)
    
    def test_is_stab_state_accpets_with_global_factor(self):
        state = (1/math.sqrt(2))*(1+1j)*self.five_stab_state.get_state_vector()

        self.assertTrue(ssv.Stabiliser_From_State_Vector(state).is_stab_state)

    def test_is_stab_state_rejects_with_invalid_linear_entry(self):
        state = self.five_stab_state.get_state_vector()
        
        state[7] = 2 # invalid weight 1 term (first basis vector)

        self.assertFalse(ssv.Stabiliser_From_State_Vector(state).is_stab_state)

    def test_is_stab_state_rejects_with_invalid_quadratic_entry(self):
        state = self.five_stab_state.get_state_vector()
        
        state[14] = -2j # invalid weight 2 term (e1 + e2)

        self.assertFalse(ssv.Stabiliser_From_State_Vector(state).is_stab_state)

    def test_is_stab_state_rejects_with_inconsistent_entry(self):
        state = self.five_stab_state.get_state_vector()
        
        state[30] = -1j # change (e1+e2+e3) element (should be 1j)

        self.assertFalse(ssv.Stabiliser_From_State_Vector(state).is_stab_state)

    def test_get_stab_state_returns_correct_state(self):
        state_vector = np.array([0,-1j, 1, 0, 1j, 0, 0, 1])

        state = ssv.Stabiliser_From_State_Vector(state_vector, allow_global_factor = True).get_stab_state()

        self.assertEqual(state.shift, 1)
        self.assertEqual(state.vector_basis, [3, 5])
        self.assertEqual(state.real_linear_part, 2)
        self.assertEqual(state.imaginary_part, 1)
        self.assertEqual(state.quadratic_form, [3])
        self.assertEqual(state.global_factor, -2j)

    def test_get_stab_state_raises_error_on_invalid_state(self):
        state_vector = np.array([1, (1/math.sqrt(2))*(1+1j)])

        state = ssv.Stabiliser_From_State_Vector(state_vector, allow_global_factor = True)

        self.assertRaises(ValueError, state.get_stab_state)

    def test_is_stab_accepts_on_invalid_state_if_flagged(self):
        state_vector = np.array([1, 1, 1, 1, 1, 1, 1, (1/math.sqrt(2))*(1+1j)])/math.sqrt(8)
        
        state = ssv.Stabiliser_From_State_Vector(state_vector, assume_stab_state = True)

        self.assertTrue(state.is_stab_state)

    def test_get_stab_state_on_random_stab_states(self):
        num_qubits = 6
        
        for _ in range(NUM_REPETITIONS):
            vector = gs.random_stab_state(num_qubits)

            state = ssv.Stabiliser_From_State_Vector(vector)
            state2 = ssv.Stabiliser_From_State_Vector(vector, assume_stab_state = True)
            state3 = ssv.Stabiliser_From_State_Vector(vector, assume_stab_state = True, check_support_first = True)

            self.assertTrue(np.allclose(state.get_stab_state().get_state_vector(), vector))
            self.assertTrue(np.allclose(state2.get_stab_state().get_state_vector(), vector))
            self.assertTrue(np.allclose(state3.get_stab_state().get_state_vector(), vector))


    def test_is_stab_state_on_random_almost_stab_states(self):
        num_qubits = 6
        
        for _ in range(NUM_REPETITIONS):
            vector = gs.random_almost_stab_state(num_qubits)

            state = ssv.Stabiliser_From_State_Vector(vector)
            state2 = ssv.Stabiliser_From_State_Vector(vector, check_support_first = True)

            self.assertFalse(state.is_stab_state)
            self.assertFalse(state2.is_stab_state)
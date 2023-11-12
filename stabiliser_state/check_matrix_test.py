import unittest
from pauli.Pauli import Pauli
from stabiliser_state.Check_Matrix import Check_Matrix
import numpy as np

class Test_Check_Matrix(unittest.TestCase):
    
    @staticmethod
    def get_default_five_qubit_check_matrix():
        paulis = []
        
        paulis.append(Pauli(5, 3, 1, 0, 1))
        paulis.append(Pauli(5, 5, 2, 0, 0))
        paulis.append(Pauli(5, 6, 4, 1, 1))
        paulis.append(Pauli(5, 16, 8, 1, 0))
        paulis.append(Pauli(5, 22, 4, 0, 1))

        return Check_Matrix(paulis)
    
    @staticmethod
    def get_default_row_reduced_five_qubit_check_matrix():
        pass
        
        paulis = []
        
        paulis.append(Pauli(5, 3, 1, 0, 1))
        paulis.append(Pauli(5, 5, 2, 0, 0))
        paulis.append(Pauli(5, 0, 7, 0, 0))
        paulis.append(Pauli(5, 16, 8, 1, 0))
        paulis.append(Pauli(5, 0, 8, 0, 0))

        return Check_Matrix(paulis)

    def test_extract_zero_x_paulis(self):
        pauli1 = Pauli(3, 0, 1, 0, 0) # 1 1 Z
        pauli2 = Pauli(3, 2, 0, 0, 0) # 1 X 1
        pauli3 = Pauli(3, 4, 0, 0, 0) # X 1 1

        check_matrix = Check_Matrix([pauli1, pauli2, pauli3])

        self.assertEqual(check_matrix.non_zero_x, [pauli2, pauli3])
        self.assertEqual(check_matrix.zero_x, [pauli1])

    def test_put_into_reduced_form(self):
        check_matrix = self.get_default_five_qubit_check_matrix()
        expected_reduced_check_matrix = self.get_default_row_reduced_five_qubit_check_matrix()

        check_matrix._Check_Matrix__put_into_reduced_form()

        self.assertEqual(check_matrix.paulis, expected_reduced_check_matrix.paulis)
        self.assertEqual(check_matrix.non_zero_x, expected_reduced_check_matrix.non_zero_x)
        self.assertEqual(check_matrix.zero_x, expected_reduced_check_matrix.zero_x)

    def test_get_stabiliser_state(self):
        check_matrix = self.get_default_five_qubit_check_matrix()
        stab_state = check_matrix.get_stabiliser_state()

        state_vector = stab_state.generate_state_vector()

        for pauli in check_matrix.paulis:
            self.assertTrue(np.array_equal(state_vector, pauli.generate_matrix()@state_vector))
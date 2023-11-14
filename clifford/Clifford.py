from __future__ import annotations

import numpy as np
import pauli.Pauli as p
import clifford.clifford_from_matrix as cfm
import stabiliser_state.Check_Matrix as cm
import F2_helper.F2_helper as f2

class Clifford:

    @staticmethod
    def from_matrix(matrix: np.ndarray, assume_clifford : bool = False) -> Clifford: # TODO test
        clifford = cfm.Clifford_From_Matrix(matrix, allow_global_factor = True, assume_clifford = assume_clifford)

        if not clifford.is_clifford:
            raise ValueError('Matrix does not correspond to a Clifford')
        
        return clifford.get_clifford()

    def __init__(self, z_conjugates : list[p.Pauli], x_cojugates : list[p.Pauli], global_phase : complex = 1):
        self.number_qubits = len(z_conjugates)
        
        self.z_conjugates = z_conjugates
        self.x_conjugates = x_cojugates

        self.global_phase = global_phase

    def get_matrix(self) -> np.ndarray:
        size = 1 << self.number_qubits
        matrix = np.zeros((size, size), dtype = complex)
        
        first_col_check_matrix = cm.Check_Matrix(self.z_conjugates)
        matrix[:, 0] = first_col_check_matrix.get_state_vector() * self.global_phase

        x_vectors = [pauli.x_vector for pauli in self.x_conjugates]
        z_vectors = [pauli.z_vector for pauli in self.x_conjugates]
        
        non_zero_indices = np.nonzero(matrix[:, 0])[0] # TODO make o(n) rather than O(N)

        quad_form = []

        for i in range(self.number_qubits):
            beta_i = self.x_conjugates[i].z_vector
            
            for j in range(i):
                if f2.mod2product(beta_i, self.x_conjugates[j].x_vector):
                    quad_form.append( 1<<i | 1<<j )

        for col_index in range(1, size):
            phase = f2.sign_evaluate_poly(quad_form, col_index)
            
            loop_counter = col_index

            while loop_counter: # TODO make more efficient by using F2 and F4?
                left_most_index = f2.fast_log2(loop_counter)
                phase *= self.x_conjugates[left_most_index].phase
                loop_counter ^= (1 << left_most_index) # remove left most one

            x_vector = f2.get_vector_expansion(self.number_qubits, x_vectors, col_index)
            z_vector = f2.get_vector_expansion(self.number_qubits, z_vectors, col_index)

            for index in non_zero_indices:
                matrix[index ^ x_vector, col_index] = phase * f2.sign_mod2product(index, z_vector) * matrix[index, 0]

        return matrix
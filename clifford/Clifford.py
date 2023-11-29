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

        old_col_index = 0

        for i in range(1, size):
            new_col_index = i ^ (i >> 1) # iterate through gray code so that we only flip one bit at a time, see https://www.geeksforgeeks.org/generate-n-bit-gray-codes/
            
            bit_flipped = f2.fast_log2(new_col_index ^ old_col_index)

            matrix[:, new_col_index] = self.x_conjugates[bit_flipped].multiply_vector(matrix[:, old_col_index])

            old_col_index = new_col_index

        return matrix
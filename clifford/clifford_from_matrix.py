import numpy as np
import stabiliser_state.stabiliser_from_state_vector as ssv
import F2_helper.F2_helper as f2
import math
import pauli.Pauli as p
import clifford.Clifford as c

class Clifford_From_Matrix: # TODO currently assumes 2^n to 2^n
    def __init__(self, matrix : np.ndarray, allow_global_factor : bool = False, assume_valid_matrix : bool = False, only_testing : bool = False): # TODO test
        self.number_qubits = f2.fast_log2(matrix.shape[0])
        self.is_clifford = False

        if not self.__set_first_col_stabilisers(matrix, allow_global_factor):
            return
        
        if not assume_valid_matrix:
            if not self.__remaining_columns_consistent(matrix):
                return
            
        if not only_testing:
            self.z_conjugates = self.__get_conjugates()
    
        matrix_times_hadamard = multiply_by_hadamard_product(matrix, self.number_qubits)

        if not self.__set_first_col_stabilisers(matrix_times_hadamard, allow_global_factor):
            return
        
        if not assume_valid_matrix:
            if not self.__remaining_columns_consistent(matrix_times_hadamard):
                return
            
        if not only_testing:
            self.x_conjugates = self.__get_conjugates()
    
        self.is_clifford = True

    def get_clifford(self) -> c.Clifford: # TODO test
        if self.is_clifford:
            return c.Clifford(self.z_conjugates, self.x_conjugates)

        raise ValueError('Matrix does not correspond to a Clifford')

    def __set_first_col_stabilisers(self, matrix : np.ndarray, allow_global_factor : bool) -> bool: #  TODO test
        first_col_state = ssv.Stabiliser_From_State_Vector(matrix[:, 0], allow_global_factor)

        if not first_col_state.is_stab_state:
            return False
        
        first_col_stab_group = first_col_state.get_stab_state().get_check_matrix().paulis

        self.pauli_pattern_pairs : list[tuple[p.Pauli, Pauli_Pattern]] = [(pauli, Pauli_Pattern()) for pauli in first_col_stab_group]

        for pair in self.pauli_pattern_pairs:
            for j in range(self.number_qubits):
                phase_bit = pair[0].get_sign_eigenvalue(matrix[: , 1<<j])
                
                if phase_bit == None:
                    return False
                
                pair[1].string |= phase_bit * (1<<j)
        
        return True

    def __remaining_columns_consistent(self, matrix : np.ndarray) -> bool: # TODO test
        for pair in self.pauli_pattern_pairs:
            for col_index in range(1, 1 << self.number_qubits): # TODO Double checking the weight 1 hamming strings again
                phase_bit = pair[0].get_sign_eigenvalue(matrix[:, col_index])
                
                if phase_bit != f2.mod2product(pair[1].string, col_index):
                    return False
        
        return True
    
    def __get_conjugates(self) -> list[p.Pauli]: # TODO test
        
        for i in range(self.number_qubits):
            pauli = self.pauli_pattern_pairs[i][0]
            pattern = self.pauli_pattern_pairs[i][1].string

            pivot_index = f2.fast_log2(pattern)

            for j in range(self.number_qubits):
                if i != j:
                    if f2.get_bit_at(self.pauli_pattern_pairs[j][1].string, pivot_index):
                        self.pauli_pattern_pairs[j][1].string ^= pattern
                        self.pauli_pattern_pairs[j][0].multiply_by_pauli_on_right(pauli)

        
        self.pauli_pattern_pairs.sort(key = lambda x : x[1].string, reverse = True)
        return [pair[0] for pair in self.pauli_pattern_pairs]

            
class Pauli_Pattern:
    def __init__(self):
        self.string = 0

# TODO is this needed?
# def is_full_rank(vectors : list[int], number_vectors : int) -> bool:
#     # copy the vectors list to not alter it
#     vectors_copy = vectors.copy()

#     # row reduce the vectors to row echelon form, stopping if we ever get an all zero vector
#     for i in range(number_vectors):

#         pivot_index = f2.fast_log2(vectors_copy[i])

#         if pivot_index == -1:
#             return False
        
#         for j in range(i+1, number_vectors): # only need row echelon form (not reduced) to check LI
#             vectors_copy[j] ^= (i!=j) * f2.get_bit_at(vectors_copy[j], pivot_index) * vectors_copy[i]

#     return True

def multiply_by_hadamard_product(matrix : np.ndarray, number_qubits : int) -> np.ndarray:
    for j in range(number_qubits):
        unscaled_multiply_by_hadmard_at(matrix, j, number_qubits)

    return matrix/math.sqrt((1 << number_qubits))

def unscaled_multiply_by_hadmard_at(matrix : np.ndarray, hadamard_index : int, number_qubits : int) -> None:
    for tail in range(1 << hadamard_index):
        for head in range(1 << (number_qubits - hadamard_index - 1)):
            shifted_head = head << (hadamard_index + 1)
            
            first_col_index = shifted_head | tail
            second_col_index = shifted_head | (1 << hadamard_index) | tail

            matrix[:, first_col_index] += matrix[:, second_col_index]
            matrix[:, second_col_index] = matrix[:, first_col_index] - 2*matrix[:, second_col_index] # since first column is set first, this is now the difference
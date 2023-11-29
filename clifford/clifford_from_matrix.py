from __future__ import annotations

import numpy as np
import stabiliser_state.stabiliser_from_state_vector as ssv
import F2_helper.F2_helper as f2
import pauli.Pauli as p
import clifford.Clifford as c
import numba

class Clifford_From_Matrix: # TODO currently assumes 2^n to 2^n unitary
    def __init__(self, matrix : np.ndarray, allow_global_factor : bool = True, assume_clifford : bool = False, only_testing : bool = False): # TODO test flags somehow?
        self.number_qubits = f2.fast_log2(matrix.shape[0])
        self.is_clifford = False

        if not self.__set_first_col_stabilisers(matrix, allow_global_factor, assume_clifford = assume_clifford):
            return
            
        self.__get_z_conjugates()

        if not assume_clifford:
            if not self.__remaining_columns_consistent(matrix):
                return
        
        self.__set_L_matrix()
        self.__set_W_paulis()
        
        if not self.__set_W_pauli_patterns(matrix):
            return
        
        self.__set_x_conjugates() # TODO can we avoid doing this inversion if we are just checking?

        if not assume_clifford:
            if not (self.__x_conjugates_hermitian_and_commute() and self.__remaining_rows_consistent(matrix)):
                return 
        
        self.is_clifford = True

    def get_clifford(self) -> c.Clifford: # TODO test
        if self.is_clifford:
            return c.Clifford(self.z_conjugates, self.x_conjugates, global_phase = self.global_factor)

        raise ValueError('Matrix does not correspond to a Clifford')

    def __set_first_col_stabilisers(self, matrix : np.ndarray, allow_global_factor : bool, assume_clifford : bool) -> bool: #  TODO test
        first_col_state = ssv.Stabiliser_From_State_Vector(matrix[:, 0], allow_global_factor, assume_stab_state = assume_clifford)

        if not first_col_state.is_stab_state:
            return False
       
        self.global_factor = first_col_state.global_factor
        first_col_stab_state = first_col_state.get_stab_state()
        self.shift = first_col_stab_state.shift

        first_col_check_matrix = first_col_stab_state.get_check_matrix()
        
        first_col_stab_group = first_col_check_matrix.paulis
        self.check_matrix_pivots = first_col_check_matrix.pivot_indices

        self.pauli_pattern_pairs : list[tuple[p.Pauli, Pauli_Pattern, int]] = [(first_col_stab_group[j], Pauli_Pattern(), j) for j in range(self.number_qubits)]

        return self.__set_pauli_patterns(matrix, assume_clifford)
    
    def __set_pauli_patterns(self, matrix : np.ndarray, assume_clifford : bool) -> bool:
        for j in range(self.number_qubits):
            non_zero_index = get_first_nonzero_index(matrix[:, 1 << j])
            non_zero_val = matrix[non_zero_index, 1 << j]

            for pair in self.pauli_pattern_pairs:
                phase = matrix[non_zero_index ^ pair[0].x_vector, 1 << j]/(non_zero_val * f2.sign_mod2product(non_zero_index, pair[0].z_vector) * pair[0].phase)
                
                match round_to_5dp(phase):
                    case 1:
                        pass
                    case -1:
                        pair[1].string |= (1<<j)
                    case _:
                        return False
        
        return True

    def __remaining_columns_consistent(self, matrix : np.ndarray) -> bool: # TODO test and numbaify
        for col_index in range(1, 1 << self.number_qubits):
            col = matrix[:, col_index]

            for i in range(self.number_qubits):
                if not self.z_conjugates[i].check_sign_eigenstate(col, self.number_qubits, f2.get_bit_at(col_index, i)):
                    return False

        return True
    
    def __remaining_rows_consistent(self, matrix : np.ndarray) -> bool:        
        old_col_index = 0
        old_support = self.shift
        
        for col_index in range(1, 1 << self.number_qubits): #TODO double checking weight 2 strings   
            new_col_index = col_index ^ (col_index >> 1) # iterate through gray code so that we only flip one bit at a time, see https://www.geeksforgeeks.org/generate-n-bit-gray-codes/
            
            bit_flipped = f2.fast_log2(new_col_index ^ old_col_index)
            pauli_flip = self.x_conjugates[bit_flipped]

            new_support = old_support ^ pauli_flip.x_vector

            if matrix[old_support, old_col_index] * f2.sign_mod2product(old_support, pauli_flip.z_vector) * pauli_flip.phase != matrix[new_support, new_col_index]:
                return False
            
            old_col_index = new_col_index
            old_support = new_support

        return True

    def __x_conjugates_hermitian_and_commute(self) -> bool:
        for i in range(self.number_qubits):
            pauli = self.x_conjugates[i]

            for j in range(i+1, self.number_qubits):
                if pauli.anticommutes_with(self.x_conjugates[j]):
                    return False
                
            if not pauli.is_hermitian():
                return False

        return True

    def __get_z_conjugates(self) -> None: # TODO test
        
        self.additions : list[tuple[int,int]] = []

        for i in range(self.number_qubits):
            pauli = self.pauli_pattern_pairs[i][0]
            pattern = self.pauli_pattern_pairs[i][1].string

            pivot_index = f2.fast_log2(pattern)

            for j in range(self.number_qubits):
                if i != j and f2.get_bit_at(self.pauli_pattern_pairs[j][1].string, pivot_index):
                    self.pauli_pattern_pairs[j][1].string ^= pattern
                    self.pauli_pattern_pairs[j][0].multiply_by_pauli_on_right(pauli)

                    self.additions.append( (i , j) )

        self.pauli_pattern_pairs.sort(key = lambda x : x[1].string)
        self.z_conjugates = [pair[0] for pair in self.pauli_pattern_pairs]
    
    def __set_L_matrix(self):
        self.L = [0]*self.number_qubits

        for j in range(self.number_qubits):
            self.L[self.pauli_pattern_pairs[j][2]] = 1 << (self.number_qubits - j - 1)

        for addition in reversed(self.additions):
            self.L[addition[1]] ^= self.L[addition[0]]

    def __set_W_paulis(self):
        self.w_pauli_tuples = [(p.Pauli(self.number_qubits, 0, 0, 0, 0), Pauli_Pattern()) for _ in range(self.number_qubits)]
        
        non_pivots = [index for index in range(self.number_qubits) if index not in self.check_matrix_pivots]
        
        for i in range(self.number_qubits):
            index = 0
            pauli = self.w_pauli_tuples[self.number_qubits - 1 - i][0]

            for non_pivot in non_pivots:
                pauli.x_vector |= (1 << non_pivot) * f2.get_bit_at(self.L[index], i)
                index += 1
            
            for pivot in self.check_matrix_pivots:
                pauli.z_vector |= (1 << pivot)* f2.get_bit_at(self.L[index], i)
                index += 1

    def __set_W_pauli_patterns(self, matrix) -> bool:
        first_entry = matrix[self.shift, 0]
        
        for i in range(self.number_qubits):
            pauli = self.w_pauli_tuples[i][0]
        
            index = self.shift ^ pauli.x_vector
            
            phase = matrix[index, 1 << i] / (first_entry * f2.sign_mod2product(pauli.z_vector, self.shift))
        
            match round_to_5dp(phase):
                case 1:
                    pass
                case -1:
                    pauli.sign_bit = 1
                case -1j:
                    pauli.i_bit = 1
                case 1j:
                    pauli.sign_bit = 1
                    pauli.i_bit = 1
                case _:
                    return False

            pauli.update_phase()

            for j in range(self.number_qubits):
                j_row_index = self.shift ^ self.w_pauli_tuples[j][0].x_vector
                i_row_index = j_row_index ^ pauli.x_vector

                col_phase = matrix[i_row_index, 1 << i ^ 1<< j]/(matrix[j_row_index , 1 << j] * pauli.phase * f2.sign_mod2product(j_row_index, pauli.z_vector))

                match round_to_5dp(col_phase):
                    case 1:
                        pass
                    case -1:
                        self.w_pauli_tuples[i][1].string |= 1 << j
                    case _:
                        return False

        return True

    def __set_x_conjugates(self):
        self.x_conjugates : list[p.Pauli] = []

        for pair in self.w_pauli_tuples:
            while pair[1].string:
                u_index = f2.fast_log2(pair[1].string)
                pair[1].string ^= (1 << u_index)
                
                pair[0].multiply_by_pauli_on_right(self.z_conjugates[u_index])

            self.x_conjugates.append(pair[0])
      
class Pauli_Pattern:
    def __init__(self):
        self.string = 0

@numba.njit() # TODO put in seperate file, shared with Pauli class
def round_to_5dp(value : complex) -> complex:
    return round(value.real, 5) + 1j*round(value.imag, 5)

@numba.njit()
def get_first_nonzero_index(state_vector : np.ndarray) -> int:
    index = 0

    while not state_vector[index]:  
        index += 1
    
    return index
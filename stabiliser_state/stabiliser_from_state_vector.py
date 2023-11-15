from __future__ import annotations

import math
import numpy as np
from operator import itemgetter
import F2_helper.F2_helper as f2
import stabiliser_state.Stabiliser_State as ss

class Stabiliser_From_State_Vector: #TODO currently assumes length 2^n
    def __init__(self, state_vector : np.ndarray, allow_global_factor : bool = True, assume_stab_state : bool = False, check_support_first : bool = False):
        self.is_stab_state = False
        
        nonzero_indices = np.nonzero(state_vector)[0]
        self.support_size = len(nonzero_indices)

        if self.support_size == 0:
            return
        
        self.dimension = f2.fast_log2(self.support_size)
        self.n = f2.fast_log2(len(state_vector))

        if not self.__support_is_power_two():
            return

        vector_space_value_pairs = self.__get_vector_space_indices_and_amplitudes(nonzero_indices, state_vector)
        
        weight_one_bitstrings = [1<<j for j in range(self.dimension)]

        self.basis_vectors = [int(vector_space_value_pairs[index][0]) for index in weight_one_bitstrings]

        if not assume_stab_state and check_support_first:
            if not self.__vector_space_consistent(vector_space_value_pairs):
                return
            
        non_zero_coeffs = [pair[1] for pair in vector_space_value_pairs]
        self.first_entry = non_zero_coeffs[0]

        if not assume_stab_state:
            if self.__first_entry_not_valid(allow_global_factor):
                return
            
        if not self.__set_linear_parts(weight_one_bitstrings, non_zero_coeffs):
            return
            
        if not self.__set_quadratic_part(non_zero_coeffs):
            return

        if not assume_stab_state:
            if not self.__coefficients_consistent(non_zero_coeffs):
                return
        
        self.global_factor = self.first_entry*(math.sqrt(self.support_size))

        self.is_stab_state = True
    
    def get_stab_state(self) -> ss.Stabiliser_State:
        if self.is_stab_state:
            return ss.Stabiliser_State(self.n, self.quadratic_real_part, self.linear_real_part, self.imag_part, self.basis_vectors, self.shift, global_factor=self.global_factor, row_reduced = True)
        
        raise ValueError('State is not a stabiliser state')
    
    def __support_is_power_two(self,) -> bool:
        return 1 << self.dimension == self.support_size
    
    def __first_entry_not_valid(self, allow_global_factor) -> bool:
        return not (allow_global_factor or is_valid_stabiliser_entry(self.first_entry*(math.sqrt(self.support_size))))

    def __get_vector_space_indices_and_amplitudes(self, nonzero_indices : list[int], state_vector : np.ndarray) -> list[tuple[int, complex]]:
        self.shift = nonzero_indices[0]

        # shift affine space to a vector space, remembering the pairing of the indicies and entries in the state vector
        vector_space_value_pairs = [(index ^ self.shift, state_vector[index]) for index in nonzero_indices]
        
        # sort the vector space in increasing F_2 order; can then apply lemma on sorted vector spaces
        vector_space_value_pairs.sort(key = itemgetter(0)) # itemgetter is faster than declaring a lambda function lambda x : x[0] according to the internet

        return vector_space_value_pairs
    
    def __set_linear_parts(self, weight_one_bitstrings : list[int], non_zero_coeffs : list[complex]) -> bool:
        self.linear_real_part = 0
        self.imag_part = 0

        # get linear terms
        for index in weight_one_bitstrings:
            match non_zero_coeffs[index]/self.first_entry:
                case 1:
                    pass
                case -1:
                    self.linear_real_part |= index
                case 1j:
                    self.imag_part |= index
                case -1j:
                    self.linear_real_part |= index
                    self.imag_part |= index
                case _:
                    #print('invalid linear term')
                    return False
        
        return True
    
    def __set_quadratic_part(self, non_zero_coeffs : list[complex]) -> bool:
        self.quadratic_real_part = []

        for a in range(self.dimension-1):
            for b in range(a+1, self.dimension):
                index = (1 << a) | (1 << b)
                
                linear_part = f2.sign_mod2product(index, self.linear_real_part)*f2.imag_mod2product(index, self.imag_part)

                value = non_zero_coeffs[index]/(self.first_entry*linear_part)

                match value:
                    case 1:
                        pass
                    case -1:
                        self.quadratic_real_part.append(index)
                    case _:
                        #print('invalid quadratic term')
                        return False
        
        return True

    def __vector_space_consistent(self, vector_space_value_pairs : list[tuple[int, complex]]) -> bool:
        # Using lemma, check that the indicies form an F2 vector space. If the support is the whole space, this is not needed
        if self.dimension != self.n:
                for j in range(1<<self.dimension): # we check the basis vectors again here, fast way to not do that?
                    value = f2.get_vector_expansion(self.dimension, self.basis_vectors, j)

                    if vector_space_value_pairs[j][0] != value:
                        #print('Support not affine space')
                        self.is_stab_state = False
                        return False
                    
        return True
    
    def __coefficients_consistent(self, non_zero_coeffs : list[complex]) -> bool:
        
        for index in range(1<<self.dimension): # TODO We are repeating columns of Hamming weight 1,2 - fast way to not do this?
            value = self.first_entry*f2.imag_mod2product(index, self.imag_part)*f2.sign_mod2product(index, self.linear_real_part)*f2.sign_evaluate_poly(self.quadratic_real_part, index)

            if non_zero_coeffs[index] != value:
                #print('inconsistent remainder')
                return False
            
        return True

def is_valid_stabiliser_entry(entry : complex) -> bool:
    match entry:
        case 1:
            return True
        case -1:
            return True
        case 1j:
            return True
        case -1j:
            return True
        case _:
            return False
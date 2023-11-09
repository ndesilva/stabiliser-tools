from __future__ import annotations

import math
import numpy as np
from operator import itemgetter
import F2_helper.F2_helper as f2
import stabiliser_state.Stabiliser_State as ss

class Stabiliser_Checker:
    def __init__(self):
        self.valid_state = False

    def load_vector(self, state_vector : np.ndarray, allow_global_factor = False) -> Stabiliser_Checker:
        nonzero_indices = np.nonzero(state_vector)[0]
        support_size = len(nonzero_indices)
        
        dimension = f2.fast_log2(support_size)
        self.n = f2.fast_log2(len(state_vector))

        # check support is a power of 2
        if 1 << dimension != support_size:
            #print('support not power of 2')
            self.valid_state = False
            return self

        self.shift = nonzero_indices[0]

        # shift affine space to a vector space, remembering the pairing of the indicies and entries in the state vector
        vector_space_value_pairs = [(index ^ self.shift, state_vector[index]) for index in nonzero_indices]
        
        # sort the vector space in increasing F_2 order; can then apply lemma on sorted vector spaces
        vector_space_value_pairs.sort(key = itemgetter(0)) # itemgetter is faster than declaring a lambda function lambda x : x[0] according to the internet

        vector_space_indicies = [pair[0] for pair in vector_space_value_pairs]
        
        weight_one_bitstrings = [1<<j for j in range(dimension)]
        self.basis_vectors = [vector_space_indicies[index] for index in weight_one_bitstrings]

        # Using lemma, check that the indicies form an F2 vector space. If the support is the whole space, this is not needed
        if dimension != self.n:
            for j in range(1<<dimension): # we check the basis vectors again here, fast way to not do that?
                value = f2.get_vector_expansion(dimension, self.basis_vectors, j)

                if vector_space_indicies[j] != value:
                    #print('Support not affine space')
                    self.valid_state = False
                    return self
        
        non_zero_coeffs = [pair[1] for pair in vector_space_value_pairs]
        first_entry = non_zero_coeffs[0]

        if not (allow_global_factor or is_valid_stabiliser_entry(first_entry*(math.sqrt(support_size)))):
            #print('invalid first entry')
            self.valid_state = False
            return self

        self.linear_real_part = 0
        self.imag_part = 0

        # get linear terms
        for index in weight_one_bitstrings:
            match non_zero_coeffs[index]/first_entry:
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
                    self.valid_state = False
                    return self
        
        self.quadratic_real_part = []

        # get quadratic terms:
        for a in range(dimension-1):
            for b in range(a+1, dimension):
                index = (1 << a) | (1 << b)
                
                linear_part = f2.sign_mod2product(index, self.linear_real_part)*f2.imag_mod2product(index, self.imag_part)

                value = non_zero_coeffs[index]/(first_entry*linear_part)

                match value:
                    case 1:
                        pass
                    case -1:
                        self.quadratic_real_part.append(index)
                    case _:
                        #print('invalid quadratic term')
                        self.valid_state = False
                        return self

        for index in range(1<<dimension): # We are repeating columns of Hamming weight 1,2 - fast way to not do this?
            value = first_entry*f2.imag_mod2product(index, self.imag_part)*f2.sign_mod2product(index, self.linear_real_part)*f2.sign_evaluate_poly(self.quadratic_real_part, index)

            if non_zero_coeffs[index] != value:
                #print('inconsistent remainder')
                self.valid_state = False
                return self
        
        self.global_factor = first_entry*(math.sqrt(support_size))

        self.valid_state = True
        return self
    
    def get_stab_state(self) -> ss.Stabiliser_State:
        if self.valid_state:
            return ss.Stabiliser_State(self.n, self.quadratic_real_part, self.linear_real_part, self.imag_part, self.basis_vectors, self.shift, global_factor=self.global_factor)
        
        raise UnboundLocalError('Stabiliser state has not been loaded or attempted to load an invalid state')


def is_valid_stabiliser_entry(entry : float) -> bool:
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
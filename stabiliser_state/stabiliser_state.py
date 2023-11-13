from __future__ import annotations
import math

import numpy as np
import stabiliser_state.stabiliser_from_state_vector as ssv
import F2_helper.F2_helper as f2
from pauli.Pauli import Pauli
import stabiliser_state.Check_Matrix as cm

class Stabiliser_State():

    @staticmethod
    def from_check_matrix(check_matrix : cm.Check_Matrix) -> Stabiliser_State: # TODO test
        return check_matrix.get_stabiliser_state()
    
    @staticmethod
    def from_statevector(state_vector : np.ndarray, assume_stab_state : bool = False) -> Stabiliser_State:
        state = ssv.Stabiliser_From_State_Vector(state_vector, allow_global_factor = True, assume_stab_state = assume_stab_state)

        if not state.is_stab_state:
            raise ValueError('State vector does not describe a stabiliser state')

        return state.get_stab_state()

    def __init__(self, number_qubits: int, quadratic_form : list[int], real_linear_part : int, imaginary_part : int, vector_basis : list[int], shift : int, global_factor : complex = 1, row_reduced = False):
        self.number_qubits = number_qubits
        
        self.quadratic_form = quadratic_form
        self.real_linear_part = real_linear_part
        self.imaginary_part = imaginary_part
        
        self.vector_basis = vector_basis
        self.shift = shift

        self.global_factor = global_factor
        self.row_reduced = row_reduced

        self.dimension = len(self.vector_basis)
    
    def get_state_vector(self) -> np.ndarray:
        size = 1 << self.number_qubits
        
        state_vector = np.zeros(size, dtype=complex)
        dimension = len(self.vector_basis)
        normalisation = self.global_factor/math.sqrt(1 << dimension)

        for j in range(1 << dimension):
            index = f2.get_vector_expansion(dimension, self.vector_basis, j) ^ self.shift
            phase = self.__get_phase(j)
            state_vector[index] = normalisation*phase
        
        return state_vector
    
    def get_check_matrix(self) -> cm.Check_Matrix: 
        # needed for finding the basis of the null space
        self.__row_reduce_basis()
        
        pauli_group = []
        
        # pivot column indices, as vector_basis is now in reduced row echelon form
        pivot_indicies = [f2.fast_log2(vector) for vector in self.vector_basis]
        
        self.__add_pure_z_stabilisers(pauli_group, pivot_indicies)

        # now do X-type stabilisers
        self.__add_x_type_stabilisers(pauli_group, pivot_indicies)
        
        return cm.Check_Matrix(pauli_group, reduced_form = True)
    
    def __get_phase(self, afffine_space_index : int) -> complex:
        return f2.sign_evaluate_poly(self.quadratic_form, afffine_space_index)*f2.sign_mod2product(self.real_linear_part, afffine_space_index)*f2.imag_mod2product(self.imaginary_part, afffine_space_index)
    
    def __row_reduce_basis(self): # TODO test self.row_reduced flag is working
        if not self.row_reduced:
            quadratic_dictionary = self.__get_quadratic_form_as_dictionary()
            self.__do_row_reduction(quadratic_dictionary) # Note that as we change the basis, we also need to change the quadratic form
            self.__set_quadratic_form_from_dict(quadratic_dictionary)

        self.row_reduced = True

    def __do_row_reduction(self, quadratic_dictionary):
        for j in range(self.dimension):
            pivot_row = self.vector_basis[j]
            pivot_index = f2.fast_log2(pivot_row)

            for i in range(self.dimension):
                if f2.get_bit_at(self.vector_basis[i], pivot_index) and i != j:
                    self.__add_ej_to_ei(quadratic_dictionary, j, i)

    def __add_ej_to_ei(self, quadratic_dictionary, j, i):
        self.vector_basis[i] ^= self.vector_basis[j]
                    
        self.imaginary_part ^= (1<<i)*f2.get_bit_at(self.imaginary_part, j)

        for k in range(self.dimension):
            quadratic_dictionary[(1 << k | 1 << i)] ^= quadratic_dictionary[( 1 << k ^ 1 << j)] # rather than deal with k = j case, just put that in the 0 bin

        quadratic_dictionary[ 1<<i ] ^= quadratic_dictionary[ 1<<j ]

    def __get_quadratic_form_as_dictionary(self) -> dict[int, int]:
        quadratic_dictionary = {(1<<j)|(1<<i) : 0 for i in range(self.dimension - 1) for j in range(i+1, self.dimension)} 
        quadratic_dictionary[0] = 0

        for coeff in self.quadratic_form:
            quadratic_dictionary[coeff] = 1

        for i in range(self.dimension):
            quadratic_dictionary[1<<i] = f2.get_bit_at(self.real_linear_part, i)

        return quadratic_dictionary
    
    def __set_quadratic_form_from_dict(self, quadratic_dictionary):
        self.real_linear_part = 0

        for i in range(self.dimension):
            self.real_linear_part |= (1<<i) * quadratic_dictionary[1<<i]
            quadratic_dictionary[ 1<<i ] = 0

        self.quadratic_form = []
        
        for index, included in quadratic_dictionary.items():
            if included:
                self.quadratic_form.append(index)

    def __add_pure_z_stabilisers(self, pauli_group, pivot_indicies):
        for j in range(self.number_qubits):
            if j in pivot_indicies:
                pass
            
            else:
                alpha = 1 << j
                
                # add a one in every pivot column where the corresponding row has a one in column j (to make alpha orthogonal to all the rows)
                for l in range(self.dimension):
                    alpha |= (f2.get_bit_at(self.vector_basis[l], j)) * (1 << pivot_indicies[l])

                # sign bit fixed by c
                sign_bit = f2.mod2product(alpha, self.shift)

                pauli_group.append(Pauli(self.number_qubits, 0, alpha, sign_bit, 0))

    def __add_x_type_stabilisers(self, pauli_group, pivot_indicies):
        for i in range(self.dimension):
            imag_bit = f2.get_bit_at(self.imaginary_part, i)
            beta_vector = 0

            # Cursed for loop, as discussed with Ming
            for j in range(self.dimension): # faster way to do this? set membership query is O(1) as opposed to O(n)
                beta_vector ^= (1<<pivot_indicies[j]) * ( ( (1 <<i | 1 <<j) in self.quadratic_form ) ^ imag_bit* f2.get_bit_at(self.imaginary_part,j) )

            sign_bit = f2.get_bit_at(self.real_linear_part, i) ^ imag_bit ^ f2.mod2product(beta_vector, self.shift)
            
            pauli_group.append(Pauli(self.number_qubits, self.vector_basis[i], beta_vector, sign_bit, imag_bit))

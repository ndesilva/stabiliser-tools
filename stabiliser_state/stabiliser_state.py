import numpy as np
import functools
import stabiliser_state.stabiliser_check as sc
import F2_helper.F2_helper as f2
from pauli.Pauli import Pauli

class Stabiliser_State():

    @staticmethod
    def from_statevector(state_vector : np.ndarray): # TODO test that this works
        state = sc.is_stabiliser_state(state_vector, return_state = True, allow_global_factor = True)

        if not state:
            raise ValueError('State vector does not describe a stabiliser state')

        return state

    def __init__(self, number_qubits: int, quadratic_form : list[int], real_linear_part : int, imaginary_part : int, vector_basis : list[int], shift : int, global_factor : complex = 1):
        self.number_qubits = number_qubits
        
        self.quadratic_form = quadratic_form
        self.real_linear_part = real_linear_part
        self.imaginary_part = imaginary_part
        
        self.vector_basis = vector_basis
        self.shift = shift

        self.global_factor = global_factor

        self.dimension = len(self.vector_basis)
    
    def generate_state_vector(self) -> np.ndarray:
        size = 1 << self.number_qubits
        
        state_vector = np.zeros(size, dtype=complex)
        dimension = len(self.vector_basis)
        normalisation = self.global_factor/np.sqrt(1 << dimension)

        for j in range(1 << dimension):
            index = self.get_state_vector_index(j)
            phase = self.get_phase(j)
            state_vector[index] = normalisation*phase
        
        return state_vector
    
    def get_state_vector_index(self, affine_space_index : int) -> int:
        vectors = [self.vector_basis[index]*(f2.get_bit_at(affine_space_index, index)) for index in range(self.dimension)]
        return functools.reduce(lambda x,y : x^y, vectors) ^ self.shift
    
    def get_phase(self, afffine_space_index : int) -> complex:
        return f2.sign_evaluate_poly(self.quadratic_form, afffine_space_index)*f2.sign_mod2product(self.real_linear_part, afffine_space_index)*f2.imag_mod2product(self.imaginary_part, afffine_space_index)
    
    # WARNING: assumes that the basis vectors are in reduced row echelon form (as will be the case if they were generated from the state_vector) - TODO change this?
    def get_stabiliser_group_generators(self) -> list[Pauli]:  # TODO test this & refactor     
        pauli_group = []
        
        # pivot column indices, assumig that the vector_basis is already in reduced row echelon form
        pivot_indicies = [f2.num_bin_digits(vector) for vector in self.vector_basis]
        
        # get basis of null space of matrix with basis vectors as rows by iterating through the non-pivot columns. These correspond to Z-type stabilisers
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

        # now do X-type stabilisers
        for i in range(self.dimension):
            
            imag_bit = f2.get_bit_at(self.imaginary_part, i)

            # add l_i * l_j part
            beta = imag_bit*self.imaginary_part

            for j in range(self.number_qubits): #faster way to do this?
                beta ^= (1<<j)*( (1<<i | 1<<j) in self.quadratic_form) # add Q + Q^t part

            sign_bit = f2.get_bit_at(self.real_linear_part, i) ^ imag_bit ^ f2.mod2product(beta, self.shift)
            
            pauli_group.append(Pauli(self.number_qubits, self.vector_basis[i], beta, sign_bit, imag_bit))
        
        return pauli_group
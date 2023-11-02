import numpy as np
import functools
from F2_helper.F2_helper import sign_mod2product, imag_mod2product, sign_evaluate_poly

class Stabiliser_State():
    def __init__(self, number_qubits: int, quadratic_form : list[int], real_linear_part : int, imaginary_part : int, vector_basis : list[int], shift : int):
        self.number_qubits = number_qubits
        
        self.quadratic_form = quadratic_form
        self.real_linear_part = real_linear_part
        self.imaginary_part = imaginary_part
        
        self.vector_basis = vector_basis
        self.shift = shift

        self.dimension = len(self.vector_basis)
    
    def generate_state_vector(self) -> np.ndarray:
        size = 1 << self.number_qubits
        
        state_vector = np.zeros(size, dtype=complex)
        dimension = len(self.vector_basis)
        normalisation = 1/np.sqrt(1 << dimension)

        for j in range(1 << dimension):
            index = self.get_state_vector_index(j)
            phase = self.get_phase(j)
            state_vector[index] = normalisation*phase
        
        return state_vector
    
    def get_state_vector_index(self, affine_space_index : int) -> int:
        vectors = [self.vector_basis[index]*(affine_space_index & (1<<index) == (1 << index)) for index in range(self.dimension)]
        return functools.reduce(lambda x,y : x^y, vectors) ^ self.shift
    
    def get_phase(self, afffine_space_index : int) -> complex:
        return sign_evaluate_poly(self.quadratic_form, afffine_space_index)*sign_mod2product(self.real_linear_part, afffine_space_index)*imag_mod2product(self.imaginary_part, afffine_space_index)
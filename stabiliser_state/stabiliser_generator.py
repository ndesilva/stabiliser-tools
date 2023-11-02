import numpy as np
import functools
from F2_helper.F2_helper import evaluate_poly

class Stabiliser_Generator():
    def __init__(self, number_qubits: int, quadratic_form : list[int], imaginary_part : list[int], vector_basis : list[int], shift : int):
        self.number_qubits = number_qubits
        self.quadratic_form = quadratic_form
        self.imaginary_part = imaginary_part
        self.vector_basis = vector_basis
        self.shift = shift

        self.dimension = len(self.vector_basis)
    
    def get_state_vector(self) -> np.ndarray:
        size = 1 << self.number_qubits
        
        state_vector = np.zeros(size)
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
        return (1-2*evaluate_poly(self.quadratic_form, afffine_space_index))*(1+(1j-1)*evaluate_poly(self.imaginary_part, afffine_space_index))
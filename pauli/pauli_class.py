import numpy as np
from F2_helper.F2_helper import phase_mod2product
class Pauli:
    def __init__(self, number_qubits : int, x_vector : int, y_vector : int, sign_bit : int, i_bit : int):
        self.number_qubits = number_qubits
        self.x_vector = x_vector
        self.y_vector = y_vector
        
        self.sign_bit = sign_bit
        self.i_bit = i_bit

    def generate_matrix(self) -> np.ndarray:
        size = 2**self.number_qubits
        matrix = np.zeros((size, size), dtype=complex)

        for j in range(size):
            matrix[j^self.x_vector, j] = self.get_phase()*phase_mod2product(j, self.y_vector)

        return matrix
    
    def get_phase(self) -> complex:
        return (1-2*self.sign_bit)*(1+(1j-1)*self.i_bit)

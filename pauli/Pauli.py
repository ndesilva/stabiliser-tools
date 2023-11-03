import numpy as np
from F2_helper.F2_helper import sign_mod2product

class Pauli:
    # Pauli is (-1)^sign bit * (-i) * (i_bit) * X^x_vector * Z^z_vector
    def __init__(self, number_qubits : int, x_vector : int, z_vector : int, sign_bit : int, i_bit : int):
        self.number_qubits = number_qubits
        self.x_vector = x_vector
        self.z_vector = z_vector
        
        self.sign_bit = sign_bit
        self.i_bit = i_bit

        self.phase = (1-2*self.sign_bit)*(1+(-1j-1)*self.i_bit)

    def generate_matrix(self) -> np.ndarray:
        size = 2**self.number_qubits
        matrix = np.zeros((size, size), dtype=complex)

        for j in range(size):
            matrix[j^self.x_vector, j] = self.phase*sign_mod2product(j, self.z_vector)

        return matrix
    
    def __eq__(self, other : object) -> bool:
        if not isinstance(other, Pauli):
            return False
        
        result = True
        result &= (self.sign_bit == other.sign_bit)
        result &= (self.i_bit == other.i_bit)
        result &= (self.number_qubits == other.number_qubits)
        result &= (self.x_vector == other.x_vector)
        result &= (self.z_vector == other.z_vector)

        return result
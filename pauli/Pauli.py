import numpy as np
import F2_helper.F2_helper as f2

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
            matrix[j^self.x_vector, j] = self.phase*f2.sign_mod2product(j, self.z_vector)

        return matrix
    
    def get_sign_eigenvalue(self, state_vector : np.ndarray) -> int | None:
        n = f2.fast_log2(len(state_vector))

        if n != self.number_qubits:
            raise ValueError('State vector and Paulis have different size')
        
        index = 0

        while not state_vector[index ^ self.x_vector]:
            if state_vector[index]:
                return None
            
            index += 1
            
        factor = self.phase * state_vector[index] * f2.sign_mod2product(self.z_vector, index) / state_vector[index ^ self.x_vector]

        match factor:
            case 1:
                bit = 0
            case -1:
                bit = 1
            case _:
                return None

        for remaining_index in range(index + 1, 1 << n):
            if self.phase * state_vector[remaining_index] * f2.sign_mod2product(self.z_vector, remaining_index) != factor * state_vector[remaining_index ^ self.x_vector]:
                return None
            
        return bit
    
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
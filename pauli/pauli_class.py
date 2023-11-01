import numpy as np

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
            matrix[j^self.x_vector, j] = self.get_phase()*Pauli.phase_mod2product(j, self.y_vector)

        return matrix
    
    def get_phase(self) -> complex:
        return (1-2*self.sign_bit)*(1+(1j-1)*self.i_bit)

    # returns (-1) ** the mod 2 inner product of the binary representations of x, y
    @staticmethod
    def phase_mod2product(x : int, y : int) -> int:
        product = x & y

        pairity = 0
        
        while product:
            pairity ^= 1
            product &= product -1

        return 1-2*pairity
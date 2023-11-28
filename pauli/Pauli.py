from __future__ import annotations

import numpy as np
import F2_helper.F2_helper as f2
import numba

class Pauli: # TODO add from_matrix method
    # Pauli is (-1)^sign bit * (-i) *=^ (i_bit) * X^x_vector * Z^z_vector
    def __init__(self, number_qubits : int, x_vector : int, z_vector : int, sign_bit : int, i_bit : int):
        self.number_qubits = number_qubits
        self.x_vector = x_vector
        self.z_vector = z_vector
        
        self.sign_bit = sign_bit
        self.i_bit = i_bit

        self.update_phase()

    def update_phase(self):
        self.phase = (1-2*self.sign_bit)*(1+(-1j-1)*self.i_bit)

    def generate_matrix(self) -> np.ndarray:
        size = 2**self.number_qubits
        matrix = np.zeros((size, size), dtype=complex)

        for j in range(size):
            matrix[j^self.x_vector, j] = self.phase*f2.sign_mod2product(j, self.z_vector)

        return matrix
    
    def multiply_by_pauli_on_right(self, other_pauli : Pauli):
        self.sign_bit ^= other_pauli.sign_bit ^ f2.mod2product(self.z_vector, other_pauli.x_vector) ^ (self.i_bit & other_pauli.i_bit)
        self.i_bit ^= other_pauli.i_bit

        self.x_vector ^= other_pauli.x_vector
        self.z_vector ^= other_pauli.z_vector

        self.update_phase()

    def is_hermitian(self) -> bool:
        return self.i_bit == f2.mod2product(self.x_vector, self.z_vector)
    
    def anticommutes_with(self, other_pauli : Pauli) -> int:
        return f2.mod2product(self.x_vector, other_pauli.z_vector) ^ f2.mod2product(self.z_vector, other_pauli.x_vector)
    
    def commutes_with(self, other_pauli : Pauli) -> int:
        return 1 ^ self.anticommutes_with(other_pauli)

    def get_sign_eigenvalue(self, state_vector : np.ndarray, assume_equation_holds: bool = False) -> int | None:
        n = f2.fast_log2(len(state_vector))

        if n != self.number_qubits:
            raise ValueError('State vector and Paulis have different size')
        
        index = get_first_nonzero_index(state_vector, self.x_vector)

        if index == None:
            return None
            
        factor = state_vector[index] * f2.sign_mod2product(self.z_vector, index) / state_vector[index ^ self.x_vector]

        match round_to_5dp(factor * self.phase):
            case 1:
                bit = 0
            case -1:
                bit = 1
            case _:
                return None
            
        if assume_equation_holds:
            return bit

        if remaining_entries_consistent(n, index, self.x_vector, self.z_vector, state_vector, factor):
            return bit
            
        return None
    
    def __eq__(self, other : object) -> bool:
        if not isinstance(other, Pauli):
            return False
        
        result = True
        result &= (self.sign_bit == other.sign_bit)
        result &= (self.i_bit == other.i_bit)
        result &= (self.number_qubits == other.number_qubits)
        result &= (self.x_vector == other.x_vector)
        result &= (self.z_vector == other.z_vector)
        result &= (self.phase == other.phase)

        return result
    
@numba.njit()
def remaining_entries_consistent(number_qubits : int, start_index : int, x_vector : int,  z_vector :int, state_vector : np.ndarray, factor : complex) -> bool:
    
    for remaining_index in range(start_index + 1, 1 << number_qubits):
        if state_vector[remaining_index] * f2.sign_mod2product(z_vector, remaining_index) != factor * state_vector[remaining_index ^ x_vector]:   
            return False
    
    return True

@numba.njit()
def get_first_nonzero_index(state_vector : np.ndarray, x_vector : int) -> int:
    index = 0

    while not state_vector[index ^ x_vector]:
        if state_vector[index]:
            return None
        
        index += 1
    
    return index

@numba.njit()
def round_to_5dp(value : complex) -> complex:
    return round(value.real, 5) + 1j*round(value.imag, 5)
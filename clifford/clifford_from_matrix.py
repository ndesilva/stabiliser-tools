import numpy as np
import stabiliser_state.stabiliser_from_state_vector as ssv
import F2_helper.F2_helper as f2
import math

# class Stabiliser_From_State_Vector:
#     def __init__(self, matrix : np.ndarray, allow_global_factor : bool = False, )

def is_clifford(matrix : np.ndarray, allow_global_factor : bool = False) -> bool: 
    n = f2.fast_log2(matrix.shape[0])

    if not columns_consistent(matrix, n, allow_global_factor):
        return False
    
    matrix_times_hadamard = multiply_by_hadamard_product(matrix, n)

    if not columns_consistent(matrix_times_hadamard, n, allow_global_factor):
        return False
    
    return True

def columns_consistent(matrix : np.ndarray, number_qubits : int, allow_global_factor : bool) -> bool:
    first_col_state = ssv.Stabiliser_From_State_Vector(matrix[:, 0], allow_global_factor)

    if not first_col_state.is_stab_state:
        return False
    
    first_col_stab_group = first_col_state.get_stab_state().get_check_matrix().paulis

    pauli_patterns = [0 for _ in range(number_qubits)]

    for pauli_index in range(number_qubits):
        pauli = first_col_stab_group[pauli_index]

        for j in range(number_qubits):
            phase_bit = pauli.get_sign_eigenvalue(matrix[:, 1<<j])
            
            if phase_bit == None:
                return False
            
            pauli_patterns[pauli_index] |= phase_bit * (1<<j)


    # if not is_full_rank(pauli_patterns, number_qubits):
    #     return False

    for pauli_index in range(number_qubits):
        pauli = first_col_stab_group[pauli_index]
        pauli_pattern = pauli_patterns[pauli_index]

        for col_index in range(1, 1 << number_qubits): # double checking the weight 1 hamming strings again
            phase_bit = pauli.get_sign_eigenvalue(matrix[:, col_index])
            
            if phase_bit != f2.mod2product(pauli_pattern, col_index):
                return False
    
    return True

def is_full_rank(vectors : list[int], number_vectors : int) -> bool:
    # copy the vectors list to not alter it
    vectors_copy = vectors.copy()

    # row reduce the vectors to row echelon form, stopping if we ever get an all zero vector
    for i in range(number_vectors):

        pivot_index = f2.fast_log2(vectors_copy[i])

        if pivot_index == -1:
            return False
        
        for j in range(i+1, number_vectors): # only need row echelon form (not reduced) to check LI
            vectors_copy[j] ^= (i!=j) * f2.get_bit_at(vectors_copy[j], pivot_index) * vectors_copy[i]

    return True

def multiply_by_hadamard_product(matrix : np.ndarray, number_qubits : int) -> np.ndarray:
    for j in range(number_qubits):
        unscaled_multiply_by_hadmard_at(matrix, j, number_qubits)

    return matrix/math.sqrt((1 << number_qubits))

def unscaled_multiply_by_hadmard_at(matrix : np.ndarray, hadamard_index : int, number_qubits : int) -> None:
    for tail in range(1 << hadamard_index):
        for head in range(1 << (number_qubits - hadamard_index - 1)):
            shifted_head = head << (hadamard_index + 1)
            
            first_col_index = shifted_head | tail
            second_col_index = shifted_head | (1 << hadamard_index) | tail

            matrix[:, first_col_index] += matrix[:, second_col_index]
            matrix[:, second_col_index] = matrix[:, first_col_index] - 2*matrix[:, second_col_index] # since first column is set first, this is now the difference
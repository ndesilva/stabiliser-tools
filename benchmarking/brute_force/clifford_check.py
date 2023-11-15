import pauli.Pauli as p
import pauli.pauli_check as pc
import numpy as np
import F2_helper.F2_helper as f2

def is_clifford(matrix : np.ndarray) -> bool:
    number_qubits = f2.fast_log2(matrix.shape[0])

    for i in range(number_qubits):
        z_i = p.Pauli(number_qubits, 0, 1<<i,0 , 0).generate_matrix()
        
        if not pc.is_pauli(matrix @ z_i @ matrix.conj().T, allow_global_factor = False):
            return False
        
        x_i = p.Pauli(number_qubits, 1<<i, 0, 0, 0).generate_matrix()

        if not pc.is_pauli(matrix @ x_i @ matrix.conj().T, allow_global_factor = False):
            return False
        
    return True
from __future__ import annotations

import numpy as np
import pauli.Pauli as p
import clifford.clifford_from_matrix as cm

class Clifford:

    @staticmethod
    def from_matrix(matrix: np.ndarray, assume_clifford : bool = False) -> Clifford: # TODO test
        clifford = cm.Clifford_From_Matrix(matrix, allow_global_factor = True, assume_clifford = assume_clifford)

        if not clifford.is_clifford:
            raise ValueError('Matrix does not correspond to a Clifford')
        
        return clifford.get_clifford()

    def __init__(self, z_conjugates : list[p.Pauli], x_cojugates : list[p.Pauli]):
        self.number_qubits = len(z_conjugates)
        
        self.z_conjugates = z_conjugates
        self.x_conjugates = x_cojugates
import numpy as np
import pauli.Pauli as p

class Clifford:

    def __init__(self, z_conjugates : list[p.Pauli], x_cojugates : list[p.Pauli]):
        self.number_qubits = len(z_conjugates)
        
        self.z_conjugates = z_conjugates
        self.x_conjugates = x_cojugates
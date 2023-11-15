import numpy as np

class Benchmarking_Data:
    def __init__(self, function_description : str, generator_description : str, number_qubits : list[int], times : np.ndarray):
        self.function_description = function_description
        self.generator_description = generator_description
        self.number_qubits = number_qubits
        self.times = times
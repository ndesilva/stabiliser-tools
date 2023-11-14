import matplotlib.pyplot as plt
import numpy as np
import pickle
from benchmarking.data.Benchmarking_Data import Benchmarking_Data

data_filename = 'benchmarking/data/test_function_with_test_generator.npy'

with open(data_filename, 'rb') as fl:
    data : Benchmarking_Data = pickle.load(fl)

print(data.function_description)
print(data.generator_description)
print(data.number_qubits)
print(data.times)

# _, ax = plt.subplots()
# ax.plot(qubit_numbers, pauli_times, label = 'random Pauli')
# ax.plot(qubit_numbers, random_matrix_times, label = 'random matrix')
# ax.plot(qubit_numbers, almost_pauli_times, label = 'random almost-Pauli')

# ax.set_xlabel(f'n')
# ax.set_ylabel('execution time (s)')
# ax.set_title(f'Benchmarking is_pauli() as the number of qubits increases')
# ax.legend()

# plt.show()
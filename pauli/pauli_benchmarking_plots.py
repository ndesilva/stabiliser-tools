import matplotlib.pyplot as plt
import numpy as np

data_filename = 'pauli_benchmarking_data.py'

with open(data_filename, 'rb') as fl:
    qubit_numbers = np.load(fl)
    
    pauli_times = np.load(fl)
    random_matrix_times = np.load(fl)
    almost_pauli_times = np.load(fl)

_, ax = plt.subplots()
ax.plot(qubit_numbers, pauli_times, label = 'random Pauli')
ax.plot(qubit_numbers, random_matrix_times, label = 'random matrix')
ax.plot(qubit_numbers, almost_pauli_times, label = 'random almost-Pauli')

ax.set_xlabel(f'n')
ax.set_ylabel('execution time (s)')
ax.set_title(f'Benchmarking is_pauli() as the number of qubits increases')
ax.legend()

plt.show()
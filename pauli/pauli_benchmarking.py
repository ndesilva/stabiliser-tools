import random
import pauli_check
import time
import numpy as np
import matplotlib.pyplot as plt

data_filename = 'pauli_benchmarking_data.py'

reps = int(1e3)

max_qubits = 10
qubit_numbers = range(1, max_qubits + 1)

pauli_entries = [1, -1, 1j, -1j]

def time_with_generator(matrix_geneneration) -> np.ndarray:
    times = np.zeros(max_qubits)
    
    for n in qubit_numbers:
        print(n)
        size = 1 << n
        timer = 0

        for _ in range(reps):
            matrix = matrix_geneneration(size, n)
            
            st = time.time()
            pauli_check.is_pauli(matrix)
            et = time.time()

            timer += et-st

        times[n-1] = timer/reps
    
    return times

'''Given a randomly generated Pauli, check how long it takes to accept'''
def random_pauli(size : int, n : int) -> np.ndarray:
    s = random.choice(pauli_entries)
    p = random.randrange(size)
    q = random.randrange(size)

    return pauli_check.generate_pauli(n, s, p, q)

'''Given a radomly generated matrix, check how long it takes to reject'''
def random_matrix(size : int, n : int) -> np.ndarray:
    return np.random.rand(size, size)

'''Given a random Pauli, break it in one entry, check how long it takes to reject'''
def random_almost_pauli(size : int, n : int):
    s = random.choice(pauli_entries)
    p = random.randrange(size)
    q = random.randrange(size)

    matrix = pauli_check.generate_pauli(n, s, p, q)

    i = random.randrange(size)
    j = random.randrange(size)

    matrix[i,j] = 0

    return matrix

pauli_times = time_with_generator(random_pauli)
random_matrix_times = time_with_generator(random_matrix)
almost_pauli_times = time_with_generator(random_almost_pauli)

print('Done!')

with open(data_filename, 'wb') as fl:
    np.save(fl, qubit_numbers)

    np.save(fl, pauli_times)
    np.save(fl, random_matrix_times)
    np.save(fl, almost_pauli_times)
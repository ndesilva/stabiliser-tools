import numpy as np
import random
from pauli.Pauli import Pauli

def random_matrix(n : int) -> np.ndarray:
    return np.random.rand(1<<n , 1<<n)

def random_pauli_matrix(n : int) -> np.ndarray:
    s = random.randrange(2)
    t = random.randrange(2)
    p = random.randrange(1<<n)
    q = random.randrange(1<<n)

    return Pauli(n, p, q, s, t).generate_matrix()

def random_almost_pauli_matrix(n : int):
    s = random.randrange(2)
    t = random.randrange(2)
    p = random.randrange(1<<n)
    q = random.randrange(1<<n)

    matrix = Pauli(n, p, q, s, t).generate_matrix()

    i = random.randrange(1 << n)
    j = random.randrange(1 << n)

    if matrix[i,j]:
        if random.randrange(1):
            matrix[i,j] *= 1j
        else:
            matrix[i,j] = 0

    return matrix
import numpy as np
import random
from pauli.Pauli import Pauli
import benchmarking.generator_dependencies.randstab as rs
import scipy.stats as sts

def random_unitary(n : int) -> np.ndarray:
    return sts.unitary_group.rvs(1<<n)

def random_pauli_matrix(n : int) -> np.ndarray:
    s = random.randrange(2)
    t = random.randrange(2)
    p = random.randrange(1<<n)
    q = random.randrange(1<<n)

    return Pauli(n, p, q, s, t).generate_matrix()

def random_almost_pauli_matrix(n : int) -> np.ndarray:
    matrix = random_pauli_matrix(n)
    modify_random_matrix_entry(n, matrix)

    return matrix


def random_state(n : int) -> np.ndarray:
    unnormalised = np.random.rand(1<<n) + 1j*np.random.rand(1<<n)
    return unnormalised / np.sqrt(unnormalised @ unnormalised.conj())

def random_stab_state(n : int) -> np.ndarray:
    return rs.random_stabilizer_state(n)

def random_almost_stab_state(n : int) -> np.ndarray:
    stab_state = random_stab_state(n)

    i = random.randrange( 1 << n )

    if stab_state[i]:
        if random.randrange(1):
            stab_state[i] *= 1j
        else:
            stab_state[i] = 0
    else:
        stab_state[i] = 1

    return stab_state

def random_clifford(n : int) -> np.ndarray:
    pass

def random_almost_clifford(n : int) -> np.ndarray:
    matrix = random_clifford(n)
    modify_random_matrix_entry(n, matrix)

    return matrix

def modify_random_matrix_entry(n : int, matrix : np.ndarray) -> None:
    i = random.randrange(1 << n)
    j = random.randrange(1 << n)

    if matrix[i,j]:
        if random.randrange(1):
            matrix[i,j] *= 1j
        else:
            matrix[i,j] = 0
    else:
        matrix[i,j] = 1
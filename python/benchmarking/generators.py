import numpy as np
import random
import sys
import benchmarking.generator_dependencies.randstab as rs
import scipy.stats as sts
import qiskit.quantum_info as qi

PATH_TO_LIBRARY = './build/ninja-multi-vcpkg/cpp/src/Release'
sys.path.append(PATH_TO_LIBRARY)
from fast import Pauli, Stabiliser_State

from typing import Tuple


def random_unitary(n: int) -> np.ndarray:
    return sts.unitary_group.rvs(1 << n)


def random_pauli_matrix(n: int) -> np.ndarray:
    s = random.randrange(2)
    t = random.randrange(2)
    p = random.randrange(1 << n)
    q = random.randrange(1 << n)

    return Pauli(n, p, q, s, t).generate_matrix()


def random_almost_pauli_matrix(n: int) -> np.ndarray:
    matrix = random_pauli_matrix(n)
    modify_random_matrix_entry(n, matrix)

    return matrix


def random_vector(n: int) -> np.ndarray:
    unnormalised = np.random.rand(1 << n) + 1j*np.random.rand(1 << n)
    return unnormalised / np.sqrt(unnormalised @ unnormalised.conj())


def random_stab_state(n: int) -> np.ndarray:
    return rs.random_stabilizer_state(n)


def random_stab_state_with_assump(n: int) -> Tuple[np.ndarray, bool]:
    return rs.random_stabilizer_state(n), True


def random_almost_stab_state(n: int) -> np.ndarray:
    stab_state = rs.random_stabilizer_state(n)

    i = random.randrange(1 << n)

    if stab_state[i]:
        if random.randrange(2):
            stab_state[i] *= 1j
        else:
            stab_state[i] = 0
    else:
        stab_state[i] = 1

    return stab_state


def worst_case_stab_state(n: int) -> np.ndarray:
    stab = Stabiliser_State(n)
    stab.basis_vectors = [1 << k for k in range(n)]
    stab.dim = 1 << n
    stab.imaginary_part = 1 << (n+1) - 1
    stab.real_linear_part = 1 << (n+1) - 1
    stab.quadratic_form = {}
    for i in range(n):
        for j in range(i+1, n):
            stab.quadratic_form[(1 << i) ** (1 << j)] = 1
    return np.array(stab.get_state_vector)


def worst_case_almost_stab_state(n: int) -> np.ndarray:
    stab_vector = worst_case_stab_state(n)
    i = random.randrange(1 << n)

    if stab_vector[i]:
        if random.randrange(2):
            stab_vector[i] *= 1j
        else:
            stab_vector[i] = 0
    else:
        stab_vector[i] = 1

    return stab_vector


def best_case_stab_state(n: int) -> np.ndarray:
    e_1 = np.zeros(1 << n)
    e_1[0] = 1
    return e_1


def best_case_stab_state_with_assump(n: int) -> Tuple[np.ndarray, bool]:
    return best_case_stab_state(n), True


def random_clifford(n: int) -> np.ndarray:
    return qi.random_clifford(n).to_matrix()


def random_almost_clifford(n: int) -> np.ndarray:
    matrix = random_clifford(n)

    if random.randrange(2):
        modify_random_matrix_entry(n, matrix)
    else:
        modify_random_column(n, matrix)

    return matrix


def modify_random_matrix_entry(n: int, matrix: np.ndarray) -> None:
    i = random.randrange(1 << n)
    j = random.randrange(1 << n)

    if matrix[i, j]:
        if random.randrange(2):
            matrix[i, j] *= 1j
        else:
            matrix[i, j] = 0
    else:
        matrix[i, j] = 1


def modify_random_column(n: int, matrix: np.ndarray) -> None:
    j = random.randrange(1 << n)
    matrix[:, j] *= 1j

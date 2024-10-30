import random
import stim
import sys
import numpy as np
import generator_dependencies.randstab as rs
import scipy.stats as sts
import qiskit.quantum_info as qi
from math import sqrt
from typing import Tuple

PATH_TO_LIBRARY = './build/ninja-multi-vcpkg/cpp/src/Release'
sys.path.append(PATH_TO_LIBRARY)
import fast as fst # type: ignore

def random_unitary(n: int) -> np.ndarray:
    return sts.unitary_group.rvs(1 << n)


def random_pauli_matrix(n: int) -> np.ndarray:
    s = random.randrange(2)
    t = random.randrange(2)
    p = random.randrange(1 << n)
    q = random.randrange(1 << n)

    return fst.Pauli(n, p, q, s, t).generate_matrix()


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


def random_full_support_stab_state(n: int) -> np.ndarray:
    stab = fst.Stabiliser_State(n)
    stab.basis_vectors = [1 << k for k in range(n)]
    stab.dim = n
    stab.imaginary_part = random.randrange(1 << n)
    stab.real_linear_part = random.randrange(1 << n)
    quadratic_form = {}
    quadratic_form[0] = 0
    for i in range(n):
        for j in range(i+1, n):
            quadratic_form[(1 << i) ^ (1 << j)] = random.randrange(2)
    stab.quadratic_form = quadratic_form # for now, the quadratic_form is not opaque (see https://pybind11.readthedocs.io/en/stable/advanced/cast/stl.html), so this is a workaround
    return np.array(stab.get_state_vector())


def random_full_support_almost_stab_state(n: int) -> np.ndarray:
    stab_vector = random_full_support_stab_state(n)
    i = random.randrange(1 << n)

    if stab_vector[i]:
        if random.randrange(2):
            stab_vector[i] *= 1j
        else:
            stab_vector[i] = 0
    else:
        stab_vector[i] = 1

    return stab_vector


def computational_zero(n: int) -> np.ndarray:
    e_1 = np.zeros(1 << n)
    e_1[0] = 1
    return e_1


def random_clifford(n: int) -> np.ndarray:
    matrix = qi.random_clifford(n).to_matrix()
    first_col_norm = matrix[0] @ matrix[0].conjugate()
    return matrix / sqrt(first_col_norm)


def random_clifford_with_assumption(n : int) -> Tuple[np.ndarray, bool]:
    return random_clifford(n), True

def get_identity_matrix(n : int) -> np.ndarray:
    return np.eye( 1<< n )

def get_Hadamard_matrix(n : int) -> np.ndarray:
    N = 1 << n
    factor = 1/sqrt(N)
    matrix = [ [factor * (1 - 2 * ((i & j).bit_count() & 1)) for i in range(N)] for j in range(N)]

    return np.array(matrix)

def get_anti_identiy_matrix(n : int) -> np.ndarray:
    return np.eye(1 << n)[::-1]

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


def rand_s_v(n: int) -> np.ndarray:
    func = random.sample(
        [random_stab_state, random_almost_stab_state], 1)[0]
    return func(n)


def rand_s_v_to_succinct(n: int) -> np.ndarray | Tuple[np.ndarray, bool]:
    # func = random.sample(
    #     [random_stab_state_with_assump, computational_zero, random_full_support_stab_state], 1)[0]
    func = random_stab_state_with_assump
    return func(n)


def rand_succinct(n: int) -> Tuple[fst.Stabiliser_State, stim.Tableau]:
    var = rand_s_v_to_succinct(n)
    if type(var) is tuple:
        statevector = var[0]
    else:
        statevector = var
    our_succinct = fst.stabiliser_state_from_statevector(statevector, assume_valid=True)
    stim_succinct = stim.Tableau.from_state_vector(statevector, endian='big')
    return our_succinct, stim_succinct


def rand_our_succinct(n: int) -> fst.Stabiliser_State:
    var = rand_s_v_to_succinct(n)
    if type(var) is tuple:
        statevector = var[0]
    else:
        statevector = var
    our_succinct = fst.stabiliser_state_from_statevector(statevector, assume_valid=True)
    return our_succinct


def rand_check_matrix(n: int) -> Tuple[fst.Check_Matrix, stim.Tableau]:
    statevector = random_stab_state_with_assump(n)[0]
    our_check_matrix = fst.Check_Matrix(fst.stabiliser_state_from_statevector(statevector, assume_valid=True))
    stim_check_matrix = stim.Tableau.from_state_vector(statevector, endian='big')
    return our_check_matrix, stim_check_matrix


def rand_our_check_matrix(n: int) -> fst.Check_Matrix:
    statevector = random_stab_state(n)
    return fst.Check_Matrix(fst.stabiliser_state_from_statevector(statevector, assume_valid=True))


def rand_clifford_test(n: int) -> np.ndarray:
    func = random.sample(
        [random_clifford, random_almost_clifford], 1)[0]
    return func(n)


def rand_clifford_matrix(n: int) -> np.ndarray | Tuple[np.ndarray, bool]:
    # func = random.sample(
    #     [random_clifford_with_assumption, get_identity_matrix, 
    #      get_Hadamard_matrix, get_anti_identiy_matrix], 1)[0]
    func = random_clifford_with_assumption
    return func(n)


def rand_clifford_succinct(n: int) -> Tuple[fst.Clifford, qi.Clifford, stim.Tableau]:
    qiskit_clifford = qi.random_clifford(n)
    mat = qiskit_clifford.to_matrix()
    our_clifford = fst.clifford_from_matrix(mat)
    stim_clifford = stim.Tableau.from_unitary_matrix(mat, endian='big')
    return our_clifford, qiskit_clifford, stim_clifford
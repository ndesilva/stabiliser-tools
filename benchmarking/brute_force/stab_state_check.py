import math
import itertools as it
import numpy as np
from numpy.linalg import norm

eps = 1e-5
rnd_dec = 5


def state_eq(state1, state2):
    diff = norm(state1 - state2)
    return diff < eps


def paulistring_to_matrix(paulistring):
    X = np.array([[0, 1],
                  [1, 0]])
    Y = np.array([[0, -1j],
                  [1j, 0]])
    Z = np.array([[1, 0],
                  [0, -1]])
    I = np.eye(2)
    st_to_pauli_matrices = {
        'I': I,
        'Z': Z,
        'X': X,
        'Y': Y
    }

    Pm = []
    for x in paulistring:
        if len(Pm) == 0:
            Pm = st_to_pauli_matrices[x]
        else:
            Pm = np.kron(Pm, st_to_pauli_matrices[x])

    return Pm


def is_stab_from_paulis(state: np.ndarray, n: int):
    state = np.round(state, decimals=rnd_dec)

    # find all stabilisers of the state:

    pauli_labels = ['I', 'X', 'Y', 'Z']

    pauli_combinations = tuple(it.product(pauli_labels, repeat=n))

    num_stabilisers = 0

    for labels in pauli_combinations:
        p_st = ''.join(labels)
        Pm = paulistring_to_matrix(p_st)

        resp = np.dot(Pm, state)
        resm = np.dot(-1 * Pm, state)

        if state_eq(resp, state) or state_eq(resm, state):
            num_stabilisers += 1

        if num_stabilisers >= (1 << n):
            return True

    return False


if __name__ == '__main__':
    # Is stab state
    state_vector = np.array([1] * 8)
    print(is_stab_from_paulis(state_vector, 3))

    # Not stab state
    state_vector = np.array(
        [1, 1, 1, 1, 1, 1, 1, (1/math.sqrt(2))*(1+1j)])/math.sqrt(8)
    print(is_stab_from_paulis(state_vector, 3))

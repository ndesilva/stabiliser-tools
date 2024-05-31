import math
import itertools as it
import numpy as np
import F2_helper.F2_helper as f2
from numpy.linalg import norm

eps = 1e-5
rnd_dec = 5


def state_eq(state1, state2):
    diff = norm(state1 - state2)
    return diff < eps


def state_eq_up_to_phase(state1, state2):
    # print(f'Comparing {state1 = }, {state2 = }')
    if (state_eq(state1, np.zeros_like(state1))
            and not state_eq(state2, np.zeros_like(state2))) or \
       (state_eq(state2, np.zeros_like(state2))
            and not state_eq(state1, np.zeros_like(state1))):
        return False, np.nan
    elif state_eq(state1, state2):
        return True, 1
    else:
        state2_rndd = np.round(state2, decimals=rnd_dec)
        i = np.nonzero(state2_rndd)[0][0]
        state2_renorm = (state1[i] / state2[i]) * state2
        if state_eq(state1, state2_renorm):
            return True, state2[i] / state1[i]
        else:
            return False, np.nan


def get_matrix(check_vector, phase):
    """
    Converts a Pauli element's check vector into its matrix form.

    Parameters
    ----------
    check_vector : numpy.ndarray
        The check vector.
    phase : complex
        The phase of the element: should be :math:`\pm 1`.

    Returns
    -------
    matrix : numpy.ndarray
        The matrix representation of the Pauli element.
    operator_str : string
        The Pauli element in a readable format.

    """

    X = np.array([[0, 1],
                  [1, 0]])
    Y = np.array([[0, -1j],
                  [1j, 0]])
    Z = np.array([[1, 0],
                  [0, -1]])
    I = np.eye(2)

    pauli_matrices = {
        '00': I,
        '01': Z,
        '10': X,
        '11': Y
    }
    pauli_matrix_strs = {
        '00': 'I',
        '01': 'Z',
        '10': 'X',
        '11': 'Y'
    }

    n = int(len(check_vector) / 2)
    matrix = 1
    operator_str = '' if phase == 1 else '-'
    for i in range(n):
        # Figure out what the operator on the ith qubit is
        ith_operator = pauli_matrices[str(check_vector[i]) +
                                      str(check_vector[n+i])]
        matrix = np.kron(matrix, ith_operator)
        operator_str \
            += pauli_matrix_strs[str(check_vector[i]) +
                                 str(check_vector[n+i])]

    return (phase * matrix, operator_str)


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


def is_stab_from_paulis(state: np.ndarray,):
    n = f2.fast_log2(state.shape[0])

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


def rref_binary(matrix):
    num_rows, num_cols = matrix.shape
    matrix = np.copy(matrix)

    row_to_comp = 0
    for j in range(num_cols):
        col = matrix[row_to_comp:, j]
        if np.count_nonzero(col) > 0:
            i = np.nonzero(col)[0][0] + row_to_comp
            matrix[(row_to_comp, i), :] = matrix[(i, row_to_comp), :]

            for ii in range(num_rows):
                if ii != row_to_comp and matrix[ii, j] != 0:
                    matrix[ii, :] ^= matrix[row_to_comp, :]

            row_to_comp += 1
            if row_to_comp == num_rows:
                break

    return matrix


def pauli_to_symplectic(pauli_string):
    x = []
    z = []
    for pauli in pauli_string:
        if pauli == 'I':
            x.append(0)
            z.append(0)
        elif pauli == 'X':
            x.append(1)
            z.append(0)
        elif pauli == 'Y':
            x.append(1)
            z.append(1)
        elif pauli == 'Z':
            x.append(0)
            z.append(1)
    return x + z


def stab_to_xmatr(state: np.ndarray):
    """
    Returns
    -------
    ind_matrix

    ind_signs
    """

    n = f2.fast_log2(state.shape[0])

    state = np.round(state, decimals=rnd_dec)

    # find all stabilizers of the state:

    pauli_labels = ['I', 'X', 'Y', 'Z']

    pauli_combinations = tuple(it.product(pauli_labels, repeat=n))

    stabset = []  # list of all stabilizers
    for labels in pauli_combinations:
        p_st = ''.join(labels)
        Pm = paulistring_to_matrix(p_st)

        resp = np.dot(Pm, state)
        resm = np.dot(-1 * Pm, state)

        if state_eq(resp, state):
            stabset.append(p_st)
        elif state_eq(resm, state):
            stabset.append('-'+p_st)

        if len(stabset) == (1 << n):
            break

    # convert stabilizers into binary symplectic format
    signs = []  # 0 if +, 1 if -
    symp_matr = []
    for stab in stabset:
        symp_str = pauli_to_symplectic(stab)
        symp_int = [int(bit) for bit in symp_str]
        symp_matr.append(symp_int)
        if stab[0] == '-':
            signs.append(1)
        else:
            signs.append(0)

    symp_matr = np.array(symp_matr)
    # print(symp_matr)

    # Convert check matrix to rref, and find new phases
    ind_matrix = rref_binary(symp_matr)[:n, :]

    signs = []
    for row in ind_matrix:
        pauli, _ = get_matrix(row, 1)
        is_it, phase = state_eq_up_to_phase(state, pauli @ state)
        signs.append(0 if phase == 1 else 1)
    ind_signs = np.array(signs)

    return ind_matrix, ind_signs


if __name__ == '__main__':
    # Is stab state
    state_vector = np.array([1] * 8)
    print(is_stab_from_paulis(state_vector))

    # Not stab state
    state_vector = np.array(
        [1, 1, 1, 1, 1, 1, 1, (1/math.sqrt(2))*(1+1j)])/math.sqrt(8)
    print(is_stab_from_paulis(state_vector))

    print(stab_to_xmatr(state_vector))

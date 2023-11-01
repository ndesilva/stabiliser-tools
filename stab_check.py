#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code for checking if a given state is a stabilizer state
"""

import numba

import numpy as np

# TODO Currently constraint matrices have only been generated for k = 1,...,6.
# For higher numbers of qubits, we can either generate and save
# more constraint matrices or just generate them on the fly.
with open("ConstraintMatrices", "r") as datafile:
    CM = eval(datafile.read().lower())
ConstraintMatrices = [np.array(x, dtype=np.int8) for x in CM]

# TODO This code is faster than the code below, but it needs to be fixed.
# E.g. the following matrix doesn't get row reduced correctly:
# mat = np.array([[1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0],
#                 [0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0],
#                 [0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1],
#                 [0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1],
#                 [0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0],
#                 [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0]], dtype=np.int8)
# @numba.njit(parallel=False)
# def gf2elim(M):
#     m, n = M.shape
#     i = 0
#     j = 0
#     while i < m and j < n:
#         k = np.argmax(M[i:, j]) + i
#         temp = np.copy(M[k])
#         M[k] = M[i]
#         M[i] = temp
#         aijn = M[i, j:]
#         # make a copy otherwise M will be directly affected
#         col = np.copy(M[:, j])
#         col[i] = 0  # avoid xoring pivot row with itself
#         flip = np.outer(col, aijn)
#         M[:, j:] = M[:, j:] ^ flip
#         i += 1
#         j += 1
#     return M


def gf2elim(M):
    num_rows, num_cols = M.shape
    M = np.copy(M)

    row_to_comp = 0
    for j in range(num_cols):
        col = M[row_to_comp:, j]
        if np.count_nonzero(col) > 0:
            i = np.nonzero(col)[0][0] + row_to_comp
            M[(row_to_comp, i), :] = M[(i, row_to_comp), :]

            for ii in range(num_rows):
                if ii != row_to_comp and M[ii, j] != 0:
                    M[ii, :] ^= M[row_to_comp, :]

            row_to_comp += 1
            if row_to_comp == num_rows:
                break

    return M


def is_stab(state: np.ndarray) -> bool:
    # Check if the support is an affine space by checking if the shifted
    # support is a vector space
    exponents = dict((2**k, k) for k in range(12))

    nonzero_indices = np.nonzero(state)[0]
    k = exponents[len(nonzero_indices)]

    supp = np.array([[int(c) for c in format(i, '06b')]
                    for i in nonzero_indices], dtype=np.int8)
    shift = supp[0]
    supp = shift ^ supp
    supp = gf2elim(supp)
    supp = supp[~np.all(supp == 0, axis=1)]

    if supp.shape[0] != k:
        return False

    # Now check if the amplitudes can be described by a quadratic polynomial
    # (quadratic form + constant)
    st = (1 - state[nonzero_indices]) // 2
    st = st[0] ^ st

    mat = gf2elim(np.column_stack(
        (ConstraintMatrices[k-1], st)).astype(np.int8))

    return mat[~np.all(mat == 0, axis=1)].shape[0] == k*(k+1) // 2

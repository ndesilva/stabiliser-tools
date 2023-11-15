import math
import time
import random
import functools
import numpy as np
from stabiliser_state.Stabiliser_State import Stabiliser_State
from F2_helper.F2_helper import *
import benchmarking.generator_dependencies.randstab as rs
import benchmarking.generators as gs

pauli_entries = [1, -1, 1j, -1j]

def list_check(float):
    return float in pauli_entries

def match_check(float):
    match float:
        case 1:
            return True
        case -1:
            return True
        case 1j:
            return True
        case _:
            return False
        
def or_check(float):
    if float == 1 or float == -1 or float == 1j or float == -1j:
        return True
    else:
        return False
    
l = [32,16,4,2]

def f2_linear_function_bitwise(integer):
    b = [ integer & index == index for index in l]
    return int(functools.reduce(lambda x,y : x^y, b))

lvec = np.array([1,1,0,1,1,0])

def f2_linear_function_matrix(integer):
    vector = np.array([int(c) for c in format(integer, '06b')])
    return np.dot(lvec, vector) % 2

l_int = 54

def f2_linear_function_int(integer):
    prod = l_int & integer
    wt = 0

    while prod:
        prod &= prod - 1
        wt ^=1

    return wt

def imaginary_part_linear(integer):
    return 1+(1j-1)*integer

def imaginary_part_power(integer):
    return (1j)**integer

def log2(integer):
    return integer & - integer

def get_bit_at_using_bin(int, i):
    return bin(int)[2+i]

def math_sqrt(x):
    return math.sqrt(x)

def np_sqrt(x):
    return np.sqrt(x)

def rref_binary(xmatr_aug):
    """
    'rref' function specifically for augmented check matrices.

    """
    num_rows, num_cols = xmatr_aug.shape
    xmatr_aug = np.copy(xmatr_aug)

    row_to_comp = 0
    for j in range(num_cols):
        col = xmatr_aug[row_to_comp:, j]
        if np.count_nonzero(col) > 0:
            i = np.nonzero(col)[0][0] + row_to_comp
            temp = np.copy(xmatr_aug[i, :])
            xmatr_aug[i, :] = xmatr_aug[row_to_comp, :]
            xmatr_aug[row_to_comp, :] = temp

            for ii in range(num_rows):
                if ii != row_to_comp and xmatr_aug[ii, j] != 0:
                    xmatr_aug[ii, :] ^= xmatr_aug[row_to_comp, :]

            row_to_comp += 1
            if row_to_comp == num_rows:
                break

    return xmatr_aug

def row_reduce(state : Stabiliser_State):
    state.__row_reduce_basis()

def pairity_1(integer):
    return bin(integer).count('1') & 1

def pairity_2(integer):
    p = 0

    while integer:
        p ^= 1
        integer &= integer - 1
    
    return p

def np_round(value : complex) -> complex:
    return np.round(value, 5)

@numba.njit()
def py_round(value : complex) -> complex:
    return round(value.real, 5) + 1j*round(value.imag, 5)

functions_to_time = [np_round, py_round]
reps = int(1e6)

for function in functions_to_time:

    timer = 0

    for i in range(reps):
        r = random.random() + 1j*random.random()

        st = time.perf_counter()       
        
        function(r)
        
        et = time.perf_counter()
        timer += et - st

    print(timer/reps)
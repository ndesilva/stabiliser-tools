import math
import time
import random
import functools
import numpy as np
from F2_helper.F2_helper import num_bin_digits, get_bit_at

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

def log2(power2):
    return int(np.log2(power2))

def get_bit_at_using_bin(int, i):
    return bin(int)[2+i]

def math_sqrt(x):
    return math.sqrt(x)

def np_sqrt(x):
    return np.sqrt(x)
    
functions_to_time = [math_sqrt, np_sqrt]
reps = int(1e6)

for function in functions_to_time:

    timer = 0

    for i in range(reps):
        r = random.randrange(1<<10)
       
        st = time.time()       
        
        function(r)
        
        et = time.time()
        timer += et - st

    print(timer/reps)
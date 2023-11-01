import time
import random
import functools
import numpy as np

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
    
l = [32,16,8,4,2]

def f2_linear_function_bitwise(integer):
    b = [ integer & index == index for index in l]
    return int(functools.reduce(lambda x,y : x^y, b))

lvec = np.array([1,1,0,1,1,0])

def f2_linear_function_matrix(integer):
    vector = np.array([int(c) for c in format(integer, '06b')])
    return np.dot(lvec, vector) % 2

print(f2_linear_function_bitwise(6), f2_linear_function_bitwise(3), f2_linear_function_bitwise(50))
print(f2_linear_function_matrix(6), f2_linear_function_matrix(3), f2_linear_function_matrix(50))

    
functions_to_time = [f2_linear_function_bitwise, f2_linear_function_matrix]
reps = int(1e7)

for function in functions_to_time:

    timer = 0

    for i in range(reps):
        r = random.randrange(2**6)
       
        st = time.time()       
        
        function(r)
        
        et = time.time()
        timer += et - st

    print(timer/reps)
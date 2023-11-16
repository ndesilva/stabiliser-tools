import numpy as np
import time
import pickle
from benchmarking.Benchmarking_Data import Benchmarking_Data
import benchmarking.generators as gs
import cProfile

import clifford.Clifford as c
import benchmarking.brute_force.clifford_check as bf_cc
import clifford.clifford_from_matrix as cc

import stabiliser_state.stabiliser_from_state_vector as ssc

reps = int(5)

min_qubits = 1
max_qubits = 12
qubit_numbers = list(range(min_qubits, max_qubits + 1))

profile = cProfile.Profile()

def time_function_with_generator(function_to_time, generator) -> np.ndarray:
    times = np.zeros(max_qubits)
    
    for n in qubit_numbers:
        print(f'n is {n}')
        timer = 0

        for j in range(reps):
            print(j)
            input = generator(n)

            st = time.perf_counter()
            # profile.enable()
            
            function_to_time(input)
            
            # profile.disable()
            et = time.perf_counter()

            timer += et-st

        times[n-1] = timer/reps
    
    return times

def our_method_clifford(matrix : np.ndarray) -> bool:
    return cc.Clifford_From_Matrix(matrix, only_testing = True).is_clifford

def our_method_stabiliser_state(state : np.ndarray) -> bool:
    return ssc.Stabiliser_From_State_Vector(state).is_stab_state

functions_to_time = [bf_cc.is_clifford, our_method_clifford] 
function_strings = ['brute force', 'our method']

generation_types = [gs.random_unitary, gs.random_clifford, gs.random_almost_clifford]
generation_strings = ['random unitary', 'random clifford', 'random perturbed clifford']

pre_string = 'testing C1'

num_functions = len(functions_to_time)
num_generators = len(generation_types)

assert num_functions == len(function_strings)
assert num_generators == len(generation_strings) 

base_filestring = './benchmarking/data'

for function_index in range(num_functions):
    for generation_index in range(num_generators):
        function_string = function_strings[function_index]
        generation_string = generation_strings[generation_index]
        
        filename = f'{base_filestring}/{pre_string} {function_string} on {generation_string}.npy'

        print(f'timing {function_string} with {generation_string}')

        times = time_function_with_generator(functions_to_time[function_index], generation_types[generation_index])

        data = Benchmarking_Data(function_string, generation_string, qubit_numbers, times)

        with open(filename, 'wb') as fl:
            pickle.dump(data, fl)

# profile.dump_stats('./logs/numba_with_compound_ifs.log')
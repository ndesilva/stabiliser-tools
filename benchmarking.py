import numpy as np
import time
import pickle
from benchmarking.Benchmarking_Data import Benchmarking_Data
import benchmarking.generators as gs

import clifford.Clifford as c
import benchmarking.brute_force.clifford_check as bf_cc

reps = int(1)

max_qubits = 12
qubit_numbers = list(range(1, max_qubits + 1))

def time_function_with_generator(function_to_time, generator) -> np.ndarray:
    times = np.zeros(max_qubits)
    
    for n in qubit_numbers:
        print(n)
        timer = 0

        for _ in range(reps):
            input = generator(n)
            
            st = time.perf_counter()
            function_to_time(input)
            et = time.perf_counter()

            timer += et-st

        times[n-1] = timer/reps
    
    return times

def our_method_clifford(matrix : np.ndarray) -> c.Clifford:
    return c.Clifford.from_matrix(matrix, assume_clifford = True)

functions_to_time = [bf_cc.is_clifford, our_method_clifford] 
function_strings = ['brute force', 'our method']

generation_types = [gs.random_clifford]
generation_strings = ['random clifford']

pre_string = 'converting C1 to C2'

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
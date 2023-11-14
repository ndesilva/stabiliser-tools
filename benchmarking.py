import numpy as np
import time
import pickle
from benchmarking.data.Benchmarking_Data import Benchmarking_Data

reps = int(1e3)

max_qubits = 10
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

test_lambda = lambda x : x^2

functions_to_time = [test_lambda] 
function_strings = ['test_function']

test_generator = lambda x : 1

generation_types = [test_generator]
generation_strings = ['test_generator']

num_functions = len(functions_to_time)
num_generators = len(generation_types)

assert num_functions == len(function_strings)
assert num_generators == len(generation_strings) 

base_filestring = './benchmarking/data'

for function_index in range(num_functions):
    for generation_index in range(num_generators):
        function_string = function_strings[function_index]
        generation_string = generation_strings[generation_index]
        
        filename = f'{base_filestring}/{function_string}_with_{generation_string}.npy'

        times = time_function_with_generator(functions_to_time[function_index], generation_types[generation_index])

        data = Bechmarking_Data(function_string, generation_string, qubit_numbers, times)

        with open(filename, 'wb') as fl:
            pickle.dump(data, fl)
import numpy as np
import time
import pickle
from Benchmarking_Data import Benchmarking_Data
from benchmarking_config import configs
import cProfile

reps = int(1e3)

min_qubits = 3
max_qubits = 11
qubit_numbers = list(range(min_qubits, max_qubits + 1))

profile = cProfile.Profile()

def time_function_with_generator(function_to_time, generator) -> np.ndarray:
    times = np.zeros(max_qubits - min_qubits + 1)

    for n in qubit_numbers:
        print(f'n is {n}')
        timer = 0

        for _ in range(reps):
            input = generator(n)

            if type(input) is tuple:
                st = time.perf_counter()
                # profile.enable()
                function_to_time(*input)
                # profile.disable()
                et = time.perf_counter()
            else:
                st = time.perf_counter()
                # profile.enable()
                function_to_time(input)
                # profile.disable()
                et = time.perf_counter()

            timer += et-st

        times[n - min_qubits] = timer/reps

    return times


base_filestring = './python/benchmarking/data'


def append_benchmarking_data(pre_string='', title='', functions_to_time=[], 
                             function_strings=[], generation_types=[], generation_strings=[]):
    num_functions = len(functions_to_time)
    num_generators = len(generation_types)

    assert num_functions == len(function_strings)
    assert num_generators == len(generation_strings)

    for function_index in range(num_functions):
        for generation_index in range(num_generators):
            function_string = function_strings[function_index]
            generation_string = generation_strings[generation_index]

            filename = f'{base_filestring}/{pre_string} {function_string} on {generation_string}.npy'

            print(f'timing {function_string} with {generation_string}')

            times = time_function_with_generator(
                functions_to_time[function_index], generation_types[generation_index])

            data = Benchmarking_Data(
                function_string, generation_string, qubit_numbers, times, title)

            with open(filename, 'wb') as fl:
                pickle.dump(data, fl)


if __name__ == '__main__':
    for config in configs:
        append_benchmarking_data(**config)

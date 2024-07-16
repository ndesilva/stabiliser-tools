import numpy as np
import time
import pickle
from benchmarking.Benchmarking_Data import Benchmarking_Data
from benchmarking_config import configs
import cProfile

reps = int(5)

min_qubits = 5
max_qubits = 10
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


base_filestring = './benchmarking/data'


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
                function_string, generation_string, qubit_numbers, times)

            with open(filename, 'wb') as fl:
                pickle.dump(data, fl)


if __name__ == '__main__':
    for config in configs:
        append_benchmarking_data(**config)

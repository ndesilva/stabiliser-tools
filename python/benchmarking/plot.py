import matplotlib.pyplot as plt
import numpy as np
import pickle
from Benchmarking_Data import Benchmarking_Data
from benchmarking_config import configs
from enum import Enum

base_data_path = './python/benchmarking/data/'
base_output_path = './python/benchmarking/plots/'

def make_plots(pre_string='', title='', function_strings='', generation_strings='', **kwargs):
    data_to_plot = []

    for function_string in function_strings:
        for generation_string in generation_strings:
            file_name = f'{base_data_path}{pre_string} {function_string} on {generation_string}.npy'

            with open(file_name, 'rb') as fl:
                data : Benchmarking_Data = pickle.load(fl)
                data_to_plot.append(data)

    fig, ax = plt.subplots()

    for data in data_to_plot:
        ax.plot(data.number_qubits, data.times, label = f'{data.function_description} with {data.generator_description}')

    ax.set_xlabel(f'n')
    ax.set_ylabel('execution time (s)')
    ax.set_title(title)
    ax.legend()
    ax.set_xticks(data.number_qubits)

    fig.savefig(f'{base_output_path}{title}.pdf')

class Time_type(Enum):
    BEST = 0
    AVERAGE = 1
    WORST = 0

class Line():
    def __init__(self, pre_string : str, function_string : str, generation_string : str, label : str, time_type : Time_type):
        self.pre_string = pre_string
        self.function_string = function_string
        self.generation_string = generation_string
        self.label = label
        self.time_type = time_type

    def get_filename(self):
        return f'{base_data_path}{self.pre_string} {self.function_string} on {self.generation_string}.npy'
    
    def add_data(self, data : Benchmarking_Data):
        self.qubit_numbers = data.number_qubits

        match self.time_type:
            case Time_type.WORST:
                self.times = data.worst_times
            case Time_type.AVERAGE:
                self.times = data.times
            case Time_type.BEST:
                self.times = data.best_times

def make_alternative_plots():
    lines = [
        Line("converting S1 to efficient rep", "our method", "random stab state without assump", "Our algorithm, best", Time_type.BEST),
        Line("converting S1 to efficient rep", "our method", "random stab state without assump", "Our algorithm, average", Time_type.AVERAGE),
        Line("converting S1 to efficient rep", "our method", "random stab state without assump", "Our algorithm, worst", Time_type.WORST),
        Line("converting S1 to efficient rep", "our method", "best case stab state", "Our algorithm, alt best state", Time_type.BEST),
        # Line("converting S1 to efficient rep", "stim", "random stab state without assump", "Stim, average", Time_type.AVERAGE)
    ]

    fig, ax = plt.subplots()

    for line in lines:
        with open(line.get_filename(), 'rb') as fl:
            data : Benchmarking_Data = pickle.load(fl)
            line.add_data(data)

        ax.plot(line.qubit_numbers, line.times, label = line.label)
    
    ax.set_xlabel(f'n')
    ax.set_ylabel('execution time (s)')
    ax.set_title("Stabiliser state benchmarking")
    ax.legend()
    ax.set_xticks(lines[0].qubit_numbers)

    fig.savefig(f'{base_output_path}Stabiliser.pdf')
    print(lines[0].times)


if __name__ == '__main__':
    for config in configs:
        make_plots(**config)

    make_alternative_plots()

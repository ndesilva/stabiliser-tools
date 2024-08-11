import matplotlib.pyplot as plt
import numpy as np
import pickle
from Benchmarking_Data import Benchmarking_Data
from benchmarking_config import configs

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

class Line():
    def __init__(self, pre_string : str, function_string : str, generation_string : str, label : str):
        self.pre_string = pre_string
        self.function_string = function_string
        self.generation_string = generation_string
        self.label = label

    def get_filename(self):
        return f'{base_data_path}{self.pre_string} {self.function_string} on {self.generation_string}.npy'
    
    def add_data(self, data : Benchmarking_Data):
        self.qubit_numbers = data.number_qubits
        self.times = data.times

def make_alternative_plots():
    lines = [
        Line("converting S1 to efficient rep", "our method", "best case stab state", "Our algorithm, best case"),
        Line("converting S1 to efficient rep", "our method", "random stab state without assump", "Our algorithm, average"),
        Line("converting S1 to efficient rep", "our method", "worst case stab state without assump", "Our algorithm, worst case"),
        Line("converting S1 to efficient rep", "stim", "random stab state without assump", "Stim, average")
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


if __name__ == '__main__':
    for config in configs:
        make_plots(**config)

    make_alternative_plots()

import pickle
import matplotlib.pyplot as plt
from matplotlib import rc, use
from Benchmarking_Data import Benchmarking_Data
from math import floor, log10
from enum import Enum
from typing import List, Tuple
from benchmarking_config import configs

base_data_path = './python/benchmarking/data/'
base_output_path = './python/benchmarking/plots/'

use('pgf')
rc('font', family='serif')

class Numerical_Type(Enum):
    LINEAR = 0
    LOG = 1

class Plot_Type(Enum):
    PERCENTILE = 0
    AVERAGE = 1

class Line():
    def __init__(self, pre_string : str, function_string : str, generation_string : str, label : str, percentile : int = 50, numerical_type : Numerical_Type = Numerical_Type.LINEAR, plot_type : Plot_Type = Plot_Type.PERCENTILE):
        self.pre_string = pre_string
        self.function_string = function_string
        self.generation_string = generation_string
        
        self.label = label
        
        self.percentile = percentile
        self.numerical_type = numerical_type
        self.plot_type = plot_type

        with open(self.get_filename(), 'rb') as fl:
            data : Benchmarking_Data = pickle.load(fl)
            self.add_data(data)

    def get_filename(self):
        return f'{base_data_path}{self.pre_string} {self.function_string} on {self.generation_string}.npy'
    
    def add_data(self, data : Benchmarking_Data):
        self.qubit_numbers = data.number_qubits

        index = floor(self.percentile * (data.reps - 1) /100)
        
        raw_times = [data.times[qubit_index][index] for qubit_index in range(len(data.number_qubits))]

        match self.numerical_type:
            case Numerical_Type.LINEAR:
                self.times = raw_times
            case Numerical_Type.LOG:
                self.times = [log10(time) for time in raw_times]

    def plot(self, axis, low_percentile : int = 10, high_percentile : int = 90):
        match self.plot_type:
            case Plot_Type.PERCENTILE:
                self.plot_percentile(axis, low_percentile, high_percentile)
            case Plot_Type.AVERAGE:
                self.plot_average_line(axis)

    def plot_average_line(self, axis):
        match self.function_string:
            case 'our method':
                c = '#4AAFD5'
            case 'stim':
                c = '#91B187'
            case 'Qiskit':
                c = '#E7A339'
        return axis.plot(self.qubit_numbers, self.times, label = self.label, color=c)

    def plot_percentile(self, axis, low_percentile : int = 10, high_percentile : int = 90):
        low_line = Line(self.pre_string, self.function_string, self.generation_string, self.label, numerical_type = self.numerical_type, percentile = high_percentile)
        high_line = Line(self.pre_string, self.function_string, self.generation_string, self.label, numerical_type = self.numerical_type, percentile = low_percentile)

        line_list = self.plot_average_line(axis)
        last_color = line_list[-1].get_color()
        axis.fill_between(self.qubit_numbers, low_line.times, high_line.times, alpha = 0.2, color = last_color)

def make_plot(pre_string : str, line_specs : List[Tuple[str, str, str, Plot_Type]], 
              title : str, low_percentile : int = 33, high_percentile : int = 66):
    fig, ax = plt.subplots()
    
    for spec in line_specs:
        line = Line(pre_string, spec[0], spec[1], spec[2], numerical_type = Numerical_Type.LOG, plot_type = spec[3])
        line.plot(ax, low_percentile, high_percentile)
        ax.set_xticks(line.qubit_numbers)
    
    ax.set_xlabel(r'$n$')
    ax.set_ylabel(r'$\log_{10}$(execution time in s)')
    ax.set_title(title)
    if len(line_specs) > 1:
        ax.legend()

    print(f'Saving figure {base_output_path}{title}.pdf')
    fig.savefig(f'{base_output_path}{title}.pdf')

    print(f'Saving figure {base_output_path}{title}.pgf')
    fig.savefig(f'{base_output_path}{title}.pgf')

if __name__ == '__main__':
    recapitalisation = {'our method': 'Our method',
                        'stim': 'Stim',
                        'Qiskit': 'Qiskit'}

    for config in configs:
        make_plot(
            config['pre_string'],
            [(fs, config['generation_strings'][0], recapitalisation[fs], Plot_Type.AVERAGE) for fs in config['function_strings']],
            config['title']
        )

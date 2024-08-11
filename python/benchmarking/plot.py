import matplotlib.pyplot as plt
import numpy as np
import pickle
from Benchmarking_Data import Benchmarking_Data
from benchmarking_config import configs
from math import floor, log10
from enum import Enum

base_data_path = './python/benchmarking/data/'
base_output_path = './python/benchmarking/plots/'

class Plot_Type(Enum):
    LINEAR = 0
    LOG = 1

class Line():
    def __init__(self, pre_string : str, function_string : str, generation_string : str, label : str, percentile : int = 50, plot_type : Plot_Type = Plot_Type.LINEAR):
        self.pre_string = pre_string
        self.function_string = function_string
        self.generation_string = generation_string
        
        self.label = label
        
        self.percentile = percentile
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

        match self.plot_type:
            case Plot_Type.LINEAR:
                self.times = raw_times
            case Plot_Type.LOG:
                self.times = [log10(time) for time in raw_times]

    def plot_line(self, axis):
        axis.plot(self.qubit_numbers, self.times, label = self.label)

    def plot_percentile(self, axis, low_percentile : int = 10, high_percentile : int = 90):
        low_line = Line(self.pre_string, self.function_string, self.generation_string, self.label, plot_type = Plot_Type.LOG, percentile = high_percentile)
        high_line = Line(self.pre_string, self.function_string, self.generation_string, self.label, plot_type = Plot_Type.LOG, percentile = low_percentile)

        axis.plot(self.qubit_numbers, self.times, label = self.label)
        axis.fill_between(self.qubit_numbers, low_line.times, high_line.times, alpha = 0.2)

def make_plot(pre_string : str, line_specs : list[tuple[str, str, str]], title : str, low_percentile : int = 10, high_percentile : int = 90):
    fig, ax = plt.subplots()
    
    for spec in line_specs:
        line = Line(pre_string, spec[0], spec[1], spec[2], plot_type = Plot_Type.LOG)
        line.plot_percentile(ax, low_percentile, high_percentile)
        ax.set_xticks(line.qubit_numbers)
    
    ax.set_xlabel(f'n')
    ax.set_ylabel('log execution time (s)')
    ax.set_title(title)
    ax.legend()

    fig.savefig(f'{base_output_path}{title}.pdf')

if __name__ == '__main__':
    make_plot(
        "converting S1 to efficient rep",
        [
            ("our method", "random stab state without assump", "Our method, random stabiliser state (without assumption)"),
            ("our method", "random stab state with assump", "Our method, random stabiliser state (with assumption)"),
            ("stim", "random stab state with assump", "Stim, random stabiliser state")
        ],
        "(S1) to efficient representation, comparison with stim"
    )

    make_plot(
        "converting S1 to efficient rep",
        [
            ("our method", "random stab state without assump", "Random stabiliser state (without assumption)"),
            ("our method", "random stab state with assump", "Random stabiliser state (with assumption)"),
            ("our method", "computational zero", "Computational 0 state"),
            ("our method", "random full support stabiliser state", "Random full support stabiliser state")
        ],
        "(S1) to efficient representation, different inputs"
    )

    make_plot(
        "testing (S1)",
        [
            ("our method", "random stabiliser state", "Our method, random stabiliser state"),
            ("our method", "almost stab state", 'Our method, random "almost" stabiliser state'),
            ("stim", "random stabiliser state", "Stim, random stabiliser state"),
            ("stim", "almost stab state", 'Stim, random "almost" stabiliser state')
        ],
        "Testing (S1), comparison with stim"
    )

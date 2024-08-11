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

def make_plots(pre_string='', title='', function_strings='', generation_strings='', **kwargs):
    lines_to_plot : list[Line] = []

    for function_string in function_strings:
        for generation_string in generation_strings:
            line = Line(pre_string, function_string, generation_string, f'{function_string} with {generation_string}')

            with open(line.get_filename(), 'rb') as fl:
                data : Benchmarking_Data = pickle.load(fl)
                line.add_data(data)
            
            lines_to_plot.append(line)

    fig, ax = plt.subplots()

    for line in lines_to_plot:
        line.plot_line(ax)

    ax.set_xlabel(f'n')
    ax.set_ylabel('execution time (s)')
    ax.set_title(title)
    ax.legend()
    ax.set_xticks(lines_to_plot[0].qubit_numbers)

    fig.savefig(f'{base_output_path}{title}.pdf')
    plt.close()

def stim_comparison_plot(pre_string : str, generation_string : str, title : str,  low_percentile = 10, high_percentile = 90, plot_type = Plot_Type.LOG):
    lines = [
        Line(pre_string, "our method", generation_string, f"Our algorithm, {low_percentile}th percentile", low_percentile, plot_type = Plot_Type.LOG),
        Line(pre_string, "our method", generation_string, f"Our algorithm", plot_type = Plot_Type.LOG),
        Line(pre_string, "our method", generation_string, f"Our algorithm, {high_percentile}th percentile", high_percentile, plot_type = Plot_Type.LOG),
        Line(pre_string, "stim", generation_string, f"Stim, {low_percentile}th percentile", low_percentile, plot_type = Plot_Type.LOG),
        Line(pre_string, "stim", generation_string, f"Stim", plot_type = Plot_Type.LOG),
        Line(pre_string, "stim", generation_string, f"Stim, {high_percentile}th percentile", high_percentile, plot_type = Plot_Type.LOG),
    ]

    for line in lines:
        with open(line.get_filename(), 'rb') as fl:
            data : Benchmarking_Data = pickle.load(fl)
            line.add_data(data)

    fig, ax = plt.subplots()

    lines[1].plot_line(ax)
    lines[4].plot_line(ax)
    ax.fill_between(lines[0].qubit_numbers, lines[0].times, lines[2].times, alpha = 0.2)
    ax.fill_between(lines[0].qubit_numbers, lines[3].times, lines[5].times, alpha = 0.2)
    
    ax.set_xlabel(f'n')
    ax.set_ylabel('log execution time (s)')
    ax.set_title(title)
    ax.legend()
    ax.set_xticks(lines[0].qubit_numbers)

    fig.savefig(f'{base_output_path}{title}.pdf')
    plt.show()

if __name__ == '__main__':
    for config in configs:
        make_plots(**config)

    stim_comparison_plot("converting S1 to efficient rep", "random stab state without assump", "(S1) to efficient representation, comparison with stim")
    stim_comparison_plot("non-stab reject", "worst case non-stab state", "Rejecting almost (S1), comparison with stim")
    stim_comparison_plot("non-stab reject", "random stabiliser state", "Accepting (S1), comparison with stim")

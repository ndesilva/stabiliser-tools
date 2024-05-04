#include <iostream>
#include <vector>
#include <complex>
#include <chrono>
#include <fstream>
#include <numeric>

#include "./stim/circuit_vs_amplitudes.h"
#include "stabiliser_state/stabiliser_state_from_statevector.h"

using namespace std::chrono;

std::vector<std::complex<float>> get_last_basis_vector(int number_qubits) {
    size_t vector_length = 1 << number_qubits;
    
    std::vector<std::complex<float>> statevector (vector_length, 0);

    statevector.at(vector_length - 1) = 1;

    return statevector;
}

std::vector<std::complex<float>> get_uniform_superposition(int number_qubits) {
    size_t vector_length = 1 << number_qubits;
    float normalisation_factor = 1/sqrt(vector_length);
    
    std::vector<std::complex<float>> statevector (vector_length, normalisation_factor);

    return statevector;
}

std::vector<std::complex<float>> get_non_stabiliser_state(int number_qubits) {
    size_t vector_length = 1 << number_qubits;
    float normalisation_factor = 1/sqrt(vector_length);
    
    std::vector<std::complex<float>> statevector (vector_length, normalisation_factor);

    statevector.at(0) *= -1;

    return statevector;
}

int main() {
    const std::string FILE_NAME = "./benchmarking/stabiliser_bechmark_data.txt";
    const int REPETITIONS = 1000;

    const int MIN_N = 3;
    const int MAX_N = 10;

    const int STATEVECTOR_TYPES = 3;

    int number_ns = MAX_N - MIN_N + 1;

    std::vector<int> qubit_numbers (number_ns);
    std::iota(qubit_numbers.begin(), qubit_numbers.end(), MIN_N);

    std::vector<std::vector<double>> ratios (STATEVECTOR_TYPES);

    std::vector<double> basis_vector_ratios (number_ns);
    std::vector<double> uniform_superposition_ratios (number_ns);
    std::vector<double> non_stabiliser_ratios (number_ns);

    std::vector<std::vector<std::complex<float>>> statevectors;

    std::cout << "Beginning benchmarking" << std::endl;

    for (const auto number_qubits : qubit_numbers) {
        statevectors.clear();

        statevectors.push_back(get_last_basis_vector(number_qubits));
        statevectors.push_back(get_uniform_superposition(number_qubits));
        statevectors.push_back(get_non_stabiliser_state(number_qubits));

        bool result;
    
        for (int i = 0; i < STATEVECTOR_TYPES; i++) {
            std::vector<std::complex<float>> statevector = statevectors[i];
            
            long long stim_total_time = 0;
            for (int j = 0; j < REPETITIONS; j++) {
                auto stim_start_time = high_resolution_clock::now();

                result = stim::stabilizer_state_vector_to_circuit(statevector);

                auto stim_end_time = high_resolution_clock::now();
                stim_total_time += duration_cast<microseconds>(stim_end_time - stim_start_time).count();
            }

            double stim_average_time = double(stim_total_time)/REPETITIONS;

            long long our_total_time = 0;
            for (int j = 0; j < REPETITIONS; j++) {
                auto our_start_time = high_resolution_clock::now();

                fst::Stabiliser_From_Vector_Convertor convertor (statevector);
                result = convertor.is_stabiliser_state;

                auto our_end_time = high_resolution_clock::now();
                our_total_time += duration_cast<microseconds>(our_end_time - our_start_time).count();
            }

            double our_average_time = double(our_total_time)/REPETITIONS;

            double time_improvement = double(stim_average_time)/our_average_time;
            ratios.at(i).push_back(time_improvement);
        }
    }
    
    std::cout << "done getting data, saving to file" << std::endl;

    std::ofstream data_file (FILE_NAME);
    
    for (int i = 0; i < STATEVECTOR_TYPES; i++){
        data_file << "[";

        for (const auto ratio : ratios.at(i)) {
            data_file << ratio << ", ";
        }

        data_file << "]\n";
    }

    return 0;
}
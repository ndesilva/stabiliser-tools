#include <iostream>
#include <vector>
#include <complex>
#include <chrono>
#include <fstream>
#include <numeric>

#include "stim/circuit_vs_amplitudes.h"
#include "stabiliser_state_from_statevector.h"
#include "f2_helper.h"

using namespace std::chrono;

std::vector<std::complex<float>> get_last_basis_vector(const std::size_t number_qubits) {
    const std::size_t vector_length = fst::integral_pow_2( number_qubits );
    
    std::vector<std::complex<float>> statevector (vector_length, 0);

    statevector.back() = 1;

    return statevector;
}

std::vector<std::complex<float>> get_uniform_superposition( const std::size_t number_qubits) {
    const std::size_t vector_length = fst::integral_pow_2( number_qubits );
	const float normalisation_factor = static_cast<float>( 1.0f / std::sqrt( vector_length ) );
    
    return std::vector<std::complex<float>>(vector_length, normalisation_factor);
}

std::vector<std::complex<float>> get_non_stabiliser_state( const std::size_t number_qubits) {
    const std::size_t vector_length = fst::integral_pow_2( number_qubits );
    const float normalisation_factor = static_cast<float>( 1.0f / std::sqrt( vector_length ) );
    
    std::vector<std::complex<float>> statevector (vector_length, normalisation_factor);

	statevector.back() *= -1;

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

    std::vector<bool> statevector_is_stabiliser (STATEVECTOR_TYPES);
    std::vector<std::vector<double>> ratios (STATEVECTOR_TYPES);

    std::vector<double> basis_vector_ratios (number_ns);
    std::vector<double> uniform_superposition_ratios (number_ns);
    std::vector<double> non_stabiliser_ratios (number_ns);

    std::vector<std::vector<std::complex<float>>> statevectors;

    std::cout << "Beginning benchmarking" << std::endl;

    for (const auto number_qubits : qubit_numbers) {
        statevectors.clear();
        statevector_is_stabiliser.clear();

        statevectors.push_back(get_last_basis_vector(number_qubits));
        statevector_is_stabiliser.push_back(true);

        statevectors.push_back(get_uniform_superposition(number_qubits));
        statevector_is_stabiliser.push_back(true);

        statevectors.push_back(get_non_stabiliser_state(number_qubits));
        statevector_is_stabiliser.push_back(false);

        bool result;
    
        for (int i = 0; i < STATEVECTOR_TYPES; i++) {
            std::vector<std::complex<float>> statevector = statevectors[i];
            bool expected_result = statevector_is_stabiliser[i];
            
            long long stim_total_time = 0;
            for (int j = 0; j < REPETITIONS; j++) {
                auto stim_start_time = high_resolution_clock::now();

                result = stim::stabilizer_state_vector_to_circuit(statevector);

                auto stim_end_time = high_resolution_clock::now();
                stim_total_time += duration_cast<microseconds>(stim_end_time - stim_start_time).count();

                if (result != expected_result)
                {
                    std::cout << "STIM error, expected: " << expected_result << "got: " << result << "on index: " << i << std::endl;
                }
            }

            double stim_average_time = double(stim_total_time)/REPETITIONS;

            long long our_total_time = 0;
            for (int j = 0; j < REPETITIONS; j++) {
                auto our_start_time = high_resolution_clock::now();

                result = fst::is_stabiliser_state(statevector);

                auto our_end_time = high_resolution_clock::now();
                our_total_time += duration_cast<microseconds>(our_end_time - our_start_time).count();

                if (result != expected_result)
                {
                    std::cout << "OUR error, expected: " << expected_result << " got " << result << " on index " << i << " qubit number: " << number_qubits << std::endl;
                }
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
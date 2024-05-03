#include <iostream>
#include <vector>
#include <complex>
#include <chrono>
#include <fstream>

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
    std::cout << "###\nBenchmarking Stim\n###\n";

    std::vector<std::vector<std::complex<float>>> statevectors;

    statevectors.push_back(get_last_basis_vector(10));
    statevectors.push_back(get_uniform_superposition(10));
    statevectors.push_back(get_non_stabiliser_state(10));

    for (const auto &vector : statevectors) {
        auto start_time = high_resolution_clock::now();

        bool result = stim::stabilizer_state_vector_to_circuit(vector);

        auto end_time = high_resolution_clock::now();

        auto time_taken = duration_cast<microseconds>(end_time - start_time);

        std::cout << "Stim Time taken: " << time_taken.count() << " microseconds" << std::endl;
        std::cout << "Result was: " << result << std::endl;
    }

    std::cout << "###\nBenchmarking Our Method\n###\n";

    for (const auto &vector : statevectors) {
        auto start_time = high_resolution_clock::now();

        fst::Stabiliser_From_Vector_Convertor convertor (vector);

        bool result = convertor.is_stabiliser_state;

        auto end_time = high_resolution_clock::now();

        auto time_taken = duration_cast<microseconds>(end_time - start_time);

        std::cout << "Our Time taken: " << time_taken.count() << " microseconds" << std::endl;
        std::cout << "Result was: " << result << std::endl;
    }
    return 0;
}
#include <iostream>
#include <vector>
#include <complex>

#include "./stim/circuit_vs_amplitudes.h"

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
    std::cout << "Testing Stim\n";

    bool basis_result = stim::stabilizer_state_vector_to_circuit(get_last_basis_vector(10));
    bool uniform_result = stim::stabilizer_state_vector_to_circuit(get_uniform_superposition(10));
    bool non_stab_result = stim::stabilizer_state_vector_to_circuit(get_non_stabiliser_state(10));

    std::cout << "Basis Vector: " << basis_result << std::endl;
    std::cout << "Uniform: " << uniform_result << std::endl;
    std::cout << "Non Stabiliser: " << non_stab_result << std::endl;
    return 0;
}
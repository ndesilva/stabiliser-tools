#include <iostream>
#include <vector>
#include <complex>

#include "./stim/circuit_vs_amplitudes.h"
#include "f2_helper/f2_helper.h"

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

	statevector.front() *= -1;

    return statevector;
}

int main() {
    std::cout << "Testing Stim\n";

    const bool basis_result = stim::stabilizer_state_vector_to_circuit(get_last_basis_vector(10));
    const bool uniform_result = stim::stabilizer_state_vector_to_circuit(get_uniform_superposition(10));
    const bool non_stab_result = stim::stabilizer_state_vector_to_circuit(get_non_stabiliser_state(10));

    std::cout << "Basis Vector: " << basis_result << std::endl;
    std::cout << "Uniform: " << uniform_result << std::endl;
    std::cout << "Non Stabiliser: " << non_stab_result << std::endl;
    return 0;
}
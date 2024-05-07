#include <iostream>
#include <vector>
#include <complex>
#include <chrono>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <array>
#include <numeric>

#include "stim/circuit_vs_amplitudes.h"
#include "stabiliser_state_from_statevector.h"
#include "f2_helper.h"

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
    static constexpr std::size_t REPETITIONS = 1000;

    static constexpr std::size_t MIN_N = 3;
    static constexpr std::size_t MAX_N = 10;

    static constexpr std::size_t STATEVECTOR_TYPES = 3;

    static constexpr std::size_t number_ns = MAX_N - MIN_N + 1;

    std::vector<std::size_t> qubit_numbers (number_ns);
    std::iota(qubit_numbers.begin(), qubit_numbers.end(), MIN_N);

    std::array<bool, STATEVECTOR_TYPES> statevector_is_stabiliser{};
    std::array<std::vector<double>, STATEVECTOR_TYPES> ratios;
    std::array<std::vector<std::complex<float>>, STATEVECTOR_TYPES> statevectors;

    std::vector<double> basis_vector_ratios (number_ns);
    std::vector<double> uniform_superposition_ratios (number_ns);
    std::vector<double> non_stabiliser_ratios (number_ns);

    std::cout << "Beginning benchmarking" << std::endl;

	for ( const std::size_t number_qubits : qubit_numbers )
	{
		statevectors[ 0 ] = get_last_basis_vector( number_qubits );
		statevector_is_stabiliser[ 0 ] = true;

		statevectors[ 1 ] = get_uniform_superposition( number_qubits );
		statevector_is_stabiliser[ 1 ] = true;

        statevectors[ 2 ] = get_non_stabiliser_state( number_qubits );
        statevector_is_stabiliser[ 2 ] = false;

        namespace chrono = std::chrono;
        using double_ms = chrono::duration<double, std::micro>;
    
        for ( std::size_t i = 0; i < STATEVECTOR_TYPES; i++) {
            const std::vector<std::complex<float>>& statevector = statevectors[i];
            const bool expected_result = statevector_is_stabiliser[i];
            
            double_ms stim_total_time{ 0 };
            for ( std::size_t j = 0; j < REPETITIONS; j++) {
                const auto stim_start_time = chrono::high_resolution_clock::now();

                const bool result = stim::stabilizer_state_vector_to_circuit(statevector);

                const auto stim_end_time = chrono::high_resolution_clock::now();
                stim_total_time += duration_cast<chrono::microseconds>(stim_end_time - stim_start_time);

                if (result != expected_result)
                {
                    std::cout << "STIM error, expected: " << expected_result << "got: " << result << "on index: " << i << std::endl;
                }
            }

			const double_ms stim_average_time = stim_total_time / REPETITIONS;

            double_ms our_total_time{ 0 };
            for ( std::size_t j = 0; j < REPETITIONS; j++) {
                const auto our_start_time = chrono::high_resolution_clock::now();

                const bool result = fst::is_stabiliser_state(statevector);

                const auto our_end_time = chrono::high_resolution_clock::now();
                our_total_time += duration_cast<chrono::microseconds>(our_end_time - our_start_time);

                if (result != expected_result)
                {
                    std::cout << "OUR error, expected: " << expected_result << " got " << result << " on index " << i << " qubit number: " << number_qubits << std::endl;
                }
            }

            const double_ms our_average_time = our_total_time/REPETITIONS;

            const double time_improvement = stim_average_time / our_average_time;
			ratios[ i ].push_back( time_improvement );
        }
    }

    static constexpr std::string_view file_name = "./benchmarking/stabiliser_bechmark_data.txt";
    
    std::cout << "done getting data, saving to file: " << file_name << std::endl;

    std::ofstream data_file ( file_name.data() );

	for ( const std::vector<double>& ratio : ratios )
	{
		data_file << "[";
		std::ranges::copy( ratio, std::ostream_iterator<double>( data_file, ", " ) );
		data_file << "]\n";
	}

    return 0;
}
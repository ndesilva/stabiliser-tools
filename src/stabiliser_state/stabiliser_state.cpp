#include "stabiliser_state.h"
#include "f2_helper/f2_helper.h"

#include <math.h>
#include <iostream>

using namespace fst;

Stabiliser_State::Stabiliser_State(const int number_qubits, const int dim) 
    : number_qubits(number_qubits), dim(dim) {}

Stabiliser_State::Stabiliser_State(const int number_qubits) 
    : number_qubits(number_qubits) {
        dim = number_qubits;
    }

std::vector<std::complex<float>> Stabiliser_State::get_state_vector() const {
    int support_size = 1 << dim;
    std::complex<float> factor = global_phase / float(sqrt(support_size));
    
    std::vector<std::complex<float>> state_vector (1 << number_qubits, 0);

    for(int vector_index = 0; vector_index < support_size; vector_index++) {
        int whole_space_index = shift ^ evaluate_basis_expansion(vector_index);
        std::complex<float> phase = get_phase(vector_index);
        state_vector[whole_space_index] = factor*phase;
    }

    return state_vector;  
}

size_t Stabiliser_State::evaluate_basis_expansion(const size_t vector_index) const {
    int result = 0;

    for(int j = 0; j < dim; j++){
        result ^= basis_vectors[j]*( (1<<j & vector_index) == 1 << j );
    }

    return result;
}

std::complex<float> Stabiliser_State::get_phase(const size_t vector_index) const {
    float real_linear = sign_f2_dot_product(vector_index, real_linear_part);
    std::complex<float> imag_linear = imag_f2_dot_product(vector_index, imaginary_part); 
    float real_quadratic = evaluate_quadratic_form(vector_index, quadratic_form);

    return real_linear*imag_linear*real_quadratic;
}
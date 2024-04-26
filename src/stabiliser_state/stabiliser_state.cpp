#include "stabiliser_state.h"
#include "f2_helper/f2_helper.h"

#include <math.h>
#include <iostream>

using namespace fst;

Stabiliser_State::Stabiliser_State(int number_qubits, int dim) 
    : number_qubits(number_qubits), dim(dim) {}

Stabiliser_State::Stabiliser_State(int number_qubits) 
    : number_qubits(number_qubits) {
        dim = number_qubits;
    }

std::vector<std::complex<float>> Stabiliser_State::get_state_vector() const {
    int support_size = 1 << dim;
    std::complex<float> factor = global_factor / float(sqrt(support_size));
    
    std::vector<std::complex<float>> state_vector (1 << number_qubits, 0);

    for(int vector_index = 0; vector_index < support_size; vector_index++) {
        int whole_space_index = shift ^ evaluate_basis_expansion(vector_index);
        std::complex<float> phase = get_phase(vector_index);
        state_vector[whole_space_index] = factor*phase;
    }

    return state_vector;  
}

int Stabiliser_State::evaluate_basis_expansion(int vector_index) const {
    int result = 0;

    for(int j = 0; j < dim; j++){
        result ^= basis_vectors[j]*( (1<<j & vector_index) == 1 << j );
    }

    return result;
}

std::complex<float> Stabiliser_State::get_phase(int vector_index) const {
    float real_linear = sign_f2_dot_product(vector_index, real_linear_part);
    std::complex<float> imag_linear = imag_f2_dot_product(vector_index, imaginary_part); 
    float real_quadratic = evaluate_quadratic_form(vector_index);

    return real_linear*imag_linear*real_quadratic;
}

int Stabiliser_State::evaluate_quadratic_form(int vector_index) const {    
    int mod2_result = 0;
    
    for(const auto &term : quadratic_form) {
        mod2_result ^= ((term & vector_index) == term);
    }

    return 1 - 2*mod2_result;
}
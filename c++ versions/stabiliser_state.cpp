#include "stabiliser_state.h"
#include "f2_helper.h"
#include <math.h>

using namespace fst;

template <int dimension>
Stabiliser_State<dimension>::Stabiliser_State(int number_qubits) 
    : number_qubits(number_qubits){
        dim = dimension;
    }

template <int dimension>
std::vector<std::complex<float>> Stabiliser_State<dimension>::get_state_vector() const {
    int support_size = 1 << dim;
    std::complex<float> factor = global_factor / sqrt(support_size);
    
    std::vector<std::complex<float>> state_vector (1 << number_qubits, 0);

    for(int vector_index = 0; vector_index < support_size; vector_index++) {
        int whole_space_index = shift ^ evaluate_basis_expansion(vector_index);
        std::complex<int> phase = get_phase(vector_index);
        state_vector[whole_space_index] = factor*phase;
    }

    return state_vector;  
}

template <int dimension>
int Stabiliser_State<dimension>::evaluate_basis_expansion(int vector_index) const {
    int result = 0;

    for(int j = 0; j < dim; j++){
        result ^= vector_basis[j]*( (1<<j & vector_index) == 1 << j )
    }

    return result;
}

template <int dimension>
std::complex<int> Stabiliser_State<dimension>::get_phase(int vector_index) const {
    int real_linear = sign_f2_dot_product(vector_index, real_linear_part);
    std::complex<int> imag_linear = imag_f2_dot_product(vector_index, imaginary_part); 
    int real_quadratic = evaluate_quadratic_form(vector_index);

    return real_linear*imag_linear*real_quadratic;
}

template <int dimension>
int Stabiliser_State<dimension>::evaluate_quadratic_form(int vector_index) const {
    int mod2_result = 0;
    
    for(const auto &term : quadratic_form) {
        mod2_result ^= (term & vector_index == term);
    }

    return 1 - 2*mod2_result;
}
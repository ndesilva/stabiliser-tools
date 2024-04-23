#include "stabiliser_state.h"

using namespace fst;

template <int dimension>
Stabiliser_State<dimension>::Stabiliser_State(int number_qubits, int dim) 
    : number_qubits(number_qubits), dim(dim) {}

template <int dimension>
std::vector<std::complex<float>> Stabiliser_State<dimension>::get_state_vector() const {
    std::vector<std::complex<float>> state_vector (1 << number_qubits, 0);

    return state_vector;  
}

template <int dimension>
std::complex<float> Stabiliser_State<dimension>::get_phase(int vector_index) const {
    return 0
}

template <int dimension>
int Stabiliser_State<dimension>::evaluate_quadratic_form(int vector_index) const {
    int result = 0;
    
    for(const auto &term : quadratic_form){
        result += (term & vector_index == term);
    }

    return result;
}

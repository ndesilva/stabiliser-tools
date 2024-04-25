#ifndef _FAST_STABILISER_STABILISER_STATE_H
#define _FAST_STABILISER_STABILISER_STATE_H

#include <complex>
#include <vector>

namespace fst {

/// The class used to represent a stabiliser state
/// 
/// The state is stored using the ideas of Dehaene & De Moore, as a 
/// an affine space, and a quadratic and linear form over that space.
/// More precisely, it is stored as a list of basis vectors for a vector space,
/// a constant vector that is added to every element of the vector space to reach,
/// the affine space, and a quadratic and linear form defined on the vector space.
/// the template is the dimension of the vector space
template <int dimension>
struct Stabiliser_State {
    int number_qubits;
    std::array<int, dimension> vector_basis;
    int dim;
    int shift;
    
    int real_linear_part;
    int imaginary_part;
    std::vector<int> quadratic_form;
    std::complex<float> global_factor;

    bool row_reduced;

    Stabiliser_State(int number_qubits, int dim);
    
    /// Return the state vector of length 2^n of the stabiliser state
    std::vector<std::complex<float>> get_state_vector() const;

private:
    std::complex<int> get_phase(int vector_index) const;
    int evaluate_quadratic_form(int vector_index) const;
    int evaluate_basis_expansion(int vector_index) const;
};

}

#include "stabiliser_state.cpp"

#endif
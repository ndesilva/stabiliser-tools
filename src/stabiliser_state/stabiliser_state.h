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
struct Stabiliser_State {
    int number_qubits;
    std::vector<size_t> basis_vectors;
    int dim;
    size_t shift;
    
    size_t real_linear_part;
    size_t imaginary_part;
    std::vector<size_t> quadratic_form;
    std::complex<float> global_phase = 1;

    bool row_reduced = false;

    Stabiliser_State(const int number_qubits, const int dim);
    explicit Stabiliser_State(const int number_qubits);
    
    /// Return the state vector of length 2^n of the stabiliser state
    std::vector<std::complex<float>> get_state_vector() const;

    /// Given vector_index, the column vector of an element of the vector
    /// space (represented as an integer) with respect to the vector basis,
    /// find its representation in the computational basis. 
    size_t evaluate_basis_expansion(const size_t vector_index) const;

    /// Given vector_index, the column vector of an element of the vector
    /// space (represented as an integer) with respect to the vector basis,
    /// find the value of the phase assigned to the corresponding state
    std::complex<float> get_phase(const size_t vector_index) const;
};

}

#endif
#ifndef _FAST_STABILISER_STABILISER_STATE_H
#define _FAST_STABILISER_STABILISER_STATE_H

#include <complex>
#include <vector>

namespace fst
{
	/// The class used to represent a stabiliser state
	///
	/// The state is stored using the ideas of Dehaene & De Moore, as a
	/// an affine space, and a quadratic and linear form over that space.
	/// More precisely, it is stored as a list of basis vectors for a vector space,
	/// a constant vector that is added to every element of the vector space to reach,
	/// the affine space, and a quadratic and linear form defined on the vector space.
	struct Stabiliser_State
	{
		std::size_t number_qubits = 0;
		std::vector<std::size_t> basis_vectors;
		std::size_t dim = 0;
		std::size_t shift = 0;

		std::size_t real_linear_part = 0;
		std::size_t imaginary_part = 0;
		std::vector<std::size_t> quadratic_form;
		std::complex<float> global_phase = 1.0;

		bool row_reduced = false;

		Stabiliser_State() = default;
		Stabiliser_State(const std::size_t number_qubits, const std::size_t dim);
		explicit Stabiliser_State(const std::size_t number_qubits);

		/// Return the state vector of length 2^n of the stabiliser state (with respect)
		/// to the computational basis)
		std::vector<std::complex<float>> get_state_vector() const;

		/// Given vector_index, the column vector of an element of the vector
		/// space (represented as an integer) with respect to the vector basis,
		/// find its representation in the computational basis.
		std::size_t evaluate_basis_expansion(const std::size_t vector_index) const;

		/// Given vector_index, the column vector of an element of the vector
		/// space (represented as an integer) with respect to the vector basis,
		/// find the value of the phase assigned to the corresponding state
		std::complex<float> get_phase(const std::size_t vector_index) const;

		bool operator==(const Stabiliser_State &other) const = default;
	};
}

#endif
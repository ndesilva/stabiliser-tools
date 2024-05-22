#ifndef _FAST_STABILISER_STABILISER_STATE_H
#define _FAST_STABILISER_STABILISER_STATE_H

#include "pauli.h"

#include <vector>
#include <complex>
#include <unordered_map>

namespace fst
{
	struct Check_Matrix;

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
		explicit Stabiliser_State(Check_Matrix &check_matrix);

		// TODO add documentation
		void set_support_from_cm(const Check_Matrix &check_matrix);
		void set_linear_and_quadratic_forms_from_cm(const Check_Matrix &check_matrix);

		/// Return the state vector of length 2^n of the stabiliser state (with respect
		/// to the computational basis)
		std::vector<std::complex<float>> get_state_vector() const;

		std::unordered_map<std::size_t, bool> get_quadratic_form_as_map() const;
		void update_real_parts_from_map(std::unordered_map<std::size_t, bool> &m_quadratic_form);

        void row_reduce_basis(std::unordered_map<std::size_t, bool> &m_quadratic_form);

        void add_vi_to_vj(const std::size_t i, const std::size_t j, const std::size_t v_i, std::unordered_map<size_t, bool> &m_quadratic_form);

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
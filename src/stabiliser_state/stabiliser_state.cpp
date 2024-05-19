#include "stabiliser_state.h"
#include "f2_helper.h"
#include "pauli.h"

#include <cmath>

namespace fst
{
	Stabiliser_State::Stabiliser_State(const std::size_t number_qubits, const std::size_t dim)
		: number_qubits(number_qubits), dim(dim)
	{
	}

	Stabiliser_State::Stabiliser_State(const std::size_t number_qubits)
		: Stabiliser_State(number_qubits, number_qubits)
	{
	}

	Check_Matrix Stabiliser_State::get_check_matrix()
	{
		std::map<std::size_t, bool> m_quadratic_form = get_quadratic_form_as_map();
		row_reduce_basis(m_quadratic_form);
		
		std::vector<Pauli> paulis;
		paulis.reserve(number_qubits);

		std::vector<int> pivot_indices; // TODO

		add_z_only_stabilisers(paulis, pivot_indices);
		add_x_stabilisers(paulis, pivot_indices, m_quadratic_form);

		return Check_Matrix(paulis);
	}

	std::vector<std::complex<float>> Stabiliser_State::get_state_vector() const
	{
		const std::size_t support_size = integral_pow_2(dim);
		const std::complex<float> factor = global_phase / float(std::sqrt(support_size));

		std::vector<std::complex<float>> state_vector(integral_pow_2(number_qubits), 0);

		for (std::size_t vector_index = 0; vector_index < support_size; vector_index++)
		{
			const std::size_t whole_space_index = shift ^ evaluate_basis_expansion(vector_index);
			const std::complex<float> phase = get_phase(vector_index);
			state_vector[whole_space_index] = factor * phase;
		}

		return state_vector;
	}

	void Stabiliser_State::row_reduce_basis(std::map<std::size_t, bool> &m_quadratic_form)
	{
		if (row_reduced) {return;}

		// TODO : implement. MUST UPDATE CLASS QUADRATIC FORM

		update_quadratic_form_from_map(m_quadratic_form);

		row_reduced = true;
	}

	std::size_t Stabiliser_State::evaluate_basis_expansion(const std::size_t vector_index) const
	{
		std::size_t result = 0;

		for (std::size_t j = 0; j < dim; j++)
		{
			result ^= basis_vectors[j] * std::popcount(integral_pow_2(j) & vector_index);
		}

		return result;
	}

	std::complex<float> Stabiliser_State::get_phase(const std::size_t vector_index) const
	{
		const float real_linear = sign_f2_dot_product(vector_index, real_linear_part);
		const std::complex<float> imag_linear = imag_f2_dot_product(vector_index, imaginary_part);
		const float real_quadratic = evaluate_quadratic_form(vector_index, std::span(quadratic_form));

		return real_linear * imag_linear * real_quadratic;
	}
}

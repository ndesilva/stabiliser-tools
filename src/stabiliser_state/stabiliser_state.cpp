#include "stabiliser_state.h"
#include "check_matrix.h"
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

	Stabiliser_State::Stabiliser_State(Check_Matrix &check_matrix)
	{
		check_matrix.row_reduce();

		number_qubits = check_matrix.number_qubits;
		dim = check_matrix.x_stabilisers.size();

		set_support_from_cm(check_matrix);
		set_linear_and_quadratic_forms_from_cm(check_matrix);

		row_reduced = true;
	}

	void Stabiliser_State::set_support_from_cm(const Check_Matrix &check_matrix)
	{
		basis_vectors.reserve(dim);

		for (const auto &pauli : check_matrix.x_stabilisers)
        {
            basis_vectors.push_back((*pauli).x_vector);
        }

		for (const auto &pauli : check_matrix.z_only_stabilisers)
        {
            std::size_t pivot_index = integral_log_2((*pauli).z_vector);
            shift |= (integral_pow_2(pivot_index) * (*pauli).sign_bit);
        }
	}

	void Stabiliser_State::set_linear_and_quadratic_forms_from_cm(const Check_Matrix &check_matrix)
	{
		for (std::size_t j = 0; j < dim; j++)
        {
            std::size_t v_j = basis_vectors[j];
            Pauli *p_j = check_matrix.x_stabilisers[j];
            std::size_t beta_j = (*p_j).z_vector;
            std::size_t imag_bit = (*p_j).imag_bit;

            imaginary_part |= integral_pow_2(j) * imag_bit;
            real_linear_part |= integral_pow_2(j) * ((*p_j).sign_bit ^ f2_dot_product(beta_j, v_j ^ shift));

            for (std::size_t i = 0; i < j; i++)
            {
                std::size_t v_i = basis_vectors[i];
                std::size_t other_imag_bit = (*check_matrix.x_stabilisers[i]).imag_bit; // TODO we are accessing the imag_bits alot, optimise?

                if (f2_dot_product(beta_j, v_i) ^ imag_bit*other_imag_bit)
                {
                    quadratic_form.push_back(integral_pow_2(i) | integral_pow_2(j));
                }
            }
        }
	}

	// Check_Matrix Stabiliser_State::get_check_matrix()
	// {
	// 	std::map<std::size_t, bool> m_quadratic_form = get_quadratic_form_as_map();
	// 	row_reduce_basis(m_quadratic_form);
		
	// 	std::vector<Pauli> paulis;
	// 	paulis.reserve(number_qubits);

	// 	std::vector<int> pivot_indices; // TODO

	// 	add_z_only_stabilisers(paulis, pivot_indices);
	// 	add_x_stabilisers(paulis, pivot_indices, m_quadratic_form);

	// 	return Check_Matrix(paulis);
	// }

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

	// void Stabiliser_State::row_reduce_basis(std::map<std::size_t, bool> &m_quadratic_form)
	// {
	// 	if (row_reduced) {return;}

	// 	// TODO : implement. MUST UPDATE CLASS QUADRATIC FORM

	// 	// update_quadratic_form_from_map(m_quadratic_form);

	// 	row_reduced = true;
	// }

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

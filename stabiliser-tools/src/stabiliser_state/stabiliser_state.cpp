#include "stabiliser_state.h"
#include "check_matrix.h"
#include "util/f2_helper.h"
#include "pauli/pauli.h"

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
		number_qubits = check_matrix.number_qubits;

		check_matrix.row_reduce();
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
		quadratic_form.reserve(dim*(dim+1)/2 + 1);
		quadratic_form[0] = 0;

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

				quadratic_form[integral_pow_2(i) | integral_pow_2(j)] = f2_dot_product(beta_j, v_i) ^ imag_bit*other_imag_bit;
            }
        }
	}

	std::vector<std::complex<float>> Stabiliser_State::get_state_vector() const
	{
		const std::size_t support_size = integral_pow_2(dim);
		std::vector<std::complex<float>> state_vector(integral_pow_2(number_qubits), 0);

		std::size_t vector_index = 0;
		bool imag_exponent = 0;
		std::size_t total_index = shift;
		std::complex<float> phase = global_phase / float(std::sqrt(support_size));

		state_vector[total_index] = phase;

		for (std::size_t iterate = 1; iterate < support_size; iterate++)
		{
			// Iterate through the Gray code
			std::size_t new_vector_index = iterate ^ (iterate >> 1);
			std::size_t flipped_bit = integral_log_2(vector_index ^ new_vector_index);

			total_index ^= basis_vectors[flipped_bit];
			float real_linear_phase_update = f_min1_pow(bit_set_at(real_linear_part, flipped_bit));

			bool new_imag_exponent = bit_set_at(imaginary_part, flipped_bit) ^ imag_exponent;
			// multiply by i if going from 1 to i, multiply by -i if going from i to 1
			std::complex<float> imaginary_phase_update {(float) 1-(imag_exponent^new_imag_exponent), (float) (imag_exponent^new_imag_exponent)*(1-2*imag_exponent)};

			bool quadratic_update_exponent = 0;

			for (std::size_t j = 0; j < dim; j++)
			{
				quadratic_update_exponent ^= (quadratic_form.at(integral_pow_2(flipped_bit) ^ integral_pow_2(j)) & bit_set_at(vector_index, j));
			}

			float quadratic_phase_update = f_min1_pow(quadratic_update_exponent);

			phase *= real_linear_phase_update * imaginary_phase_update * quadratic_phase_update;

			state_vector[total_index] = phase;
			vector_index = new_vector_index;
			imag_exponent = new_imag_exponent;
		}
		
		return state_vector;
	}

	void Stabiliser_State::row_reduce_basis()
	{
		if (row_reduced) {return;}

		for(std::size_t i = 0; i < dim; i++)
		{
			std::size_t v_i = basis_vectors[i];
			std::size_t pivot_index = integral_log_2(v_i);

			for(std::size_t j = 0; j < dim; j++)
			{
				if (i != j && bit_set_at(basis_vectors[j], pivot_index))
				{
                    add_vi_to_vj(i, j, v_i);
                }
			}
		}
		row_reduced = true;
    }

    void Stabiliser_State::add_vi_to_vj(const std::size_t i, const std::size_t j, const std::size_t v_i)
    {
        basis_vectors[j] ^= v_i;

        imaginary_part ^= integral_pow_2(j) * bit_set_at(imaginary_part, i);
		real_linear_part ^= integral_pow_2(j) * bit_set_at(real_linear_part, j);

        for (std::size_t k = 0; k < dim; k++)
        {
            quadratic_form[integral_pow_2(k) ^ integral_pow_2(j)] ^= (quadratic_form[integral_pow_2(k) ^ integral_pow_2(i)] & (k!=j));
        }

        quadratic_form[integral_pow_2(j)] ^= quadratic_form[integral_pow_2(i)];
    }
}

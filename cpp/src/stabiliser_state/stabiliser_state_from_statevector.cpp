#include "stabiliser_state_from_statevector.h"

#include "util/f2_helper.h"

#include <algorithm>
#include <optional>
#include <vector>

namespace
{
	template <bool assume_valid, bool return_state>
	auto stabiliser_from_statevector_internal(const std::span<const std::complex<float>> statevector)
		-> std::conditional_t<return_state, std::optional<fst::Stabiliser_State>, bool>
	{
		using namespace fst;

		const std::size_t state_vector_size = statevector.size();

		if (!is_power_of_2(state_vector_size))
		{
			return {};
		}

		const std::size_t number_qubits = integral_log_2(state_vector_size);
		std::size_t shift = 0;

		while (shift < state_vector_size && statevector[shift] == float(0))
		{
			++shift;
		}

		if (shift == state_vector_size)
		{
			return {};
		}

		std::vector<std::size_t> vector_space_indices;
		vector_space_indices.reserve(state_vector_size + 1);
		vector_space_indices.push_back(0);

		for (std::size_t index = shift + 1; index < state_vector_size; index++)
		{
			if (statevector[index] != float(0))
			{
				vector_space_indices.push_back(shift ^ index);
			}
		}

		const std::size_t support_size = vector_space_indices.size();

		if (!is_power_of_2(support_size))
		{
			return {};
		}

		const std::size_t dimension = integral_log_2(support_size);
		const float normalisation_factor = static_cast<float>(std::sqrt(support_size));
		const std::complex<float> first_entry = statevector[shift];
		const std::complex<float> global_phase = normalisation_factor * first_entry;

		if (std::abs(std::norm(global_phase) - 1) >= 0.125)
		{
			return {};
		}

		std::ranges::sort(vector_space_indices);

		std::vector<std::size_t> basis_vectors;
		basis_vectors.reserve(dimension);

		std::size_t real_linear_part = 0;
		std::size_t imaginary_part = 0;

		for (std::size_t j = 0; j < dimension; j++)
		{
			const std::size_t weight_one_string = integral_pow_2(j);

			const std::size_t basis_vector = vector_space_indices[weight_one_string];
			basis_vectors.push_back(basis_vector);

			std::complex<float> phase = statevector[basis_vector ^ shift] / first_entry;

			if (std::norm(phase - float(-1)) < 0.125)
			{
				real_linear_part ^= weight_one_string;
			}
			else if (std::norm(phase - std::complex<float>{0, 1}) < 0.125)
			{
				imaginary_part ^= weight_one_string;
			}
			else if (std::norm(phase - std::complex<float>{0, -1}) < 0.125)
			{
				real_linear_part ^= weight_one_string;
				imaginary_part ^= weight_one_string;
			}
			else if (std::norm(phase - float(1)) >= 0.125)
			{
				return {};
			}
		}

		std::unordered_map<std::size_t, bool> quadratic_form;
		quadratic_form.reserve(dimension * (dimension + 1)/2 + 1);

		quadratic_form[0] = 0;

		for (std::size_t j = 0; j < dimension; j++)
		{
			for (std::size_t i = j + 1; i < dimension; i++)
			{
				const std::size_t vector_index = integral_pow_2(i) | integral_pow_2(j);

				const float real_linear_eval = sign_f2_dot_product(vector_index, real_linear_part);
				const std::complex<float> imag_linear_eval = imag_f2_dot_product(vector_index, imaginary_part);
				const std::complex<float> linear_eval = real_linear_eval * imag_linear_eval;

				const std::size_t total_index = vector_space_indices[vector_index] ^ shift;

				const std::complex<float> quadratic_form_eval = statevector[total_index] / (first_entry * linear_eval);

				if (std::norm(quadratic_form_eval - float(-1)) < 0.125)
				{
					quadratic_form[vector_index] = 1;
				}
				else if(std::norm(quadratic_form_eval - float(1)) < 0.125)
				{
					quadratic_form[vector_index] = 0;
				}
				else
				{
					return {};
				}
			}
		}

		if constexpr (!assume_valid)
		{
			// TODO this duplicates alot of code from the get_state_vector() method of stabiliser_state. Fix
			std::size_t vector_index = 0;
			bool imag_exponent = 0;
			std::size_t total_index = shift;
			std::complex<float> phase = global_phase / float(std::sqrt(support_size));

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

				for (std::size_t j = 0; j < dimension; j++)
				{
					quadratic_update_exponent ^= (quadratic_form.at(integral_pow_2(flipped_bit) ^ integral_pow_2(j)) & bit_set_at(vector_index, j));
				}

				float quadratic_phase_update = f_min1_pow(quadratic_update_exponent);

				phase *= real_linear_phase_update * imaginary_phase_update * quadratic_phase_update;

				if (std::norm(phase - statevector[total_index]) >= 0.001)
				{
					return {};
				}

				vector_index = new_vector_index;
				imag_exponent = new_imag_exponent;
			}
		}

		if constexpr (return_state)
		{
			Stabiliser_State state(number_qubits, dimension);
			state.shift = shift;
			state.basis_vectors = std::move(basis_vectors);
			state.real_linear_part = real_linear_part;
			state.imaginary_part = imaginary_part;
			state.quadratic_form = std::move(quadratic_form);
			state.global_phase = global_phase;
			state.row_reduced = true;
			return state;
		}
		else
		{
			return true;
		}
	}
}

fst::Stabiliser_State fst::stabiliser_from_statevector(const std::span<const std::complex<float>> statevector, bool assume_valid)
{
	std::optional<Stabiliser_State> state = assume_valid
												? stabiliser_from_statevector_internal<true, true>(statevector)
												: stabiliser_from_statevector_internal<false, true>(statevector);

	if (!state)
	{
		throw std::invalid_argument("State was not a stabiliser state");
	}

	return *std::move(state);
}

bool fst::is_stabiliser_state(const std::span<const std::complex<float>> statevector)
{
	return stabiliser_from_statevector_internal<false, false>(statevector);
}

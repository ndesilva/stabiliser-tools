#include "stabiliser_state_from_statevector.h"

#include "f2_helper/f2_helper.h"

#include <algorithm>
#include <vector>

namespace
{
	template<bool AssumeValid>
	std::optional<fst::Stabiliser_State> make_stabalizer_internal( const std::span<const std::complex<float>> statevector )
	{
		using namespace fst;

		const std::size_t state_vector_size = statevector.size();
		const std::size_t number_qubits = integral_log_2( state_vector_size );

		if ( 1z << number_qubits != state_vector_size )
		{
			return std::nullopt;
		}

		std::size_t shift = 0;

		while ( shift < state_vector_size && statevector[ shift ] == float( 0 ) )
		{
			++shift;
		}

		if ( shift == state_vector_size )
		{
			return std::nullopt;
		}

		std::vector<std::size_t> vector_space_indices;
		vector_space_indices.reserve( state_vector_size + 1 );
		vector_space_indices.push_back( 0 );

		for ( std::size_t index = shift + 1; index < state_vector_size; index++ )
		{
			if ( statevector[ index ] != float( 0 ) )
			{
				vector_space_indices.push_back( shift ^ index );
			}
		}

		const std::size_t support_size = vector_space_indices.size();
		const std::size_t dimension = integral_log_2( support_size );

		if ( 1z << dimension != support_size )
		{
			return std::nullopt;
		}

		const float normalisation_factor = static_cast<float>( std::sqrt( support_size ) );
		const std::complex<float> first_entry = statevector[ shift ];
		const std::complex<float> global_phase = normalisation_factor * first_entry;

		if ( std::abs( std::norm( global_phase ) - 1 ) >= 0.125 )
		{
			return std::nullopt;
		}

		std::ranges::sort( vector_space_indices );

		std::vector<std::size_t> basis_vectors;
		basis_vectors.reserve( dimension );

		std::size_t real_linear_part = 0;
		std::size_t imaginary_part = 0;

		for ( std::size_t j = 0; j < dimension; j++ )
		{
			const std::size_t weight_one_string = 1z << j;

			const std::size_t basis_vector = vector_space_indices[ weight_one_string ];
			basis_vectors.push_back( basis_vector );

			std::complex<float> phase = statevector[ basis_vector ^ shift ] / first_entry;

			if ( std::norm( phase - float( -1 ) ) < 0.125 )
			{
				real_linear_part ^= weight_one_string;
			}
			else if ( std::norm( phase - std::complex<float>{0, 1} ) < 0.125 )
			{
				imaginary_part ^= weight_one_string;
			}
			else if ( std::norm( phase - std::complex<float>{0, -1} ) < 0.125 )
			{
				real_linear_part ^= weight_one_string;
				imaginary_part ^= weight_one_string;
			}
			else if ( std::norm( phase - float( 1 ) ) >= 0.125 )
			{
				return std::nullopt;
			}
		}

		std::vector<std::size_t> quadratic_form;
		quadratic_form.reserve( dimension * dimension );

		for ( std::size_t j = 0; j < dimension; j++ )
		{
			for ( std::size_t i = j + 1; i < dimension; i++ )
			{
				const std::size_t vector_index = ( 1z << i ) | ( 1z << j );

				const float real_linear_eval = static_cast<float>( sign_f2_dot_product( vector_index, real_linear_part ) );
				const std::complex<float> imag_linear_eval = imag_f2_dot_product( vector_index, imaginary_part );
				const std::complex<float> linear_eval = real_linear_eval * imag_linear_eval;

				const std::size_t total_index = vector_space_indices[ vector_index ] ^ shift;

				const std::complex<float> quadratic_form_eval = statevector[ total_index ] / ( first_entry * linear_eval );

				if ( std::norm( quadratic_form_eval - float( -1 ) ) < 0.125 )
				{
					quadratic_form.push_back( vector_index );
				}
				else if ( std::norm( quadratic_form_eval - float( 1 ) ) >= 0.125 )
				{
					return std::nullopt;
				}
			}
		}

		if constexpr ( !AssumeValid )
		{			
			const std::span<const std::size_t> form_span( quadratic_form );
			std::size_t old_vector_index = 0;
			std::size_t total_index = shift;

			for ( std::size_t i = 1; i < support_size; i++ )
			{
				// iterate through the gray code
				const std::size_t new_vector_index = i ^ ( i >> 1 );

				const std::size_t flipped_bit = integral_log_2( new_vector_index ^ old_vector_index );
				total_index ^= basis_vectors[ flipped_bit ];

				const std::complex<float> actual_phase = statevector[ total_index ];
				const float real_linear_eval = static_cast<float>( sign_f2_dot_product( new_vector_index, real_linear_part ) );
				const std::complex<float> imag_linear_eval = imag_f2_dot_product( new_vector_index, imaginary_part );
				const std::complex<float> quadratic_eval = static_cast<float>( evaluate_quadratic_form( new_vector_index, form_span ) );
				const std::complex<float> phase_eval = real_linear_eval * imag_linear_eval * quadratic_eval;

				if ( std::norm( phase_eval * first_entry - actual_phase ) >= 0.125 )
				{
					return std::nullopt;
				}

				old_vector_index = new_vector_index;
			}
		}

		Stabiliser_State state( number_qubits, dimension );
		state.shift = shift;
		state.basis_vectors = std::move( basis_vectors );
		state.real_linear_part = real_linear_part;
		state.imaginary_part = imaginary_part;
		state.quadratic_form = quadratic_form;
		state.global_phase = global_phase;
		state.row_reduced = true;

		return std::move( state );
	}
}

std::optional<fst::Stabiliser_State> fst::make_stabalizer_assume_valid( const std::span<const std::complex<float>> statevector )
{
	return make_stabalizer_internal<true>( statevector );
}

std::optional<fst::Stabiliser_State> fst::make_stabalizer_and_validate( const std::span<const std::complex<float>> statevector )
{
	return make_stabalizer_internal<false>( statevector );
}

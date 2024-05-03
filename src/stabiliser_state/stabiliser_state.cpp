#include "stabiliser_state.h"
#include "f2_helper/f2_helper.h"

#include <cmath>

namespace fst
{
	Stabiliser_State::Stabiliser_State( const std::size_t number_qubits, const std::size_t dim )
		: number_qubits( number_qubits )
		, dim( dim )
	{
	}

	Stabiliser_State::Stabiliser_State( const std::size_t number_qubits )
		: Stabiliser_State( number_qubits, number_qubits )
	{
	}

	std::vector<std::complex<float>> Stabiliser_State::get_state_vector() const
	{
		const std::size_t support_size = 1z << dim;
		const std::complex<float> factor = global_phase / float( std::sqrt( support_size ) );

		std::vector<std::complex<float>> state_vector( 1z << number_qubits, 0 );

		for ( std::size_t vector_index = 0; vector_index < support_size; vector_index++ )
		{
			const std::size_t whole_space_index = shift ^ evaluate_basis_expansion( vector_index );
			const std::complex<float> phase = get_phase( vector_index );
			state_vector[ whole_space_index ] = factor * phase;
		}

		return state_vector;
	}

	std::size_t Stabiliser_State::evaluate_basis_expansion( const std::size_t vector_index ) const
	{
		std::size_t result = 0;

		for ( std::size_t j = 0; j < dim; j++ )
		{
			result ^= basis_vectors[ j ] * ( ( 1z << j & vector_index ) == 1z << j );
		}

		return result;
	}

	std::complex<float> Stabiliser_State::get_phase( const std::size_t vector_index ) const
	{
		const float real_linear = static_cast<float>( sign_f2_dot_product( vector_index, real_linear_part ) );
		const std::complex<float> imag_linear = imag_f2_dot_product( vector_index, imaginary_part );
		const float real_quadratic = static_cast<float>( evaluate_quadratic_form( vector_index, std::span( quadratic_form ) ) );

		return real_linear * imag_linear * real_quadratic;
	}
}

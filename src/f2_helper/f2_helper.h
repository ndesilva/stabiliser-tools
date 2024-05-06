#ifndef _FAST_STABILISER_F2_HELPER_H
#define _FAST_STABILISER_F2_HELPER_H

#include <bit>
#include <complex>
#include <concepts>
#include <span>

namespace fst
{
	/// Return the number of binary digits of number, i.e. the ceiling of log_2(number)
	template<std::unsigned_integral T>
	constexpr int integral_log_2( const T number ) noexcept
	{
		return std::bit_width( number ) - 1;
	}

	/// Returns 2^exponent
	template<std::unsigned_integral T>
	constexpr T integral_pow_2( const T exponent ) noexcept
	{
		return T( 1 ) << exponent;
	}

	template<std::unsigned_integral T>
	constexpr bool is_power_of_2( const T value ) noexcept
	{
		return std::has_single_bit( value );
	}

	/// Gives the F_2 inner product between 2 F_2 vectors (represented
	/// as integers)
	template<std::unsigned_integral T>
	constexpr int f2_dot_product( const T x, const T y ) noexcept
	{
		const T product = x & y;
		const int hamming_weight = std::popcount( product );
		return hamming_weight % 2;
	}

	/// Gives (-1)^(x.y), where . is the F_2 inner product between 2
	/// F_2 vectors (represented as integers)
	template<std::unsigned_integral T>
	constexpr int sign_f2_dot_product( const T x, const T y ) noexcept
	{
		return 1 - 2 * f2_dot_product( x, y );
	}

	/// Gives (i)^(x.y), where . is the F_2 inner product between 2
	/// F_2 vectors (represented as integers)
	template<std::unsigned_integral T>
	constexpr std::complex<float> imag_f2_dot_product( const T x, const  T y ) noexcept
	{
		const int dot_product = f2_dot_product( x, y );
		return { static_cast<float>( 1 - dot_product ), static_cast<float>( dot_product ) };
	}

	/// Given vector_index, the column vector of an element of the vector
	/// space (represented as an integer) with respect to the vector basis, and
	/// Q a quadratic form with respect to the same basis (represented as a list of 
	/// coefficients, i.e. 101 corresponds to x_0 x_1), find the value of (-1)^Q(vector_index).
	template<std::unsigned_integral T, std::size_t Extent>
	constexpr int evaluate_quadratic_form( const T vector_index, const std::span<const T, Extent> quadratic_form ) noexcept
	{
		int mod2_result = 0;

		for ( const T term : quadratic_form )
		{
			mod2_result ^= ( ( term & vector_index ) == term );
		}

		return 1 - 2 * mod2_result;
	}
}

#endif
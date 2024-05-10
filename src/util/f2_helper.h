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

	/// Returns 2^exponent for an integer exponent
	template<std::unsigned_integral T>
	constexpr T integral_pow_2( const T exponent ) noexcept
	{
		return T( 1 ) << exponent;
	}

	/// Checks if the integer is a power of 2
	template<std::unsigned_integral T>
	constexpr bool is_power_of_2( const T number ) noexcept
	{
		return std::has_single_bit( number );
	}

	/// for an integer number (which should be 0 or 1), returns the
	/// negation of the number (i.e. 1-it) as a float
	template<std::unsigned_integral T>
	constexpr float float_not(const T number ) noexcept
	{
		return static_cast<float>(number ^ 1);
	}

	/// for an integer number (which should be 0 or 1), returns
	/// (-1)^number as a float
	template<std::unsigned_integral T>
	constexpr float min1_pow(const T number ) noexcept
	{
		return 1-2*static_cast<float>(number);
	}

	/// Gives the F_2 inner product between 2 F_2 vectors (represented
	/// as integers)
	template<std::unsigned_integral T>
	constexpr unsigned int f2_dot_product( const T x, const T y ) noexcept
	{
		const T product = x & y;
		const int hamming_weight = std::popcount( product );
		return hamming_weight % 2;
	}

	/// Gives (-1)^(x.y), where . is the F_2 inner product between 2
	/// F_2 vectors (represented as integers)
	template<std::unsigned_integral T>
	constexpr float sign_f2_dot_product( const T x, const T y ) noexcept
	{
		return min1_pow( f2_dot_product( x, y ) );
	}

	/// Gives (i)^(x.y), where . is the F_2 inner product between 2
	/// F_2 vectors (represented as integers)
	template<std::unsigned_integral T>
	constexpr std::complex<float> imag_f2_dot_product( const T x, const  T y ) noexcept
	{
		const unsigned int dot_product = f2_dot_product( x, y );
		return { float_not( dot_product ), static_cast<float>( dot_product ) };
	}

	/// Given vector_index, the column vector of an element of the vector
	/// space (represented as an integer) with respect to the vector basis, and
	/// Q a quadratic form with respect to the same basis (represented as a list of 
	/// coefficients, i.e. 101 corresponds to x_0 x_1), find the value of (-1)^Q(vector_index).
	template<std::unsigned_integral T, std::size_t Extent>
	constexpr float evaluate_quadratic_form( const T vector_index, const std::span<const T, Extent> quadratic_form ) noexcept
	{
		unsigned int mod2_result = 0;

		for ( const T term : quadratic_form )
		{
			mod2_result ^= ( ( term & vector_index ) == term );
		}

		return min1_pow( mod2_result );
	}
}

#endif
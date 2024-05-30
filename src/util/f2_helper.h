#ifndef _FAST_STABILISER_F2_HELPER_H
#define _FAST_STABILISER_F2_HELPER_H

#include <bit>
#include <complex>
#include <concepts>
#include <span>

namespace fst
{
	/// Return the number of binary digits of number, i.e. the ceiling of log_2(number)
	template <std::unsigned_integral T>
	constexpr int integral_log_2(const T number) noexcept
	{
		return std::bit_width(number) - 1;
	}

	/// Returns 2^exponent for an integer exponent
	template <std::unsigned_integral T>
	constexpr T integral_pow_2(const T exponent) noexcept
	{
		return T(1) << exponent;
	}

	/// Checks if the integer is a power of 2
	template <std::unsigned_integral T>
	constexpr bool is_power_of_2(const T number) noexcept
	{
		return std::has_single_bit(number);
	}

	/// Returns the binary digit of number at index index (zero indexed)
	/// (as a bool)
	template <std::unsigned_integral T, std::unsigned_integral U>
	constexpr bool bit_set_at(const T number, const U index) noexcept
	{
		return (number >> index) & 1;
	}

	/// For an integer number (which should be 0 or 1), returns the
	/// negation of the number (i.e. 1-it) as a float
	template <std::unsigned_integral T>
	constexpr float float_not(const T number) noexcept
	{
		return static_cast<float>(number ^ 1);
	}

	/// For an integer number (which should be 0 or 1), returns
	/// (-1)^number as an integer
	template <std::unsigned_integral T>
	constexpr int min1_pow(const T number) noexcept
	{
		return static_cast<int>(1 - 2 * number);
	}

	/// For an integer number (which should be 0 or 1), returns
	/// (-1)^number as a float
	template <std::unsigned_integral T>
	constexpr float f_min1_pow(const T number) noexcept
	{
		return static_cast<float>(min1_pow(number));
	}

	/// Gives the F_2 inner product between 2 F_2 vectors (represented
	/// as integers)
	template <std::unsigned_integral T>
	constexpr unsigned int f2_dot_product(const T x, const T y) noexcept
	{
		const T product = x & y;
		const int hamming_weight = std::popcount(product);
		return hamming_weight % 2;
	}

	/// Gives (-1)^(x.y), where . is the F_2 inner product between 2
	/// F_2 vectors (represented as integers)
	template <std::unsigned_integral T>
	constexpr float sign_f2_dot_product(const T x, const T y) noexcept
	{
		return f_min1_pow(f2_dot_product(x, y));
	}

	/// Gives (i)^(x.y), where . is the F_2 inner product between 2
	/// F_2 vectors (represented as integers)
	template <std::unsigned_integral T>
	constexpr std::complex<float> imag_f2_dot_product(const T x, const T y) noexcept
	{
		const unsigned int dot_product = f2_dot_product(x, y);
		return {float_not(dot_product), static_cast<float>(dot_product)};
	}
}

#endif
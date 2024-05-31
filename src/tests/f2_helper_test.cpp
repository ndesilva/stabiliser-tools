#include <catch2/catch_test_macros.hpp>

#include "test_util.h"
#include "f2_helper.h"
#include "platform.h"

#include <array>

using namespace fst;
using namespace test;

TEST_CASE("integral log 2", "[f2 helper]")
{
	const std::size_t power_of_2 = 0b10000;
	const std::size_t non_power = 0b1001;

	REQUIRE(integral_log_2(power_of_2) == 4);
	REQUIRE(integral_log_2(non_power) == 3);
	REQUIRE(integral_log_2(0uz) == -1);
}

TEST_CASE("integral pow 2", "[f2 helper]")
{
	const std::size_t test_exponent = 3;

	REQUIRE(integral_pow_2(test_exponent) == 8);
	REQUIRE(integral_pow_2(0uz) == 1);
}

TEST_CASE("is power of 2", "[f2 helper]")
{
	const std::size_t power_of_2 = 0b10000;
	const std::size_t non_power = 0b1001;

	REQUIRE(is_power_of_2(power_of_2));
	REQUIRE(is_power_of_2(1uz));
	REQUIRE_FALSE(is_power_of_2(non_power));
	REQUIRE_FALSE(is_power_of_2(0uz));
}

TEST_CASE("bit set at", "[f2 helper]")
{
	const std::size_t x = 0b10101;

	REQUIRE(bit_set_at(x, 0uz));
	REQUIRE(bit_set_at(x, 4uz));
	REQUIRE_FALSE(bit_set_at(x, 1uz));
	REQUIRE_FALSE(bit_set_at(x, 3uz));
}

TEST_CASE("float not", "[f2 helper]")
{
	REQUIRE(float_not(0uz) == 1);
	REQUIRE(float_not(1uz) == 0);
	
	std::size_t non_standard_input = 0b10;

	REQUIRE(float_not(non_standard_input) == 3);

	// This should run without error
	std::complex<float> product = I * float_not(0uz);
	REQUIRE(product == I);
}

TEST_CASE("minus one to the power", "[f2 helper]")
{
	std::size_t non_standard_input = 0b10;

	SECTION("integer version")
	{
		REQUIRE(min1_pow(0uz) == 1);
		REQUIRE(min1_pow(1uz) == -1);

		REQUIRE(min1_pow(non_standard_input) == -3);
	}

	SECTION("float version")
	{
		REQUIRE(f_min1_pow(0uz) == 1);
		REQUIRE(f_min1_pow(1uz) == -1);

		REQUIRE(f_min1_pow(non_standard_input) == -3);

		// This should run without error
		std::complex<float> product = I * f_min1_pow(0uz);
		REQUIRE(product == I);
	}
}

TEST_CASE("mod 2 product", "[f2 helper]")
{
	const std::size_t x = 0b11010;
	const std::size_t y = 0b10011;
	const std::size_t z = 0b00011;

	REQUIRE(f2_dot_product(x, y) == 0);
	REQUIRE(f2_dot_product(x, z) == 1);
}

TEST_CASE("signed mod 2 product", "[f2 helper]")
{
	const std::size_t x = 0b110011;
	const std::size_t y = 0b000011;
	const std::size_t z = 0b110101;

	REQUIRE(sign_f2_dot_product(x, y) == 1);
	REQUIRE(sign_f2_dot_product(x, z) == -1);
}

TEST_CASE("imaginary mod 2 product", "[f2 helper]")
{
	const std::size_t x = 0b101011;
	const std::size_t y = 0b000011;
	const std::size_t z = 0b101101;

	REQUIRE(imag_f2_dot_product(x, y) == 1.0f);
	REQUIRE(imag_f2_dot_product(x, z) == I);
}

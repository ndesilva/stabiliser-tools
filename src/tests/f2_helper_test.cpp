#include <catch2/catch_test_macros.hpp>

#include "f2_helper.h"

#include <array>

using namespace fst;

static constexpr std::complex<float> i = {0, 1};

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
	REQUIRE(imag_f2_dot_product(x, z) == i);
}

TEST_CASE("evaluate quadratic form", "[f2 helper]")
{
	SECTION("dimension 3")
	{
		const std::array<std::size_t, 2> quadratic_form = {3, 6}; // x_0 x_1 + x_1 x_2

		const std::size_t point_1 = 0b011;
		const std::size_t point_2 = 0b111;

		REQUIRE(evaluate_quadratic_form(point_1, std::span(quadratic_form)) == -1);
		REQUIRE(evaluate_quadratic_form(point_2, std::span(quadratic_form)) == 1);
	}

	SECTION("dimension 4")
	{
		const std::array<std::size_t, 3> quadratic_form = {3, 6, 9}; // x_0 x_1 + x_1 x_2 + x_0 x_3
		const std::size_t point_1 = 0b1111;
		const std::size_t point_2 = 0b1011;

		REQUIRE(evaluate_quadratic_form(point_1, std::span(quadratic_form)) == -1);
		REQUIRE(evaluate_quadratic_form(point_2, std::span(quadratic_form)) == 1);
	}
}
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include "stabiliser_state.h"
#include "check_matrix.h"
#include "f2_helper.h"

#include <array>

using namespace Catch::Matchers;
using namespace fst;

static constexpr std::complex<float> i = {0, 1};

namespace
{
	// TODO extract to test util file
	std::unordered_map<std::size_t, bool> get_quadratic_from_from_vector(std::size_t dimension, std::vector<std::size_t> non_zero_coeffs)
	{
		std::unordered_map<std::size_t, bool> quadratic_form;
		quadratic_form.reserve(dimension * (dimension + 1)/2);
		quadratic_form[0] = 0;

		for (std::size_t j = 0; j < dimension; j++)
		{
			for (std::size_t k = 0; k < dimension; k++)
			{
				quadratic_form[integral_pow_2(j) ^ integral_pow_2(k)] = 0;
			}
		}

		for (const auto coeff : non_zero_coeffs)
		{
			quadratic_form[coeff] = 1;
		}

		return quadratic_form;
	}

	Check_Matrix get_check_matrix()
	{
		std::vector<Pauli> paulis;

		paulis.push_back(Pauli(5, 0b00011, 0b00001, 0, 1));
		paulis.push_back(Pauli(5, 0b00101, 0b00010, 0, 0));
		paulis.push_back(Pauli(5, 0b00110, 0b00100, 1, 1));
		paulis.push_back(Pauli(5, 0b10000, 0b01000, 1, 0));
		paulis.push_back(Pauli(5, 0b10110, 0b00100, 0, 1));

		return Check_Matrix(paulis);
	}

	Stabiliser_State get_stabiliser_state()
	{
		Stabiliser_State state(5, 3);
		state.basis_vectors = {3, 5, 16};
		state.shift = 0;
		
		state.real_linear_part = 5;
		state.imaginary_part = 1;
		state.quadratic_form = get_quadratic_from_from_vector(3, {3});

		state.row_reduced = true;

		return state;
	}

	TEST_CASE("stabiliser state from check matrix", "[stabiliser state]")
	{
		Check_Matrix check_matrix = get_check_matrix();
        Stabiliser_State expected_stabiliser_state = get_stabiliser_state();

        Stabiliser_State stabiliser_state(check_matrix);

        REQUIRE(stabiliser_state == expected_stabiliser_state);
	}

	TEST_CASE("generate state vector", "[stabiliser state]")
	{
		SECTION("dimension 0, 1 qubit")
		{
			Stabiliser_State state(1, 0);

			state.quadratic_form = get_quadratic_from_from_vector(0, {});
			state.real_linear_part = 1;
			state.imaginary_part = 0;

			state.basis_vectors = {0};
			state.shift = 1;

			const std::array<std::complex<float>, 2> expected_vector{0, 1};

			REQUIRE_THAT(state.get_state_vector(), RangeEquals(expected_vector));
		}

		SECTION("dimension 1, 1 qubit")
		{
			Stabiliser_State state(1);

			state.quadratic_form = get_quadratic_from_from_vector(1, {});
			state.real_linear_part = 1;
			state.imaginary_part = 0;

			state.basis_vectors = {1};
			state.shift = 1;

			std::array<std::complex<float>, 2> expected_vector{-1 / std::sqrt(2.0f), 1 / std::sqrt(2.0f)};

			REQUIRE_THAT(state.get_state_vector(), RangeEquals(expected_vector));
		}

		SECTION("dimension 2, 3 qubits")
		{
			Stabiliser_State state(3, 2);

			state.quadratic_form = get_quadratic_from_from_vector(2, {3}); // x_0 x_1
			state.real_linear_part = 1; // x_0
			state.imaginary_part = 2;	// x_1
			state.global_phase = {0, 1};

			state.basis_vectors = {0b110, 0b001};

			state.shift = 0b100;

			std::array<std::complex<float>, 8> expected_vector{0, 0, -0.5f * i, -.5, 0.5f * i, -.5, 0, 0};

			REQUIRE_THAT(state.get_state_vector(), RangeEquals(expected_vector));
		}
	}

	TEST_CASE("row reduce", "[stabiliser state]")
	{
		Stabiliser_State state(5, 3);
		
		state.quadratic_form = get_quadratic_from_from_vector(3, {3});
		state.real_linear_part = 1;
		state.imaginary_part = 1;
		state.basis_vectors = {0b01101, 0b10110, 0b11010};
		state.shift = 0b00010;

		std::vector<std::size_t> row_reduced_basis {0b01100, 0b10110, 0b00001};
		std::vector<std::complex<float>> statevector = state.get_state_vector();

		state.row_reduce_basis();

		REQUIRE_THAT(state.basis_vectors, RangeEquals(row_reduced_basis));
		REQUIRE_THAT(statevector, RangeEquals(state.get_state_vector()));
	}
}
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include "stabiliser_state/stabiliser_state.h"
#include "stabiliser_state/check_matrix.h"
#include "util/f2_helper.h"
#include "test_util.h"

#include <array>

using namespace Catch::Matchers;
using namespace fst;
using namespace test;

namespace
{
	Check_Matrix get_unreduced_check_matrix()
	{
		std::vector<Pauli> paulis;

		paulis.push_back(Pauli(5, 0b00011, 0b00001, 0, 1));
		paulis.push_back(Pauli(5, 0b00101, 0b00010, 0, 0));
		paulis.push_back(Pauli(5, 0b00110, 0b00100, 1, 1));
		paulis.push_back(Pauli(5, 0b10000, 0b01000, 1, 0));
		paulis.push_back(Pauli(5, 0b10110, 0b00100, 0, 1));

		return Check_Matrix(paulis);
	}

	TEST_CASE("stabiliser state from check matrix", "[stabiliser state]")
	{
		Check_Matrix check_matrix = get_unreduced_check_matrix();
        Stabiliser_State stabiliser_state(check_matrix);

		std::vector<std::complex<float>> statevector = stabiliser_state.get_state_vector();

		for (const auto pauli : check_matrix.paulis)
        {
            REQUIRE_THAT(matrix_vector_mult(pauli.get_matrix(), statevector), RangeEquals(statevector));
        }
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

			const std::array<std::complex<float>, 2> expected_vector {0, 1};

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

			std::array<std::complex<float>, 8> expected_vector{0, 0, -0.5f * I, -.5, 0.5f * I, -.5, 0, 0};

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

		std::vector<std::complex<float>> old_statevector = state.get_state_vector();

		state.row_reduce_basis();

		REQUIRE(state.row_reduced);
		REQUIRE(vectors_row_reduced(state.basis_vectors));
		REQUIRE_THAT(old_statevector, RangeEquals(state.get_state_vector()));
	}
}
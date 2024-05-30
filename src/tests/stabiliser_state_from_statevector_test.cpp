#include <catch2/catch_test_macros.hpp>

#include "stabiliser_state_from_statevector.h"
#include "f2_helper.h"
#include <vector>
#include <array>

using namespace fst;

static constexpr std::complex<float> i = {0, 1};

namespace
{
	std::unordered_map<std::size_t, bool> get_quadratic_from_from_vector(std::size_t dimension, std::vector<std::size_t> non_zero_coeffs)
	{
		std::unordered_map<std::size_t, bool> quadratic_form;
		quadratic_form.reserve(dimension * (dimension + 1)/2);
		quadratic_form[0] = 0;

		for (std::size_t i = 0; i < dimension; i++)
		{
			for (std::size_t j = 0; j < dimension; j++)
			{
				quadratic_form[integral_pow_2(i) ^ integral_pow_2(j)] = 0;
			}
		}

		for (const auto coeff : non_zero_coeffs)
		{
			quadratic_form[coeff] = 1;
		}

		return quadratic_form;
	}

	std::array<std::complex<float>, 8> get_three_qubit_stabiliser_statevector()
	{
		return {0, 0, -0.5, 0.5f * i, 0.5, 0.5f * i, 0, 0};
	}

	std::array<std::complex<float>, 32> get_five_qubit_stabiliser_statevector()
	{
		std::array<std::complex<float>, 32> statevector;
		statevector.fill(0);

		const float root_8 = std::sqrt(8.0f);

		statevector[1] = 1 / root_8;
		statevector[7] = -i / root_8;
		statevector[8] = 1 / root_8;
		statevector[14] = -i / root_8;
		statevector[17] = -1 / root_8;
		statevector[23] = i / root_8;
		statevector[24] = 1 / root_8;
		statevector[30] = -i / root_8;

		return statevector;
	}

	std::array<std::complex<float>, 32> get_five_qubit_stabiliser_statevector_with_phase(std::complex<float> phase)
	{
		std::array<std::complex<float>, 32> statevector = get_five_qubit_stabiliser_statevector();
		
		for (auto &elt : statevector)
		{
			elt *= phase;
		}

		return statevector;
	}
}

TEST_CASE("testing correct stabiliser states", "[statevector -> stabiliser state]")
{
	SECTION("1 qubit, dimension 0")
	{
		const std::vector<std::complex<float>> statevector = {0, 1};

		REQUIRE(is_stabiliser_state(statevector));
	}

	SECTION("3 qubits, dimension 2")
	{
		const std::array statevector = get_three_qubit_stabiliser_statevector();

		REQUIRE(is_stabiliser_state(statevector));
	}

	SECTION("5 qubits, dimension 3")
	{
		const std::array statevector = get_five_qubit_stabiliser_statevector();

		REQUIRE(is_stabiliser_state(statevector));
	}

	SECTION("5 qubits, dimension 3, global factor")
	{
		const std::complex<float> global_phase(1 / std::sqrt(2.0f), 1 / std::sqrt(2.0f));
		const std::array statevector = get_five_qubit_stabiliser_statevector_with_phase(global_phase);

		REQUIRE(is_stabiliser_state(statevector));
	}
}

TEST_CASE("testing incorrect stabiliser states", "[statevector -> stabiliser state]")
{
	SECTION("all zero state (dimension 2)")
	{
		const std::array<std::complex<float>, 4> statevector{0, 0, 0, 0};

		REQUIRE_FALSE(is_stabiliser_state(statevector));
	}

	SECTION("non power of 2 sized vector (dimension 25)")
	{
		std::array<std::complex<float>, 25> statevector;
		statevector.fill(0.2f);

		REQUIRE_FALSE(is_stabiliser_state(statevector));
	}

	SECTION("non-normalised")
	{
		std::array statevector = get_three_qubit_stabiliser_statevector();

		for (auto &elt : statevector)
		{
			elt *= 2;
		}

		REQUIRE_FALSE(is_stabiliser_state(statevector));
	}

	SECTION("incorrect support size")
	{
		std::array statevector = get_three_qubit_stabiliser_statevector();

		statevector[0] = 1; // add an element to the support

		REQUIRE_FALSE(is_stabiliser_state(statevector));
	}

	SECTION("support not affine space")
	{
		std::array statevector = get_three_qubit_stabiliser_statevector();

		statevector[1] = 0; // remove element from the support
		statevector[0] = 1; // keep support size the same, but now not affine

		REQUIRE_FALSE(is_stabiliser_state(statevector));
	}

	SECTION("invalid basis vector entry")
	{
		std::array statevector = get_five_qubit_stabiliser_statevector();

		statevector[7] = 2; // invalid entry for the first basis vector

		REQUIRE_FALSE(is_stabiliser_state(statevector));
	}

	SECTION("invalid weight two vector entry")
	{
		std::array statevector = get_five_qubit_stabiliser_statevector();

		statevector[14] = -2.0f * i; // invalid entry for e_1 + e_2

		REQUIRE_FALSE(is_stabiliser_state(statevector));
	}

	SECTION("inconsistent remaining entry")
	{
		std::array statevector = get_five_qubit_stabiliser_statevector();

		statevector[30] = -1.0f * i; // inconsistent entry for e_1 + e_2 + e_3 (should be i)

		REQUIRE_FALSE(is_stabiliser_state(statevector));
	}

	SECTION("uniform with incorrect phase")
	{
		const std::size_t vector_length = fst::integral_pow_2(5u);
		const float normalisation_factor = static_cast<float>(1.0f / std::sqrt(vector_length));

		std::vector<std::complex<float>> statevector(vector_length, normalisation_factor);

		statevector.back() *= -1;

		REQUIRE_FALSE(is_stabiliser_state(statevector));
	}
}

TEST_CASE("incorrect stabiliser state flagged as stabiliser state", "[statevector -> stabiliser state]")
{
	std::array statevector = get_five_qubit_stabiliser_statevector();

	statevector[30] = -i; // inconsistent entry for e_1 + e_2 + e_3 (should be i)

	// check this doesn't throw invalid argument, despite not being a stabiliser state
	stabiliser_from_statevector(statevector, true);
}

Stabiliser_State get_expected_five_qubit_stabiliser_state(std::complex<float> global_phase)
{
	Stabiliser_State state(5, 3);

	state.basis_vectors = {6,9,16};
	state.shift = 1;

	state.real_linear_part = 5;
	state.imaginary_part = 1;
	state.quadratic_form = get_quadratic_from_from_vector(3, {6});
	state.global_phase = global_phase;

	state.row_reduced = true;
	return state;
}

TEST_CASE("get stabiliser state", "[statevector -> stabiilser state]")
{
	SECTION("stabiliser input, dimension 5")
	{
		const std::complex<float> global_phase(1 / std::sqrt(2.0f), 1 / std::sqrt(2.0f));
		const std::array statevector = get_five_qubit_stabiliser_statevector_with_phase(global_phase);

		const Stabiliser_State expected_state = get_expected_five_qubit_stabiliser_state(global_phase);

		Stabiliser_State state = stabiliser_from_statevector(statevector);

		// allow for some tolerance
		if (std::abs(state.global_phase - global_phase) <= 0.01)
		{
			state.global_phase = global_phase;
		}

		REQUIRE(expected_state == state);
	}

	SECTION("non-stabiliser input")
	{
		std::array statevector = get_five_qubit_stabiliser_statevector();

		statevector[1] = {1, 1}; // invalid entry

		REQUIRE_THROWS_AS(stabiliser_from_statevector(statevector), std::invalid_argument);
	}
}
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include <vector>
#include <array>

#include "stabiliser_state/stabiliser_state_from_statevector.h"
#include "util/f2_helper.h"
#include "test_util.h"

using namespace Catch::Matchers;
using namespace fst;
using namespace test;

namespace
{
	std::vector<std::complex<float>> get_three_qubit_stabiliser_statevector()
	{
		return {0, 0, -0.5, 0.5f * I, 0.5, 0.5f * I, 0, 0};
	}

	std::vector<std::complex<float>> get_five_qubit_stabiliser_statevector()
	{
		std::vector<std::complex<float>> statevector(32, 0);

		const float root_8 = std::sqrt(8.0f);

		statevector[1] = 1 / root_8;
		statevector[7] = -I / root_8;
		statevector[8] = 1 / root_8;
		statevector[14] = -I / root_8;
		statevector[17] = -1 / root_8;
		statevector[23] = I / root_8;
		statevector[24] = 1 / root_8;
		statevector[30] = -I / root_8;

		return statevector;
	}

	std::vector<std::complex<float>> get_five_qubit_stabiliser_statevector_with_phase(const std::complex<float> phase)
	{
		std::vector<std::complex<float>> statevector = get_five_qubit_stabiliser_statevector();
		
		for (auto &elt : statevector)
		{
			elt *= phase;
		}

		return statevector;
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
			const std::vector<std::complex<float>> statevector = get_three_qubit_stabiliser_statevector();

			REQUIRE(is_stabiliser_state(statevector));
		}

		SECTION("5 qubits, dimension 3")
		{
			const std::vector<std::complex<float>> statevector = get_five_qubit_stabiliser_statevector();

			REQUIRE(is_stabiliser_state(statevector));
		}

		SECTION("5 qubits, dimension 3, global factor")
		{
			const std::complex<float> global_phase(1 / std::sqrt(2.0f), 1 / std::sqrt(2.0f));
			const std::vector<std::complex<float>> statevector = get_five_qubit_stabiliser_statevector_with_phase(global_phase);

			REQUIRE(is_stabiliser_state(statevector));
		}
	}

	TEST_CASE("testing incorrect stabiliser states", "[statevector -> stabiliser state]")
	{
		SECTION("all zero state (dimension 2)")
		{
			const std::vector<std::complex<float>> statevector {0, 0, 0, 0};

			REQUIRE_FALSE(is_stabiliser_state(statevector));
		}

		SECTION("non power of 2 sized vector (dimension 25)")
		{
			std::vector<std::complex<float>> statevector (25, 0.2f);

			REQUIRE_FALSE(is_stabiliser_state(statevector));
		}

		SECTION("non-normalised")
		{
			std::vector<std::complex<float>> statevector = get_three_qubit_stabiliser_statevector();

			for (auto &elt : statevector)
			{
				elt *= 2;
			}

			REQUIRE_FALSE(is_stabiliser_state(statevector));
		}

		SECTION("incorrect support size")
		{
			std::vector<std::complex<float>> statevector = get_three_qubit_stabiliser_statevector();

			statevector[0] = 1; // add an element to the support

			REQUIRE_FALSE(is_stabiliser_state(statevector));
		}

		SECTION("support not affine space")
		{
			std::vector<std::complex<float>> statevector = get_three_qubit_stabiliser_statevector();

			statevector[1] = 0; // remove element from the support
			statevector[0] = 1; // keep support size the same, but now not affine

			REQUIRE_FALSE(is_stabiliser_state(statevector));
		}

		SECTION("invalid basis vector entry")
		{
			std::vector<std::complex<float>> statevector = get_five_qubit_stabiliser_statevector();

			statevector[7] = 2; // invalid entry for the first basis vector

			REQUIRE_FALSE(is_stabiliser_state(statevector));
		}

		SECTION("invalid weight two vector entry")
		{
			std::vector<std::complex<float>> statevector = get_five_qubit_stabiliser_statevector();

			statevector[14] = -2.0f * I; // invalid entry for e_1 + e_2

			REQUIRE_FALSE(is_stabiliser_state(statevector));
		}

		SECTION("inconsistent remaining entry")
		{
			std::vector<std::complex<float>> statevector = get_five_qubit_stabiliser_statevector();

			statevector[30] = -1.0f * I; // inconsistent entry for e_1 + e_2 + e_3 (should be i)

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
		std::vector<std::complex<float>> statevector = get_five_qubit_stabiliser_statevector();

		statevector[30] = -I; // inconsistent entry for e_1 + e_2 + e_3 (should be i)

		// check this doesn't throw invalid argument, despite not being a stabiliser state
		stabiliser_from_statevector(statevector, true);
	}

	TEST_CASE("returning stabiliser state", "[statevector -> stabiilser state]")
	{
		SECTION("stabiliser input, dimension 5")
		{
			const std::complex<float> global_phase(1 / std::sqrt(2.0f), 1 / std::sqrt(2.0f));
			const std::vector<std::complex<float>> statevector = get_five_qubit_stabiliser_statevector_with_phase(global_phase);

			Stabiliser_State state = stabiliser_from_statevector(statevector);

			REQUIRE_THAT(state.get_state_vector(), RangeEquals(statevector));
		}

		SECTION("non-stabiliser input")
		{
			std::vector<std::complex<float>> statevector = get_five_qubit_stabiliser_statevector();

			statevector[30] = {1, 1}; // invalid entry

			REQUIRE_THROWS_AS(stabiliser_from_statevector(statevector), std::invalid_argument);

			// Check doesn't throw
			stabiliser_from_statevector(statevector, true);
		}
	}
}
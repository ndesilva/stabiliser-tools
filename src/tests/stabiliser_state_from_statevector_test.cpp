#include <catch2/catch_test_macros.hpp>

#include "stabiliser_state_from_statevector.h"
#include "f2_helper.h"
#include <vector>
#include <array>

using namespace fst;

namespace
{
	std::array<std::complex<float>, 8> get_three_qubit_stabiliser_statevector()
	{
		return { 0, 0, -0.5, { 0, 0.5 }, 0.5, { 0, 0.5 }, 0, 0 };
	}

	std::array<std::complex<float>, 32> get_five_qubit_stabiliser_statevector()
	{
		std::array<std::complex<float>, 32> statevector;
		statevector.fill( 0 );

		const float root_8 = std::sqrt( 8.0f );

		statevector[ 1 ] = 1.0f / root_8;
		statevector[ 7 ] = { 0, -1 / root_8 };
		statevector[ 8 ] = 1 / root_8;
		statevector[ 14 ] = { 0, -1 / root_8 };
		statevector[ 17 ] = -1 / root_8;
		statevector[ 23 ] = { 0, 1 / root_8 };
		statevector[ 24 ] = 1 / root_8;
		statevector[ 30 ] = { 0, -1 / root_8 };

		return statevector;
	}
}

TEST_CASE( "testing correct stabiliser states", "[statevector -> stabiliser state]" )
{
	SECTION( "1 qubit, dimension 0" )
	{
		const std::vector<std::complex<float>> statevector = { 0, 1 };

		REQUIRE( is_stabiliser_state( statevector ) );
	}

	SECTION( "3 qubits, dimension 2" )
	{
		const std::array statevector = get_three_qubit_stabiliser_statevector();

		REQUIRE( is_stabiliser_state( statevector ) );
	}

	SECTION( "5 qubits, dimension 3" )
	{
		const std::array statevector = get_five_qubit_stabiliser_statevector();

		REQUIRE( is_stabiliser_state( statevector ) );
	}

	SECTION( "5 qubits, dimension 3, global factor" )
	{
		std::array statevector = get_five_qubit_stabiliser_statevector();

		const std::complex<float> global_phase( 1 / std::sqrt( 2.0f ), 1 / std::sqrt( 2.0f ) );

		for ( auto& elt : statevector )
		{
			elt *= global_phase;
		}

		REQUIRE( is_stabiliser_state( statevector ) );
	}
}

TEST_CASE( "testing incorrect stabiliser states", "[statevector -> stabiliser state]" )
{
	SECTION( "all zero state (dimension 2)" )
	{
		const std::array<std::complex<float>, 4> statevector{ 0,0,0,0 };

		REQUIRE_FALSE( is_stabiliser_state( statevector ) );
	}

	SECTION( "non power of 2 sized vector (dimension 25)" )
	{
		std::array<std::complex<float>, 25> statevector;
		statevector.fill( 0.2f );

		REQUIRE_FALSE( is_stabiliser_state( statevector ) );
	}

	SECTION( "non-normalised" )
	{
		std::array statevector = get_three_qubit_stabiliser_statevector();

		for ( auto& elt : statevector )
		{
			elt *= 2;
		}

		REQUIRE_FALSE( is_stabiliser_state( statevector ) );
	}

	SECTION( "incorrect support size" )
	{
		std::array statevector = get_three_qubit_stabiliser_statevector();

		statevector[ 0 ] = 1; // add an element to the support

		REQUIRE_FALSE( is_stabiliser_state( statevector ) );
	}

	SECTION( "support not affine space" )
	{
		std::array statevector = get_three_qubit_stabiliser_statevector();

		statevector[ 1 ] = 0; // remove element from the support
		statevector[ 0 ] = 1; // keep support size the same, but now not affine

		REQUIRE_FALSE( is_stabiliser_state( statevector ) );
	}

	SECTION( "invalid basis vector entry" )
	{
		std::array statevector = get_five_qubit_stabiliser_statevector();

		statevector[ 7 ] = 2; //invalid entry for the first basis vector

		REQUIRE_FALSE( is_stabiliser_state( statevector ) );
	}

	SECTION( "invalid weight two vector entry" )
	{
		std::array statevector = get_five_qubit_stabiliser_statevector();

		statevector[ 14 ] = { 0,-2 }; // invalid entry for e_1 + e_2

		REQUIRE_FALSE( is_stabiliser_state( statevector ) );
	}

	SECTION( "inconsistent remaining entry" )
	{
		std::array statevector = get_five_qubit_stabiliser_statevector();

		statevector[ 30 ] = { 0,-1 }; // inconsistent entry for e_1 + e_2 + e_3 (should be i)

		REQUIRE_FALSE( is_stabiliser_state( statevector ) );
	}

	SECTION( "uniform with incorrect phase" ) {
		const std::size_t vector_length = fst::integral_pow_2( 5u );
		const float normalisation_factor = static_cast<float>( 1.0f / std::sqrt( vector_length ) );
		
		std::vector<std::complex<float>> statevector (vector_length, normalisation_factor);

		statevector.back() *= -1;

		REQUIRE_FALSE( is_stabiliser_state( statevector ) );
	}
}

TEST_CASE( "incorrect stabiliser state flagged as stabiliser state", "[statevector -> stabiliser state]" )
{
	std::array statevector = get_five_qubit_stabiliser_statevector();

	statevector[ 30 ] = { 0,-1 }; // inconsistent entry for e_1 + e_2 + e_3 (should be i)

	// check this doesn't throw invalid argument, despite not being a stabiliser state
	stabiliser_from_statevector(statevector, true);
}

TEST_CASE( "get stabiliser state", "[statevector -> stabiilser state]" )
{
	SECTION( "stabiliser input, dimension 5" )
	{
		const std::array statevector = get_five_qubit_stabiliser_statevector();

		const Stabiliser_State state = stabiliser_from_statevector( statevector );

		const std::vector<size_t> expected_basis{ 6,9,16 };
		const std::vector<size_t> expected_quadratic_form{ 6 };

		REQUIRE( state.number_qubits == 5 );
		REQUIRE( state.basis_vectors == expected_basis );
		REQUIRE( state.shift == 1 );
		REQUIRE( state.real_linear_part == 5 );
		REQUIRE( state.imaginary_part == 1 );
		REQUIRE( state.quadratic_form == expected_quadratic_form );
		REQUIRE( std::norm( state.global_phase - float( 1 ) ) <= 0.001 );
		REQUIRE( state.row_reduced == true );
	}

	SECTION( "non-stabiliser input" )
	{
		std::array statevector = get_five_qubit_stabiliser_statevector();

		statevector[ 1 ] = { 1,1 }; // invalid entry
		
		REQUIRE_THROWS_AS( stabiliser_from_statevector( statevector ), std::invalid_argument );
	}
}
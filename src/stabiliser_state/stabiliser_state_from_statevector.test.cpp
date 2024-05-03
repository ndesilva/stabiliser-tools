#include <catch2/catch_test_macros.hpp>

#include "stabiliser_state_from_statevector.h"
#include <vector>
#include <array>
#include <iostream>

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
		const std::vector<std::complex<float>> statevector = { 0,1 };

		const Stabiliser_From_Vector_Convertor convertor( statevector );

		REQUIRE( convertor.is_stabiliser_state );
	}

	SECTION( "3 qubits, dimension 2" )
	{
		const std::array statevector = get_three_qubit_stabiliser_statevector();

		const Stabiliser_From_Vector_Convertor convertor( statevector );

		REQUIRE( convertor.is_stabiliser_state );
	}

	SECTION( "5 qubits, dimension 3" )
	{
		const std::array statevector = get_five_qubit_stabiliser_statevector();

		const Stabiliser_From_Vector_Convertor convertor( statevector );

		REQUIRE( convertor.is_stabiliser_state );
	}

	SECTION( "5 qubits, dimension 3, global factor" )
	{
		std::array statevector = get_five_qubit_stabiliser_statevector();

		const std::complex<float> global_phase( 1 / std::sqrt( 2.0f ), 1 / std::sqrt( 2.0f ) );

		for ( auto& elt : statevector )
		{
			elt *= global_phase;
		}

		const Stabiliser_From_Vector_Convertor convertor( statevector );

		REQUIRE( convertor.is_stabiliser_state );
	}
}

TEST_CASE( "testing incorrect stabiliser states", "[statevector -> stabiilser state]" )
{
	SECTION( "all zero state (dimension 2)" )
	{
		std::vector<std::complex<float>> statevector( 4, 0 );

		const Stabiliser_From_Vector_Convertor convertor( statevector );

		REQUIRE_FALSE( convertor.is_stabiliser_state );
	}

	SECTION( "non power of 2 sized vector (dimension 25)" )
	{
		std::vector<std::complex<float>> statevector( 25, 0.2f );

		const Stabiliser_From_Vector_Convertor convertor( statevector );

		REQUIRE_FALSE( convertor.is_stabiliser_state );
	}

	SECTION( "non-normalised" )
	{
		std::array statevector = get_three_qubit_stabiliser_statevector();

		for ( auto& elt : statevector )
		{
			elt *= 2;
		}

		Stabiliser_From_Vector_Convertor convertor( statevector );

		REQUIRE_FALSE( convertor.is_stabiliser_state );
	}

	SECTION( "incorrect support size" )
	{
		std::array statevector = get_three_qubit_stabiliser_statevector();

		statevector[ 0 ] = 1; // add an element to the support

		const Stabiliser_From_Vector_Convertor convertor( statevector );

		REQUIRE_FALSE( convertor.is_stabiliser_state );
	}

	SECTION( "support not affine space" )
	{
		std::array statevector = get_three_qubit_stabiliser_statevector();

		statevector[ 1 ] = 0; // remove element from the support
		statevector[ 0 ] = 1; // keep support size the same, but now not affine

		const Stabiliser_From_Vector_Convertor convertor( statevector );

		REQUIRE_FALSE( convertor.is_stabiliser_state );
	}

	SECTION( "invalid basis vector entry" )
	{
		std::array statevector = get_five_qubit_stabiliser_statevector();

		statevector[ 7 ] = 2; //invalid entry for the first basis vector

		const Stabiliser_From_Vector_Convertor convertor( statevector );

		REQUIRE_FALSE( convertor.is_stabiliser_state );
	}

	SECTION( "invalid weight two vector entry" )
	{
		std::array statevector = get_five_qubit_stabiliser_statevector();

		statevector[ 14 ] = { 0,-2 }; // invalid entry for e_1 + e_2

		const Stabiliser_From_Vector_Convertor convertor( statevector );

		REQUIRE_FALSE( convertor.is_stabiliser_state );
	}

	SECTION( "inconsistent remaining entry" )
	{
		std::array statevector = get_five_qubit_stabiliser_statevector();

		statevector[ 30 ] = { 0,-1 }; // inconsistent entry for e_1 + e_2 + e_3 (should be i)

		const Stabiliser_From_Vector_Convertor convertor( statevector );

		REQUIRE_FALSE( convertor.is_stabiliser_state );
	}
}

TEST_CASE( "testing incorrect stabiliser state flagged as stabiliser state", "[statevector -> stabiilser state]" )
{
	std::array statevector = get_five_qubit_stabiliser_statevector();

	statevector[ 30 ] = { 0,-1 }; // inconsistent entry for e_1 + e_2 + e_3 (should be i)

	const Stabiliser_From_Vector_Convertor convertor( statevector, true );

	REQUIRE( convertor.is_stabiliser_state );
}

TEST_CASE( "get stabiliser state", "[statevector -> stabiilser state]" )
{
	SECTION( "stabiliser input, dimension 5" )
	{
		const std::array statevector = get_five_qubit_stabiliser_statevector();

		const Stabiliser_From_Vector_Convertor convertor( statevector, true );
		const Stabiliser_State state = convertor.get_stabiliser_state();

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

		const Stabiliser_From_Vector_Convertor convertor( statevector, true );

		REQUIRE_THROWS_AS( convertor.get_stabiliser_state(), std::invalid_argument );
	}
}
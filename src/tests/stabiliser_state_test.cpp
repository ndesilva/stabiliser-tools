#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include "stabiliser_state.h"

#include <array>

using namespace Catch::Matchers;
using namespace fst;

TEST_CASE( "evaluate basis expansion", "[stabiliser state]" )
{
	SECTION( "dimension 3" )
	{
		Stabiliser_State state( 3, 3 );

		state.basis_vectors = { 0b001,0b010,0b100 };

		const std::size_t vector_index = 0b101;
		const std::size_t expected_vector = 0b101;

		REQUIRE( state.evaluate_basis_expansion( vector_index ) == expected_vector );
	}

	SECTION( "dimension 4" )
	{
		Stabiliser_State state( 5, 4 );

		state.basis_vectors = { 0b10001, 0b01100, 0b00001, 0b11010 };

		const std::size_t vector_index = 0b1011;
		const std::size_t expected_vector = 0b00111;

		REQUIRE( state.evaluate_basis_expansion( vector_index ) == expected_vector );
	}
}

TEST_CASE( "get phase", "[stabiliser state]" )
{
	SECTION( "dimension 3" )
	{
		Stabiliser_State state( 3 );

		state.quadratic_form = { 3 }; // x_0 x_1
		state.real_linear_part = 1; // x_0
		state.imaginary_part = 2;   // x_1

		const std::size_t point_1 = 0b011;
		const std::size_t point_2 = 0b101;

		REQUIRE( state.get_phase( point_1 ) == std::complex<float>( 0, 1 ) );
		REQUIRE( state.get_phase( point_2 ) == -1.0f );
	}

	SECTION( "dimension 4" )
	{
		Stabiliser_State state( 3 );

		state.quadratic_form = { 3, 12 }; // x_0 x_1 + x_2 x_3
		state.real_linear_part = 1;     // x_0
		state.imaginary_part = 7;       // x_0 + x_1 + x_2

		const std::size_t point_1 = 0b0011;
		const std::size_t point_2 = 0b0001;

		REQUIRE( state.get_phase( point_1 ) == 1.0f );
		REQUIRE( state.get_phase( point_2 ) == std::complex<float>( 0, -1 ) );
	}
}

TEST_CASE( "generate state vector", "[stabiliser state]" )
{
	SECTION( "dimension 0, 1 qubit" )
	{
		Stabiliser_State state( 1, 0 );

		state.quadratic_form = {};
		state.real_linear_part = 1;
		state.imaginary_part = 0;

		state.basis_vectors = { 0 };
		state.shift = 1;

		const std::array<std::complex<float>, 2> expected_vector{ 0, 1 };

		REQUIRE_THAT( state.get_state_vector(), RangeEquals( expected_vector ) );
	}

	SECTION( "dimension 1, 1 qubit" )
	{
		Stabiliser_State state( 1 );

		state.quadratic_form = {};
		state.real_linear_part = 1;
		state.imaginary_part = 0;

		state.basis_vectors = { 1 };
		state.shift = 1;

		std::array<std::complex<float>, 2> expected_vector{ -1 / std::sqrt( 2.0f ), 1 / std::sqrt( 2.0f ) };

		REQUIRE_THAT( state.get_state_vector(), RangeEquals( expected_vector ) );
	}

	SECTION( "dimension 2, 3 qubits" )
	{
		Stabiliser_State state( 3, 2 );

		state.quadratic_form = { 3 }; // x_0 x_1
		state.real_linear_part = 1; // x_0
		state.imaginary_part = 2;   // x_1
		state.global_phase = { 0, 1 };

		state.basis_vectors = { 0b110, 0b001 };

		state.shift = 0b100;

		std::array<std::complex<float>, 8> expected_vector{ 0, 0, {0,-.5}, -.5, {0, .5}, -.5, 0, 0 };

		REQUIRE_THAT( state.get_state_vector(), RangeEquals( expected_vector ) );
	}
}
#include <catch2/catch_test_macros.hpp>

#include "stabiliser_state.h"

using namespace fst;

TEST_CASE("evaluate basis expansion", "[stabiliser state]"){
    SECTION("dimension 3") {
        Stabiliser_State<3> state(3);

        state.vector_basis[0] = 1; // 001
        state.vector_basis[1] = 2; // 010
        state.vector_basis[2] = 4; // 100

        int vector_index = 5;    // 101
        int expected_vector = 5; // 101

        REQUIRE(state.evaluate_basis_expansion(vector_index) == expected_vector);
    }

    SECTION("dimension 4") {
        Stabiliser_State<4> state(5);

        state.vector_basis[0] = 17; // 10001
        state.vector_basis[1] = 12; // 01100
        state.vector_basis[2] = 1 ; // 00001
        state.vector_basis[3] = 26; // 11010

        int vector_index = 11;   // 1011
        int expected_vector = 7; // 00111

        REQUIRE(state.evaluate_basis_expansion(vector_index) == expected_vector);
    }
}
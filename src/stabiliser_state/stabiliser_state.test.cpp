#include <catch2/catch_test_macros.hpp>

#include "stabiliser_state.h"
#include <vector>

using namespace fst;

TEST_CASE("evaluate basis expansion", "[stabiliser state]") {
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

TEST_CASE("evaluate quadratic form", "[stabiliser state]") {
    SECTION("dimension 3") {
        Stabiliser_State<3> state(3);

        state.quadratic_form = {3, 6}; // x_0 x_1 + x_1 x_2

        int point_1 = 3; // 011
        int point_2 = 7; // 111

        REQUIRE(state.evaluate_quadratic_form(point_1) == -1);
        REQUIRE(state.evaluate_quadratic_form(point_2) == 1);
    }

     SECTION("dimension 4") {
        Stabiliser_State<4> state(4);

        state.quadratic_form = {3, 6, 9}; // x_0 x_1 + x_1 x_2 + x_0 x_3
        int point_1 = 15; // 1111
        int point_2 = 11; // 1011  

        REQUIRE(state.evaluate_quadratic_form(point_1) == -1);
        REQUIRE(state.evaluate_quadratic_form(point_2) == 1);
    }
}

TEST_CASE("get phase", "[stabiliser state]") {
    SECTION("dimension 3") {
        Stabiliser_State<3> state(3);

        state.quadratic_form = {3}; // x_0 x_1
        state.real_linear_part = 1; // x_0
        state.imaginary_part = 2;   // x_1

        int point_1 = 3; // 011
        int point_2 = 5; // 101

        REQUIRE(state.get_phase(point_1) == std::complex<float> (0,1));
        REQUIRE(state.get_phase(point_2) == -1);
    }

    SECTION("dimension 4") {
        Stabiliser_State<4> state(3);

        state.quadratic_form = {3, 12}; // x_0 x_1 + x_2 x_3
        state.real_linear_part = 1;     // x_0
        state.imaginary_part = 7;       // x_0 + x_1 + x_2

        int point_1 = 3;  // 0011
        int point_2 = 1;  // 0001

        REQUIRE(state.get_phase(point_1) == 1);
        REQUIRE(state.get_phase(point_2) == std::complex<float> (0,-1));
    }
}

TEST_CASE("generate state vector", "[stabiliser state]"){
    SECTION("dimension 1") {
        Stabiliser_State<1> state(1);

        state.quadratic_form = {};
        state.real_linear_part = 1;
        state.imaginary_part = 0;

        state.vector_basis[0] = 0;
        state.shift = 1;

        std::vector<std::complex<float>> expected_vector{0, -1};

        REQUIRE(state.get_state_vector() == expected_vector);
    }
}
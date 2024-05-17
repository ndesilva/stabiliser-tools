#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include <iostream>

#include "check_matrix.h"
#include "pauli.h"

using namespace Catch::Matchers;
using namespace fst;

static constexpr std::complex<float> i = {0, 1};

namespace
{
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

    TEST_CASE("categorise paulis", "[check_matrix]")
    {
        Pauli pauli1 = Pauli(2, 0b10, 0b01, 0, 0); // X Z
        Pauli pauli2 = Pauli(2, 0b00, 0b01, 0, 0); // 1 Z

        std::vector<Pauli> paulis = {pauli1, pauli2};

        Check_Matrix check_matrix(paulis);

        REQUIRE_FALSE(check_matrix.row_reduced);
        REQUIRE(*check_matrix.x_stabilisers[0] == pauli1);
        REQUIRE(*check_matrix.z_only_stabilisers[0] == pauli2);
    }

    Check_Matrix get_row_reduced_check_matrix()
    {
        std::vector<Pauli> paulis;

        paulis.push_back(Pauli(5, 0b00011, 0b00001, 0, 1));
        paulis.push_back(Pauli(5, 0b00101, 0b00010, 0, 0));
        paulis.push_back(Pauli(5, 0b00000, 0b00111, 0, 0));
        paulis.push_back(Pauli(5, 0b10000, 0b01000, 1, 0));
        paulis.push_back(Pauli(5, 0b00000, 0b01000, 0, 0));

        return Check_Matrix(paulis);
    }

    TEST_CASE("row reduce", "[check_matrix]")
    {
        Check_Matrix starting_check_matrix = get_check_matrix();
        Check_Matrix row_reduced_check_matrix = get_row_reduced_check_matrix();

        starting_check_matrix.row_reduce();

        REQUIRE_THAT(starting_check_matrix.paulis, RangeEquals(row_reduced_check_matrix.paulis));

        for(std::size_t i = 0; i < row_reduced_check_matrix.x_stabilisers.size(); i++){
            REQUIRE(*row_reduced_check_matrix.x_stabilisers[i] == *starting_check_matrix.x_stabilisers[i]);
        }

        for(std::size_t i = 0; i < row_reduced_check_matrix.z_only_stabilisers.size(); i++){
            REQUIRE(*row_reduced_check_matrix.z_only_stabilisers[i] == *starting_check_matrix.z_only_stabilisers[i]);
        }
    }

    Stabiliser_State get_stabiliser_state()
    {
        Stabiliser_State state(5, 3);
        state.basis_vectors = {3, 5, 16};
        state.shift = 0;
        
        state.real_linear_part = 5;
        state.imaginary_part = 1;
        state.quadratic_form = {3};

        state.row_reduced = true;

        return state;
    }

    TEST_CASE("get stabiliser state") {
        Check_Matrix check_matrix = get_check_matrix();
        Stabiliser_State expected_stabiliser_state = get_stabiliser_state();

        Stabiliser_State stabiliser_state = check_matrix.get_stabiliser_state();

        REQUIRE(stabiliser_state == expected_stabiliser_state);
    }

    std::vector<std::complex<float>> get_state_vector()
    {
        std::vector<std::complex<float>> state_vector (32, 0);
        const float root_8 = std::sqrt(8.0f);

        state_vector[0] = 1 / root_8;
        state_vector[3] = -i / root_8;
        state_vector[5] = 1 / root_8;
        state_vector[6] = i / root_8;
        state_vector[16] = -1 / root_8;
        state_vector[19] = i / root_8;
        state_vector[21] = -1 / root_8;
        state_vector[22] = -i / root_8;

        return state_vector;
    }

    TEST_CASE("get state vector") {
        Check_Matrix check_matrix = get_check_matrix();
        std::vector<std::complex<float>> expected_state_vector = get_state_vector();

        std::vector<std::complex<float>> state_vector = check_matrix.get_state_vector();

        REQUIRE(state_vector == expected_state_vector);
    }
}
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include <iostream>

#include "pauli/pauli.h"
#include "clifford/clifford.h"
#include "test_util.h"

using namespace Catch::Matchers;
using namespace fst;
using namespace test;

namespace
{

    std::vector<std::vector<std::complex<float>>> get_three_qubit_clifford_matrix()
    {
        return {{0.5f*I, 0.5f*I, 0.0, 0.0, -0.5f*I, -0.5f*I, 0.0, 0.0},
            {0.0, 0.0, 0.5f*I, -0.5f*I, 0.0, 0.0, -0.5f*I, 0.5f*I},
            {0.0, 0.0, 0.5f*I, 0.5f*I, 0.0, 0.0, -0.5f*I, -0.5f*I},
            {0.5f*I, -0.5f*I, 0.0, 0.0, -0.5f*I, 0.5f*I, 0.0, 0.0},
            {0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.5, 0.5},
            {0.5, -0.5, 0.0, 0.0, 0.5, -0.5, 0.0, 0.0},
            {0.5, 0.5, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0},
            {0.0, 0.0, 0.5, -0.5, 0.0, 0.0, 0.5, -0.5},
            };
    }

    TEST_CASE("get matrix", "[clifford]")
    {
        SECTION("1 qubit")
        {
            Pauli X (1, 1, 0, 0, 0);
            Pauli Z (1, 0, 1, 0, 0);

            Clifford hadamard ({X}, {Z});

            std::vector<std::vector<std::complex<float>>> hadamard_matrix = {{1/std::sqrt(2.0f), 1/std::sqrt(2.0f)},{1/std::sqrt(2.0f), -1/std::sqrt(2.0f)}};

            REQUIRE_THAT(hadamard.get_matrix(), RangeEquals(hadamard_matrix));
        }

        SECTION("3 qubits")
        {
            std::vector<Pauli> z_conjugates = {Pauli(3, 0b00011, 0b00000, 0, 0), Pauli(3, 0b00000, 0b00111, 0, 0), Pauli(3, 0b00110, 0b00100, 0, 1)};
            std::vector<Pauli> x_conjugates = {Pauli(3, 0b00000, 0b00001, 0, 0), Pauli(3, 0b00010, 0b00000, 0, 0), Pauli(3, 0b00000, 0b00100, 1, 0)};

            Clifford clifford (z_conjugates, x_conjugates, I);

            std::vector<std::vector<std::complex<float>>> expected_matrix = get_three_qubit_clifford_matrix();

            REQUIRE_THAT(clifford.get_matrix(), RangeEquals(expected_matrix));

        }
    }
}
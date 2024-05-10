#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include "pauli.h"

using namespace Catch::Matchers;
using namespace fst;

static constexpr std::complex<float> i = {0,1};

std::vector<std::vector<std::complex<float>>> get_one_qubit_pauli_matrix()
{
    std::vector<std::vector<std::complex<float>>> matrix (2, std::vector<std::complex<float>> (2, 0));
    matrix.at(0).at(1) = {0, -1};
    matrix.at(1).at(0) = {0, 1}; // (-1)(-i) X Z

    return matrix;
}

std::vector<std::vector<std::complex<float>>> get_three_qubit_pauli_matrix()
{
    std::vector<std::vector<std::complex<float>>> matrix (8, std::vector<std::complex<float>> (8, 0));

    matrix.at(0).at(6) = 1;
    matrix.at(1).at(7) = -1;
    matrix.at(2).at(4) = 1;
    matrix.at(3).at(5) = -1;
    matrix.at(4).at(2) = -1;
    matrix.at(5).at(3) = 1;
    matrix.at(6).at(0) = -1;
    matrix.at(7).at(1) = 1; // (-1) XX1 Z1Z

    return matrix;
}

TEST_CASE( "get_matrix", "[pauli]" )
{
    SECTION( "1 qubit" )
    {
        Pauli pauli(1, 1, 1, 1, 1); // Y

        auto expected_pauli_matrix = get_one_qubit_pauli_matrix();

        REQUIRE_THAT( pauli.get_matrix(), RangeEquals(expected_pauli_matrix) );
    }

    SECTION("3 qubits")
    {
        Pauli pauli(3, 0b110, 0b101, 1, 0); // i Y X Z

        auto expected_pauli_matrix = get_three_qubit_pauli_matrix();

        REQUIRE_THAT(pauli.get_matrix(), RangeEquals(expected_pauli_matrix));
    }
}

TEST_CASE("has eigenstate", "[pauli]")
{
    SECTION("+1 eigenstate")
    {
        Pauli pauli(4, 0b1100, 0b0101, 1, 1); // X Y 1 Z
        std::vector<std::complex<float>> vector {0, 1, 0, 0, 0, -i, 0, 0, 0, 1, 0, 0, 0, -i, 0, 0};

        REQUIRE( pauli.has_eigenstate( vector, 0 ) );
        REQUIRE_FALSE( pauli.has_eigenstate( vector, 1 ) );
    }

    SECTION("-1 eigenstate")
    {
        Pauli pauli(4, 0b1110, 0b0101, 1, 1); // X Y X Z
        std::vector<std::complex<float>> vector {1.0f, 0, -1.0f, 0, i, 0, -i, 0, 1.0f, 0, -1.0f, 0, i, 0, -i, 0};

        REQUIRE( pauli.has_eigenstate( vector, 1 ) );
        REQUIRE_FALSE( pauli.has_eigenstate( vector, 0 ) );
    }

    SECTION("non-sign eigenstate")
    {
        Pauli pauli(2, 0b10, 0b10, 0, 0); // -i Y 1
        std::vector<std::complex<float>> vector {1.0f, 0, i, 0};

        auto pauli_times_vector = pauli.multiply_vector(vector);
        
        for (int index = 0; index < 4; index++)
        {
            REQUIRE( pauli_times_vector.at(index) == -i*vector.at(index) );
        }

        REQUIRE_FALSE( pauli.has_eigenstate( vector, 1 ) );
        REQUIRE_FALSE( pauli.has_eigenstate( vector, 0 ) );
    }
}
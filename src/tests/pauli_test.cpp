#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include "pauli.h"
#include <iostream>
#include <numeric>

using namespace Catch::Matchers;
using namespace fst;

static constexpr std::complex<float> i = {0,1};

namespace
{

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

        SECTION("non-eigenstate")
        {
            Pauli pauli(1, 1, 0, 0, 0); // X
            std::vector<std::complex<float>> vector {1,0};

            REQUIRE_FALSE( pauli.has_eigenstate( vector, 1 ) );
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

    TEST_CASE("multiply_by_pauli_on_right", "[pauli]")
    {
        SECTION("no overlapping X and Z")
        {
            Pauli pauli (2, 0b11, 0b10, 0, 0);      // XZ X
            Pauli other_pauli(2, 0b01, 0b10, 1, 1); // Z  X

            Pauli expected_pauli(2, 0b10, 0, 1, 1); // X 1

            pauli.multiply_by_pauli_on_right(other_pauli);

            REQUIRE(pauli == expected_pauli);
        }

        SECTION("overlapping X and Z")
        {
            Pauli pauli (3, 0b110, 0b011, 1, 1);        // i X XZ Z
            Pauli other_pauli(3, 0b010, 0b101, 1, 1);   // i Z X  Z

            Pauli expected_pauli(3, 0b100, 0b110, 0, 0); // XZ Z 1

            pauli.multiply_by_pauli_on_right(other_pauli);

            REQUIRE(pauli == expected_pauli);
        }
    }

    TEST_CASE("testing commuting/anticommuting", "[pauli]")
    {
        SECTION("commuting")
        {
            Pauli pauli (3, 0b110, 0b011, 1, 1);        // i X XZ Z
            Pauli other_pauli (3, 0b010, 0b101, 1, 1);  // i Z X  Z

            REQUIRE(pauli.commutes_with(other_pauli));
            REQUIRE_FALSE(pauli.anticommutes_with(other_pauli));
        }

        SECTION("anticommuting")
        {
            Pauli pauli (3, 0b110, 0b001, 1, 1);        // i X X Z
            Pauli other_pauli (3, 0b010, 0b101, 1, 1);  // i Z X Z

            REQUIRE(pauli.anticommutes_with(other_pauli));
            REQUIRE_FALSE(pauli.commutes_with(other_pauli));
        }
    }

    TEST_CASE("testing Hermitian", "[pauli]")
    {
        SECTION("Hermitian")
        {
            Pauli pauli (3, 0b110, 0b011, 1, 1); // X Y Z

            REQUIRE(pauli.is_hermitian());
        }

        SECTION("Anti-Hermitian")
        {
            Pauli pauli (3, 0b101, 0b110, 1, 0); // i Y Z X
        }

    }
    
    TEST_CASE("multiplying a vector", "[pauli]")
    {
        SECTION("3 qubits")
        {
            Pauli pauli (3, 0b100, 0b011, 1, 0);

            std::vector<std::complex<float>> vector (8,0);
            std::iota(vector.begin(), vector.end(), 0);

            const std::vector<std::complex<float>> expected_product {-4.0f, 5.0f, 6.0f, -7.0f, .0f,  1.0f, 2.0f, -3.0f};

            auto product = pauli.multiply_vector(vector);

            REQUIRE_THAT(product, RangeEquals(expected_product));
        }

        SECTION("5 qubits")
        {
            Pauli pauli (5, 0b10000, 0b11111, 1, 1);

            std::vector<std::complex<float>> vector (32);
            
            for (int index = 0; index < 32; index++)
            {
                auto flt = (float) index;
                vector.at(index) = flt + 2.0f*flt*i;
            }

            const std::vector<std::complex<float>> expected_product {32.0f + -16.0f*i, -34.0f + 17.0f*i, -36.0f + 18.0f*i, 38.0f + -19.0f*i, -40.0f + 20.0f*i, 42.0f + -21.0f*i, 44.0f + -22.0f*i, -46.0f + 23.0f*i, -48.0f + 24.0f*i, 50.0f + -25.0f*i, 52.0f + -26.0f*i, -54.0f + 27.0f*i, 56.0f + -28.0f*i, -58.0f + 29.0f*i, -60.0f + 30.0f*i, 62.0f + -31.0f*i, 0.0f + 0.0f*i, 2.0f + -1.0f*i, 4.0f + -2.0f*i, -6.0f + 3.0f*i, 8.0f + -4.0f*i, -10.0f + 5.0f*i, -12.0f + 6.0f*i, 14.0f + -7.0f*i, 16.0f + -8.0f*i, -18.0f + 9.0f*i, -20.0f + 10.0f*i, 22.0f + -11.0f*i, -24.0f + 12.0f*i, 26.0f + -13.0f*i, 28.0f + -14.0f*i, -30.0f + 15.0f*i};

            auto product = pauli.multiply_vector(vector);

            REQUIRE_THAT(product, RangeEquals(expected_product));
        }
    }
}
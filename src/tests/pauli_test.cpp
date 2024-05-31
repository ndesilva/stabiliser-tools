#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include "pauli.h"
#include "test_util.h"

using namespace Catch::Matchers;
using namespace fst;
using namespace test;

namespace
{

    std::vector<std::vector<std::complex<float>>> get_one_qubit_pauli_matrix()
    {
        std::vector<std::vector<std::complex<float>>> matrix(2, std::vector<std::complex<float>>(2, 0));
        matrix.at(0).at(1) = {0, -1};
        matrix.at(1).at(0) = {0, 1}; // (-1)(-i) X Z

        return matrix;
    }

    std::vector<std::vector<std::complex<float>>> get_three_qubit_pauli_matrix()
    {
        std::vector<std::vector<std::complex<float>>> matrix(8, std::vector<std::complex<float>>(8, 0));

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

    TEST_CASE("get_matrix", "[pauli]")
    {
        SECTION("1 qubit")
        {
            Pauli pauli(1, 1, 1, 1, 1); // Y

            auto expected_pauli_matrix = get_one_qubit_pauli_matrix();

            REQUIRE_THAT(pauli.get_matrix(), RangeEquals(expected_pauli_matrix));
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
            std::vector<std::complex<float>> vector {0, 1, 0, 0, 0, -I, 0, 0, 0, 1, 0, 0, 0, -I, 0, 0};

            REQUIRE(pauli.has_eigenstate(vector, 0));
            REQUIRE_FALSE(pauli.has_eigenstate(vector, 1));
        }

        SECTION("-1 eigenstate")
        {
            Pauli pauli(4, 0b1110, 0b0101, 1, 1); // X Y X Z
            std::vector<std::complex<float>> vector {1.0f, 0, -1.0f, 0, I, 0, -I, 0, 1.0f, 0, -1.0f, 0, I, 0, -I, 0};

            REQUIRE(pauli.has_eigenstate(vector, 1));
            REQUIRE_FALSE(pauli.has_eigenstate(vector, 0));
        }

        SECTION("non-eigenstate")
        {
            Pauli pauli(1, 1, 0, 0, 0); // X
            std::vector<std::complex<float>> vector {1, 0};

            REQUIRE_FALSE(pauli.has_eigenstate(vector, 1));
            REQUIRE_FALSE(pauli.has_eigenstate(vector, 0));
        }

        SECTION("non-sign eigenstate")
        {
            Pauli pauli(2, 0b10, 0b10, 0, 0); // -i Y 1
            std::vector<std::complex<float>> vector {1.0f, 0, I, 0};

            auto pauli_times_vector = pauli.multiply_vector(vector);

            for (int index = 0; index < 4; index++)
            {
                REQUIRE(pauli_times_vector.at(index) == -I * vector.at(index));
            }

            REQUIRE_FALSE(pauli.has_eigenstate(vector, 1));
            REQUIRE_FALSE(pauli.has_eigenstate(vector, 0));
        }

        SECTION("invalid vector size")
        {
            Pauli pauli(5, 0b10000, 0b11111, 1, 1);

            std::vector<std::complex<float>> length_31_vector (31, 0);
            std::vector<std::complex<float>> length_33_vector (33, 0);

            REQUIRE_THROWS_AS(pauli.has_eigenstate(length_31_vector, 0), std::invalid_argument);
            REQUIRE_THROWS_AS(pauli.has_eigenstate(length_31_vector, 1), std::invalid_argument);
            REQUIRE_THROWS_AS(pauli.has_eigenstate(length_33_vector, 0), std::invalid_argument);
            REQUIRE_THROWS_AS(pauli.has_eigenstate(length_33_vector, 1), std::invalid_argument);
        }
    }

    TEST_CASE("multiply_by_pauli_on_right", "[pauli]")
    {
        SECTION("no overlapping X and Z")
        {
            Pauli pauli(2, 0b11, 0b10, 0, 0);       // XZ X
            Pauli other_pauli(2, 0b01, 0b10, 1, 1); // Z  X

            Pauli expected_pauli(2, 0b10, 0, 1, 1); // X 1

            pauli.multiply_by_pauli_on_right(other_pauli);

            REQUIRE(pauli == expected_pauli);
        }

        SECTION("overlapping X and Z")
        {
            Pauli pauli(3, 0b110, 0b011, 1, 1);       // i X XZ Z
            Pauli other_pauli(3, 0b010, 0b101, 1, 1); // i Z X  Z

            Pauli expected_pauli(3, 0b100, 0b110, 0, 0); // XZ Z 1

            pauli.multiply_by_pauli_on_right(other_pauli);

            REQUIRE(pauli == expected_pauli);
        }

        SECTION("extra case")
        {
            Pauli pauli(5, 0b00110, 0b00100, 1, 1);       // i  1 1 XZ X 1
            Pauli other_pauli(5, 0b00011, 0b00001, 0, 1); // -i 1 1 1  X XZ

            Pauli expected_pauli(5, 0b00101, 0b00101, 0, 0); // 1 1 1 XZ 1 XZ

            pauli.multiply_by_pauli_on_right(other_pauli);

            REQUIRE(pauli == expected_pauli);
        }

        SECTION("incompatible sizes")
        {
            Pauli five_qubit_pauli(5, 0b00110, 0b00100, 1, 1);
            Pauli four_qubit_pauli(4, 0, 0, 0, 0);
            Pauli six_qubit_pauli(6, 0, 0, 0, 0);

            REQUIRE_THROWS_AS(five_qubit_pauli.multiply_by_pauli_on_right(four_qubit_pauli), std::invalid_argument);
            REQUIRE_THROWS_AS(five_qubit_pauli.multiply_by_pauli_on_right(six_qubit_pauli), std::invalid_argument);
        }
    }

    TEST_CASE("testing commuting/anticommuting", "[pauli]")
    {
        SECTION("commuting")
        {
            Pauli pauli(3, 0b110, 0b011, 1, 1);       // i X XZ Z
            Pauli other_pauli(3, 0b010, 0b101, 1, 1); // i Z X  Z

            REQUIRE(pauli.commutes_with(other_pauli));
            REQUIRE_FALSE(pauli.anticommutes_with(other_pauli));
        }

        SECTION("anticommuting")
        {
            Pauli pauli(3, 0b110, 0b001, 1, 1);       // i X X Z
            Pauli other_pauli(3, 0b010, 0b101, 1, 1); // i Z X Z

            REQUIRE(pauli.anticommutes_with(other_pauli));
            REQUIRE_FALSE(pauli.commutes_with(other_pauli));
        }
    }

    TEST_CASE("testing Hermitian", "[pauli]")
    {
        SECTION("Hermitian")
        {
            Pauli pauli(3, 0b110, 0b011, 1, 1); // X Y Z

            REQUIRE(pauli.is_hermitian());
        }

        SECTION("Anti-Hermitian")
        {
            Pauli pauli(3, 0b101, 0b110, 1, 0); // i Y Z X

            REQUIRE_FALSE(pauli.is_hermitian());
        }
    }

    std::vector<std::complex<float>> get_three_qubit_vector()
    {
        std::vector<std::complex<float>> vector(8);
        
        for (int index = 0; index < 8; index++)
        {
            vector.at(index) = (float) index;
        }
        
        return vector;
    }

    std::vector<std::complex<float>> get_five_qubit_vector()
    {
        std::vector<std::complex<float>> vector(32);
        
        for (int index = 0; index < 32; index++)
        {
            auto flt = (float)index;
            vector.at(index) = flt + 2.0f * flt * I;
        }
        
        return vector;
    }

    TEST_CASE("multiplying a vector", "[pauli]")
    {
        SECTION("3 qubits")
        {
            Pauli pauli(3, 0b100, 0b011, 1, 0);
            std::vector<std::complex<float>> vector = get_three_qubit_vector();

            std::vector<std::complex<float>> expected_product = matrix_vector_mult(pauli.get_matrix(), vector);

            auto product = pauli.multiply_vector(vector);

            REQUIRE_THAT(product, RangeEquals(expected_product));
        }

        SECTION("5 qubits")
        {
            Pauli pauli(5, 0b10000, 0b11111, 1, 1);

            std::vector<std::complex<float>> vector = get_five_qubit_vector();
            
            std::vector<std::complex<float>> expected_product = matrix_vector_mult(pauli.get_matrix(), vector);

            auto product = pauli.multiply_vector(vector);

            REQUIRE_THAT(product, RangeEquals(expected_product));
        }

        SECTION("invlaid vector size")
        {
            Pauli pauli(5, 0b10000, 0b11111, 1, 1);

            std::vector<std::complex<float>> length_31_vector (31, 0);
            std::vector<std::complex<float>> length_33_vector (33, 0);

            REQUIRE_THROWS_AS(pauli.multiply_vector(length_31_vector), std::invalid_argument);
            REQUIRE_THROWS_AS(pauli.multiply_vector(length_33_vector), std::invalid_argument);
        }
    }
}
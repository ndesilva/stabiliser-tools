#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include <vector>

#include "util/f2_helper.h"
#include "clifford/clifford_from_matrix.h"
#include "test_util.h"

using namespace Catch::Matchers;
using namespace fst;
using namespace test;

namespace
{
    std::vector<std::vector<std::complex<float>>> get_three_qubit_clifford()
    {
        return {{0.5f*I, 0.5f*I, .0f, .0f, -0.5f*I, -0.5f*I, .0f, .0f, },
                {.0f, .0f, 0.5f*I, -0.5f*I, .0f, .0f, -0.5f*I, 0.5f*I, },
                {.0f, .0f, 0.5f*I, 0.5f*I, .0f, .0f, -0.5f*I, -0.5f*I, },
                {0.5f*I, -0.5f*I, .0f, .0f, -0.5f*I, 0.5f*I, .0f, .0f, },
                {.0f, .0f, 0.5f, 0.5f, .0f, .0f, 0.5f, 0.5f, },
                {0.5f, -0.5f, .0f, .0f, 0.5f, -0.5f, .0f, .0f, },
                {0.5f, 0.5f, .0f, .0f, 0.5f, 0.5f, .0f, .0f, },
                {.0f, .0f, 0.5f, -0.5f, .0f, .0f, 0.5f, -0.5f, }};
    }

    std::vector<std::vector<std::complex<float>>> get_three_qubit_clifford_with_phase(const std::complex<float> phase)
    {
        std::vector<std::vector<std::complex<float>>> matrix = get_three_qubit_clifford();

        for (auto &row: matrix)
        {
            for (auto &elt : row)
            {
                elt *= phase;
            }
        }

        return matrix;
    }

    std::vector<std::vector<std::complex<float>>> get_alternative_three_qubit_clifford()
    {
        return {{.0f, -0.5f*I, .0f, 0.5f*I, 0.5f, .0f, -0.5f, .0f},
                {.0f, -0.5f*I, .0f, 0.5f*I, -0.5f, .0f, 0.5f, .0f},
                {0.5f, .0f, -0.5f, .0f, .0f, -0.5f*I, .0f, 0.5f*I},
                {-0.5f, .0f, 0.5f, .0f, .0f, -0.5f*I, .0f, 0.5f*I},
                {.0f, 0.5f, .0f, 0.5f, 0.5f*I, .0f, 0.5f*I, .0f},
                {.0f, -0.5f, .0f, -0.5f, 0.5f*I, .0f, 0.5f*I, .0f},
                {0.5f*I, .0f, 0.5f*I, .0f, .0f, 0.5f, .0f, 0.5f},
                {0.5f*I, .0f, 0.5f*I, .0f, .0f, -0.5f, .0f, -0.5f}};

    }

    std::vector<std::vector<std::complex<float>>> get_four_qubit_clifford()
    {
        return {{0.35355339f, .0f, 0.35355339f*I, .0f, .0f, 0.35355339f, .0f, 0.35355339f*I, .0f, 0.35355339f*I, .0f, 0.35355339f, 0.35355339f*I, .0f, 0.35355339f, .0f},
                {.0f, 0.35355339f*I, .0f, 0.35355339f, -0.35355339f*I, .0f, -0.35355339f, .0f, -0.35355339f, .0f, -0.35355339f*I, .0f, .0f, 0.35355339f, .0f, 0.35355339f*I}, 
                {.0f, -0.35355339f*I, .0f, 0.35355339f, 0.35355339f*I, .0f, -0.35355339f, .0f, 0.35355339f, .0f, -0.35355339f*I, .0f, .0f, -0.35355339f, .0f, 0.35355339f*I}, 
                {-0.35355339f, .0f, 0.35355339f*I, .0f, .0f, -0.35355339f, .0f, 0.35355339f*I, .0f, -0.35355339f*I, .0f, 0.35355339f, -0.35355339f*I, .0f, 0.35355339f, .0f}, 
                {.0f, -0.35355339f, .0f, 0.35355339f*I, 0.35355339f, .0f, -0.35355339f*I, .0f, 0.35355339f*I, .0f, -0.35355339f, .0f, .0f, -0.35355339f*I, .0f, 0.35355339f}, 
                {-0.35355339f*I, .0f, 0.35355339f, .0f, .0f, -0.35355339f*I, .0f, 0.35355339f, .0f, -0.35355339f, .0f, 0.35355339f*I, -0.35355339f, .0f, 0.35355339f*I, .0f}, 
                {0.35355339f*I, .0f, 0.35355339f, .0f, .0f, 0.35355339f*I, .0f, 0.35355339f, .0f, 0.35355339f, .0f, 0.35355339f*I, 0.35355339f, .0f, 0.35355339f*I, .0f},     
                {.0f, 0.35355339f, .0f, 0.35355339f*I, -0.35355339f, .0f, -0.35355339f*I, .0f, -0.35355339f*I, .0f, -0.35355339f, .0f, .0f, 0.35355339f*I, .0f, 0.35355339f}, 
                {-0.35355339f, .0f, 0.35355339f*I, .0f, .0f, 0.35355339f, .0f, -0.35355339f*I, .0f, 0.35355339f*I, .0f, -0.35355339f, -0.35355339f*I, .0f, 0.35355339f, .0f}, 
                {.0f, -0.35355339f*I, .0f, 0.35355339f, -0.35355339f*I, .0f, 0.35355339f, .0f, -0.35355339f, .0f, 0.35355339f*I, .0f, .0f, -0.35355339f, .0f, 0.35355339f*I}, 
                {.0f, -0.35355339f*I, .0f, -0.35355339f, -0.35355339f*I, .0f, -0.35355339f, .0f, -0.35355339f, .0f, -0.35355339f*I, .0f, .0f, -0.35355339f, .0f, -0.35355339f*I},
                {-0.35355339f, .0f, -0.35355339f*I, .0f, .0f, 0.35355339f, .0f, 0.35355339f*I, .0f, 0.35355339f*I, .0f, 0.35355339f, -0.35355339f*I, .0f, -0.35355339f, .0f}, 
                {.0f, 0.35355339f, .0f, 0.35355339f*I, 0.35355339f, .0f, 0.35355339f*I, .0f, 0.35355339f*I, .0f, 0.35355339f, .0f, .0f, 0.35355339f*I, .0f, 0.35355339f},     
                {0.35355339f*I, .0f, 0.35355339f, .0f, .0f, -0.35355339f*I, .0f, -0.35355339f, .0f, -0.35355339f, .0f, -0.35355339f*I, 0.35355339f, .0f, 0.35355339f*I, .0f}, 
                {0.35355339f*I, .0f, -0.35355339f, .0f, .0f, -0.35355339f*I, .0f, 0.35355339f, .0f, -0.35355339f, .0f, 0.35355339f*I, 0.35355339f, .0f, -0.35355339f*I, .0f}, 
                {.0f, 0.35355339f, .0f, -0.35355339f*I, 0.35355339f, .0f, -0.35355339f*I, .0f, 0.35355339f*I, .0f, -0.35355339f, .0f, .0f, 0.35355339f*I, .0f, -0.35355339f}};

    }
    
    TEST_CASE("testing correct cliffords", "[matrix -> clifford]")
    {
        SECTION("1 qubit")
        {
            std::vector<std::vector<std::complex<float>>> hadamard_matrix = {{1/sqrt(2.0f), 1/sqrt(2.0f)} , {1/sqrt(2.0f), -1/sqrt(2.0f)}};

            REQUIRE(is_clifford_matrix(hadamard_matrix));
        }

        SECTION("3 qubits")
        {
            std::vector<std::vector<std::complex<float>>> clifford_matrix = get_three_qubit_clifford();

            REQUIRE(is_clifford_matrix(clifford_matrix));
        }

        SECTION("3 qubits, alternative")
        {
            std::vector<std::vector<std::complex<float>>> clifford_matrix = get_alternative_three_qubit_clifford();

            REQUIRE(is_clifford_matrix(clifford_matrix));
        }

        SECTION("3 qubits, global factor")
        {
            const std::complex<float> global_phase(1 / std::sqrt(2.0f), 1 / std::sqrt(2.0f));
            
            std::vector<std::vector<std::complex<float>>> clifford_matrix = get_three_qubit_clifford_with_phase(global_phase);

            REQUIRE(is_clifford_matrix(clifford_matrix));
        }

        SECTION("4 qubits")
        {
            std::vector<std::vector<std::complex<float>>> clifford_matrix = get_four_qubit_clifford();

            REQUIRE(is_clifford_matrix(clifford_matrix));
        }
        
    }

    void multiply_column_by_phase (std::vector<std::vector<std::complex<float>>> &matrix, const std::size_t col_index, const std::complex<float> phase)
    {
        for (std::size_t row_index = 0; row_index < matrix.size(); row_index++)
        {
            matrix[row_index][col_index] *= phase;
        }
    }

    void replace_column (std::vector<std::vector<std::complex<float>>> &matrix, const std::size_t col_index, const std::vector<std::complex<float>> replacement_col)
    {
        for (std::size_t row_index = 0; row_index < matrix.size(); row_index++)
        {
            matrix[row_index][col_index] = replacement_col[row_index];
        }
    }

    TEST_CASE("testing incorrect cliffords", "[matrix -> clifford]")
    {
        std::vector<std::vector<std::complex<float>>> matrix = get_three_qubit_clifford();

        SECTION("incorrect size")
        {
            std::vector<std::vector<std::complex<float>>> dim3_matrix (3, std::vector<std::complex<float>> (3));

            REQUIRE_FALSE(is_clifford_matrix(dim3_matrix));
        }

        SECTION("first column not stabiliser state")
        {
            matrix[0][0] = 0;

            REQUIRE_FALSE(is_clifford_matrix(matrix));
        }

        SECTION("Initial column all zero")
        {
            replace_column(matrix, 1, std::vector<std::complex<float>> (8, .0f));

            REQUIRE_FALSE(is_clifford_matrix(matrix));
        }

        SECTION("Remaining column all zero")
        {
            replace_column(matrix, 7, std::vector<std::complex<float>> (8, .0f));

            REQUIRE_FALSE(is_clifford_matrix(matrix));
        }

        SECTION("impossible P pauli phase")
        {
            multiply_column_by_phase(matrix, 1, std::complex<float> {1/sqrt(2.0f), 1/sqrt(2.0f)});
            matrix[0][1] *= std::complex<float> {1/sqrt(2.0f), -1/sqrt(2.0f)};

            REQUIRE_FALSE(is_clifford_matrix(matrix));
        }

        SECTION("inital column not stabilised")
        {
            matrix[0][1] *= -1.0f;

            REQUIRE_FALSE(is_clifford_matrix(matrix));
        }

        SECTION("remaining column not stabilised")
        {
            matrix[1][7] *= -1.0f;

            REQUIRE_FALSE(is_clifford_matrix(matrix));
        }

        SECTION("impossible initial W relative phase")
        {
            multiply_column_by_phase(matrix, 1, std::complex<float> {1/sqrt(2.0f), 1/sqrt(2.0f)});

            REQUIRE_FALSE(is_clifford_matrix(matrix));
        }

        SECTION("incorrect remaining relative phase")
        {
            multiply_column_by_phase(matrix, 7, -1.0f);

            REQUIRE_FALSE(is_clifford_matrix(matrix));
        }
    }

    TEST_CASE("returning clifford", "[matrix -> clifford]")
    {
        SECTION("clifford input, 3 qubit")
        {
            const std::vector<std::vector<std::complex<float>>> matrix = get_three_qubit_clifford_with_phase(std::complex<float> {1/sqrt(2.0f), 1/sqrt(2.0f)});

            Clifford clifford = clifford_from_matrix(matrix);

            REQUIRE_THAT(clifford.get_matrix(), RangeEquals(matrix));
        }

        SECTION("incorrect first column, 3 qubit")
        {
            std::vector<std::vector<std::complex<float>>> matrix = get_three_qubit_clifford_with_phase(std::complex<float> {1/sqrt(2.0f), 1/sqrt(2.0f)});
            matrix[6][0] *= -1.0f;

            REQUIRE_THROWS_AS(clifford_from_matrix(matrix), std::invalid_argument);
        }

        SECTION("incorrect remaining column, 3 qubit")
        {
            std::vector<std::vector<std::complex<float>>> matrix = get_three_qubit_clifford_with_phase(std::complex<float> {1/sqrt(2.0f), 1/sqrt(2.0f)});
            matrix[0][5] *= -1.0f;

            REQUIRE_THROWS_AS(clifford_from_matrix(matrix), std::invalid_argument);

            // Check doesn't throw
            clifford_from_matrix(matrix, true);
        }
    }
}
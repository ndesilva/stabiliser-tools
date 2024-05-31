#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include <iostream>

#include "check_matrix.h"
#include "stabiliser_state.h"
#include "pauli.h"
#include "f2_helper.h"

using namespace Catch::Matchers;
using namespace fst;

static constexpr std::complex<float> i = {0, 1};

namespace
{
    std::vector<std::complex<float>> matrix_vector_mult(const std::vector<std::vector<std::complex<float>>> &matrix, const std::vector<std::complex<float>> &vector)
    {
        std::size_t size = vector.size();

        std::vector<std::complex<float>> result (size, 0);

        for(std::size_t i = 0; i < size; i ++)
        {
            for (std::size_t j = 0; j < size; j++) 
            {
                result[i] += matrix[i][j] * vector[j];
            }
        }

        return result;
    }

    std::unordered_map<std::size_t, bool> get_quadratic_from_from_vector(std::size_t dimension, std::vector<std::size_t> non_zero_coeffs)
	{
		std::unordered_map<std::size_t, bool> quadratic_form;
		quadratic_form.reserve(dimension * (dimension + 1)/2);
		quadratic_form[0] = 0;

		for (std::size_t j = 0; j < dimension; j++)
		{
			for (std::size_t k = 0; k < dimension; k++)
			{
				quadratic_form[integral_pow_2(j) ^ integral_pow_2(k)] = 0;
			}
		}

		for (const auto coeff : non_zero_coeffs)
		{
			quadratic_form[coeff] = 1;
		}

		return quadratic_form;
	}

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

    TEST_CASE("check matrix from list of paulis", "[check_matrix]")
    {
        Pauli pauli1 = Pauli(2, 0b10, 0b01, 0, 0); // X Z
        Pauli pauli2 = Pauli(2, 0b00, 0b01, 0, 0); // 1 Z

        std::vector<Pauli> paulis = {pauli1, pauli2};

        Check_Matrix check_matrix(paulis);

        REQUIRE_FALSE(check_matrix.row_reduced);
        REQUIRE(*check_matrix.x_stabilisers[0] == pauli1);
        REQUIRE(*check_matrix.z_only_stabilisers[0] == pauli2);
    }

    TEST_CASE("check matrix from stabiliser state", "[check matrix]")
    {
        Stabiliser_State state (5, 3);
        state.basis_vectors = {0b10001, 0b10100, 0b11110};
        state.shift = 0b00001;

        state.real_linear_part = 0b00101;
        state.imaginary_part = 0b00101;
        state.quadratic_form = get_quadratic_from_from_vector(3, {0b110, 0b011});

        std::vector<std::complex<float>> old_statevector = state.get_state_vector();

        Check_Matrix check_matrix(state);
        std::vector<std::complex<float>> statevector = state.get_state_vector();

        REQUIRE(check_matrix.row_reduced);
        REQUIRE_THAT(old_statevector, RangeEquals(statevector));

        for (const auto pauli : check_matrix.paulis)
        {
            REQUIRE_THAT(matrix_vector_mult(pauli.get_matrix(), statevector), RangeEquals(statevector));
        }
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

        REQUIRE(starting_check_matrix.row_reduced);
        REQUIRE_THAT(starting_check_matrix.paulis, RangeEquals(row_reduced_check_matrix.paulis));

        for(std::size_t index = 0; index < row_reduced_check_matrix.x_stabilisers.size(); index++){
            REQUIRE(*row_reduced_check_matrix.x_stabilisers[index] == *starting_check_matrix.x_stabilisers[index]);
        }

        for(std::size_t index = 0; index < row_reduced_check_matrix.z_only_stabilisers.size(); index++){
            REQUIRE(*row_reduced_check_matrix.z_only_stabilisers[index] == *starting_check_matrix.z_only_stabilisers[index]);
        }
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
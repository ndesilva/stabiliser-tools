#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>

#include <unordered_set>

#include "stabiliser_state/check_matrix.h"
#include "stabiliser_state/stabiliser_state.h"
#include "pauli/pauli.h"
#include "test_util.h"

using namespace Catch::Matchers;
using namespace fst;
using namespace test;

namespace
{
    bool paulis_sorted_correctly(const Check_Matrix &check_matrix)
    {   
        if (check_matrix.number_qubits != check_matrix.paulis.size())
        {
            return false;
        }
        
        if (check_matrix.paulis.size() != check_matrix.x_stabilisers.size() + check_matrix.z_only_stabilisers.size())
        {
            return false;
        }

        std::unordered_set<Pauli, Pauli_Hasher> pauli_set (check_matrix.paulis.begin(), check_matrix.paulis.end());

        for (const auto pauli_pointer : check_matrix.x_stabilisers)
        {
            Pauli pauli = *pauli_pointer;
            if (!pauli_set.contains(pauli) || pauli.x_vector == 0)
            {
                return false;
            }
        }

        for (const auto pauli_pointer : check_matrix.z_only_stabilisers)
        {
            Pauli pauli = *pauli_pointer;
            if (!pauli_set.contains(pauli) || pauli.x_vector != 0)
            {
                return false;
            }
        }

        return true;
    }

    bool is_row_reduced(const Check_Matrix &check_matrix)
    {
        std::vector<std::size_t> x_vectors;
        std::vector<std::size_t> z_vectors;
        
        for(const auto pauli : check_matrix.paulis)
        {
            if(pauli.x_vector == 0)
            {
                z_vectors.push_back(pauli.z_vector);
            }
            else
            {
                x_vectors.push_back(pauli.x_vector);
            }
        }

        return vectors_row_reduced(x_vectors) && vectors_row_reduced(z_vectors);
    }

    std::vector<Pauli> get_pauli_list(){
        Pauli pauli1 = Pauli(2, 0b10, 0b01, 0, 0); // X Z
        Pauli pauli2 = Pauli(2, 0b00, 0b01, 0, 0); // 1 Z

        return {pauli1, pauli2};
    }

    TEST_CASE("check matrix from list of paulis", "[check_matrix]")
    {
        std::vector<Pauli> paulis = get_pauli_list();

        Check_Matrix check_matrix(paulis);

        REQUIRE_FALSE(check_matrix.row_reduced);
        REQUIRE(paulis_sorted_correctly(check_matrix));
        REQUIRE_THAT(check_matrix.paulis, RangeEquals(paulis));
    }

    Stabiliser_State get_stabiliser_state()
    {
        Stabiliser_State state (5, 3);
        state.basis_vectors = {0b10001, 0b10100, 0b11110};
        state.shift = 0b00001;

        state.real_linear_part = 0b00101;
        state.imaginary_part = 0b00101;
        state.quadratic_form = get_quadratic_from_from_vector(3, {0b110, 0b011});

        return state;
    }

    TEST_CASE("check matrix from stabiliser state", "[check matrix]")
    {
        Stabiliser_State state = get_stabiliser_state();
        std::vector<std::complex<float>> statevector = state.get_state_vector();

        Check_Matrix check_matrix(state);

        REQUIRE_THAT(statevector, RangeEquals(state.get_state_vector()));
        REQUIRE(check_matrix.row_reduced);
        REQUIRE(paulis_sorted_correctly(check_matrix));
        REQUIRE(is_row_reduced(check_matrix));

        for (const auto pauli : check_matrix.paulis)
        {
            REQUIRE_THAT(matrix_vector_mult(pauli.get_matrix(), statevector), RangeEquals(statevector));
        }
    }

    Check_Matrix get_unreduced_check_matrix()
    {
        std::vector<Pauli> paulis;

        paulis.push_back(Pauli(5, 0b00011, 0b00001, 0, 1));
        paulis.push_back(Pauli(5, 0b00101, 0b00010, 0, 0));
        paulis.push_back(Pauli(5, 0b00110, 0b00100, 1, 1));
        paulis.push_back(Pauli(5, 0b10000, 0b01000, 1, 0));
        paulis.push_back(Pauli(5, 0b10110, 0b00100, 0, 1));

        return Check_Matrix(paulis);
    }

    TEST_CASE("row reduce", "[check_matrix]")
    {
        Check_Matrix check_matrix = get_unreduced_check_matrix();
        std::vector<std::complex<float>> statevector = check_matrix.get_state_vector();

        check_matrix.row_reduce();

        REQUIRE(check_matrix.row_reduced);
        REQUIRE(paulis_sorted_correctly(check_matrix));
        REQUIRE(is_row_reduced(check_matrix));
        REQUIRE_THAT(check_matrix.get_state_vector(), RangeEquals(statevector));
    }

    std::vector<std::complex<float>> get_state_vector()
    {
        std::vector<std::complex<float>> state_vector (32, 0);
        const float root_8 = std::sqrt(8.0f);

        state_vector[0] = 1 / root_8;
        state_vector[3] = -I / root_8;
        state_vector[5] = 1 / root_8;
        state_vector[6] = I / root_8;
        state_vector[16] = -1 / root_8;
        state_vector[19] = I / root_8;
        state_vector[21] = -1 / root_8;
        state_vector[22] = -I / root_8;

        return state_vector;
    }

    TEST_CASE("get state vector") {
        Check_Matrix check_matrix = get_unreduced_check_matrix();
        std::vector<std::complex<float>> expected_state_vector = get_state_vector();

        REQUIRE(check_matrix.get_state_vector() == expected_state_vector);
    }
}
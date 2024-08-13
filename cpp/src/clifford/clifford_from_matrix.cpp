#include "clifford_from_matrix.h"

#include "util/f2_helper.h"
#include "stabiliser_state/stabiliser_state.h"
#include "stabiliser_state/check_matrix.h"
#include "stabiliser_state/stabiliser_state_from_statevector.h"

#include <optional>
#include <tuple>

using namespace fst;

namespace
{
    std::vector<std::vector<std::complex<float>>> transpose_matrix(const std::vector<std::vector<std::complex<float>>> &matrix, const std::size_t &size)
    {
        std::vector<std::vector<std::complex<float>>> transposed_matrix ( size, std::vector<std::complex<float>> (size) );
        
        for (std::size_t i = 0; i < size; i++)
        {
            for (std::size_t j = 0; j < size; j++)
            {
                transposed_matrix[i][j] = matrix[j][i];
            }
        }

        return transposed_matrix;
    }

    template <bool assume_valid, bool return_state>
    auto clifford_from_matrix_internal(const std::vector<std::vector<std::complex<float>>> &matrix)
        -> std::conditional_t<return_state, std::optional<fst::Clifford>, bool>
    {
        const std::size_t size = matrix.size();
        std::vector<std::vector<std::complex<float>>> transposed_matrix = transpose_matrix(matrix, size);

        if (!is_power_of_2(size))
        {   
            return {};
        }

        Stabiliser_State first_col_state;

        try
        {
            first_col_state = std::move(stabiliser_from_statevector(transposed_matrix[0], assume_valid));
        }
        catch (...)
        {   
            return {};
        }

        std::size_t number_qubits = first_col_state.number_qubits;

        Check_Matrix first_col_check_matrix (first_col_state);
        
        std::vector<Pauli> uncorrected_W_paulis;
        uncorrected_W_paulis.reserve(number_qubits);

        // TODO: we calculate these when constructing the check_matrix, remove calculation?
        std::unordered_set<std::size_t> pivot_indices;
        pivot_indices.reserve(number_qubits);

        for (auto pauli_p : first_col_check_matrix.x_stabilisers)
        {
            std::size_t pivot_index = integral_log_2(pauli_p->x_vector);
            pivot_indices.insert(pivot_index);            
            uncorrected_W_paulis.push_back(Pauli(number_qubits, 0, integral_pow_2(pivot_index), 0, 0));
        }

        for (std::size_t i = 0; i < number_qubits; i++)
        {
            if (!pivot_indices.contains(i))
            {
                uncorrected_W_paulis.push_back(Pauli(number_qubits, integral_pow_2(i),0, 0, 0));
            }
        }

        std::vector<std::size_t> first_col_effects (number_qubits, 0);

        for (std::size_t i = 0; i < number_qubits; i++)
        {
            std::size_t col_index = integral_pow_2(i);
            std::size_t row_index = 0;

            while (row_index < size && transposed_matrix[col_index][row_index] == .0f)
            {
                ++row_index;
            }
            
            if (row_index == size)
            {
                return {};
            }

            std::complex<float> non_zero_entry = transposed_matrix[col_index][row_index];

            for (std::size_t j = 0; j < number_qubits; j++)
            {
                //TODO make pointer?
                Pauli pauli = first_col_check_matrix.paulis[j];
                std::complex<float> phase = transposed_matrix[col_index][row_index ^ pauli.x_vector]/(non_zero_entry * sign_f2_dot_product(row_index, pauli.z_vector) * pauli.get_phase());

                if (std::norm(phase + 1.0f) < 0.125)
                {
                    first_col_effects[j] ^= col_index;
                }
                else if (std::norm(phase - 1.0f) >= 0.125)
                {
                    return {};
                }
            }
        }

        std::vector<std::size_t> pauli_ordering (number_qubits);

        for (std::size_t i = 0; i < number_qubits; i++)
        {
            if (first_col_effects[i] == 0)
            {
                return {};
            }

            std::size_t pivot_index = integral_log_2(first_col_effects[i]);
            pauli_ordering[pivot_index] = i;

            for (std::size_t j = 0; j < number_qubits; j++)
            {
                if ( i != j && bit_set_at(first_col_effects[j], pivot_index))
                {
                    first_col_effects[j] ^= first_col_effects[i];
                    first_col_check_matrix.paulis[j].multiply_by_pauli_on_right(first_col_check_matrix.paulis[i]);
                    uncorrected_W_paulis[i].multiply_by_pauli_on_right(uncorrected_W_paulis[j]);
                }
            }
        }

        std::vector<Pauli> z_conjugates (number_qubits);
        std::vector<Pauli> W_paulis (number_qubits);

        for (std::size_t i = 0; i < number_qubits; i++)
        {
            z_conjugates[i] = std::move(first_col_check_matrix.paulis[pauli_ordering[i]]);
            W_paulis[i] = std::move(uncorrected_W_paulis[pauli_ordering[i]]);
        }

        if constexpr (!assume_valid)
        {
            for (std::size_t col_index = 1; col_index < size; col_index++)
            {
                for (std::size_t i = 1; i < number_qubits; i++)
                {
                    if (! z_conjugates[i].has_eigenstate(transposed_matrix[col_index], bit_set_at(col_index, i)))
                    {
                        return {};
                    }
                }
            }
        }

        for (std::size_t i = 0; i < number_qubits; i++)
        {
            std::size_t non_zero_index = first_col_state.shift ^ W_paulis[i].x_vector;
            std::size_t col_index = integral_pow_2(i);
            std::complex<float> relative_phase = transposed_matrix[col_index][non_zero_index]/(transposed_matrix[0][first_col_state.shift]*sign_f2_dot_product(first_col_state.shift, W_paulis[i].z_vector)*W_paulis[i].get_phase());

            if (std::norm(relative_phase + 1.0f) < 0.125)
            {
                W_paulis[i].sign_bit ^= 1;
            }
            else if (std::norm(relative_phase + std::complex<float>{0,1}) < 0.125)
            {
                W_paulis[i].sign_bit ^= W_paulis[i].imag_bit;
                W_paulis[i].imag_bit ^= 1;
            }
            else if (std::norm(relative_phase - std::complex<float>{0,1}) < 0.125)
            {
                W_paulis[i].sign_bit ^= !W_paulis[i].imag_bit;
                W_paulis[i].imag_bit ^= 1;
            }
            else if (std::norm(relative_phase - 1.0f) >= 0.125)
            {
                return {};
            }
        }

        for (std::size_t i = 0; i < number_qubits; i++)
        {
            std::size_t i_non_zero_index = first_col_state.shift ^ W_paulis[i].x_vector;
            std::size_t i_col_index = integral_pow_2(i);
            std::complex<float> i_non_zero_entry = transposed_matrix[i_col_index][i_non_zero_index]; 

            for (std::size_t j = 0; j < number_qubits; j++)
            {
                std::size_t ij_non_zero_index = i_non_zero_index ^ W_paulis[j].x_vector;
                std::complex<float> relative_phase = transposed_matrix[i_col_index ^ integral_pow_2(j)][ij_non_zero_index]/(i_non_zero_entry*sign_f2_dot_product(i_non_zero_index, W_paulis[j].z_vector)*W_paulis[j].get_phase());

                if (std::norm(relative_phase + 1.0f) < 0.125)
                {
                    W_paulis[j].multiply_by_pauli_on_right(z_conjugates[i]);
                }
                else if (std::norm(relative_phase - 1.0f) >= 0.125)
                {
                    return {};
                }
            }
        }

        if constexpr (!assume_valid)
        {
            std::size_t old_col_index = 0;
            std::size_t old_support = first_col_state.shift;

            for (std::size_t i = 1; i < size ; i++)
            {   
                // Iterate through the gray code
                std::size_t new_col_index = i ^ (i >> 1);
                std::size_t flipped_bit = integral_log_2(new_col_index ^ old_col_index);

                Pauli pauli_flip = W_paulis[flipped_bit];
                std::size_t new_support = old_support ^ pauli_flip.x_vector;

                if (std::norm(transposed_matrix[new_col_index][new_support] - transposed_matrix[old_col_index][old_support]*sign_f2_dot_product(old_support, pauli_flip.z_vector)*pauli_flip.get_phase()) >= 0.001)
                {
                    return {};
                }
                
                old_col_index = new_col_index;
                old_support = new_support;
            }
        }

        if constexpr (return_state)
        {
            return Clifford (z_conjugates, W_paulis, first_col_state.global_phase);
        }
        else
        {
            return true;
        }
    }
}

fst::Clifford fst::clifford_from_matrix(const std::vector<std::vector<std::complex<float>>> &matrix, const bool assume_valid)
{
    std::optional<Clifford> clifford = assume_valid 
                                ? clifford_from_matrix_internal<true, true>(matrix)
                                : clifford_from_matrix_internal<false, true>(matrix);

    if (!clifford)
    {
        throw std::invalid_argument("Matrix was not a Clifford");
    }

    return *std::move(clifford);
}

bool fst::is_clifford_matrix(const std::vector<std::vector<std::complex<float>>> &matrix )
{
    return clifford_from_matrix_internal<false, false>(matrix);
}
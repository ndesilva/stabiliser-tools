#include "check_matrix.h"
#include "stabiliser_state.h"
#include "f2_helper.h"

namespace fst
{
    Check_Matrix::Check_Matrix(std::vector<Pauli> paulis, bool row_reduced)
        : paulis(paulis), row_reduced(row_reduced)
    {
        number_qubits = paulis.size();
        categorise_paulis();
    }

    void Check_Matrix::categorise_paulis()
    {
        for (auto &pauli : paulis)
        {
            if (pauli.x_vector == 0)
            {
                z_only_stabilisers.push_back(&pauli);
            }
            else{
                x_stabilisers.push_back(&pauli);
            }
        }
    }

    Check_Matrix::Check_Matrix(Stabiliser_State &stabiliser_state)
    {
        number_qubits = stabiliser_state.number_qubits;

        stabiliser_state.row_reduce_basis();

        paulis.reserve(number_qubits);

        std::vector<std::size_t> pivot_vectors;
        std::unordered_set<int> pivot_indices_set;

        for(std::size_t i = 0; i < stabiliser_state.dim; i++)
        {
            int pivot_index = integral_log_2(stabiliser_state.basis_vectors[i]);
            pivot_indices_set.insert(pivot_index);
            pivot_vectors[i] = integral_pow_2((std::size_t) pivot_index);
        }

        add_z_only_stabilisers(pivot_vectors, pivot_indices_set, stabiliser_state);
        add_x_stabilisers(pivot_vectors, stabiliser_state);

        row_reduced = true;
    }

    void Check_Matrix::add_z_only_stabilisers(std::vector<std::size_t> &pivot_vectors, std::unordered_set<int> pivot_indices_set, Stabiliser_State &state)
    {
        for(std::size_t i = 0; i < number_qubits; i++)
        {
            if (!pivot_indices_set.contains(i))
            {
                std::size_t alpha = integral_pow_2(i);

                // Make alpha perpendicular to the basis vectors
                for (std::size_t j = 0; j < state.dim; j++)
                {
                    alpha |= bit_set_at(state.basis_vectors[j], i) * pivot_vectors[j];
                }

            bool sign_bit = f2_dot_product(alpha, state.shift);
            
            Pauli pauli(number_qubits, 0, alpha, sign_bit, 0);
            paulis.push_back(pauli);
            z_only_stabilisers.push_back(&pauli);
            }

        }
    }

    void Check_Matrix::add_x_stabilisers(std::vector<std::size_t> &pivot_vectors, Stabiliser_State &state)
    {
        for (std::size_t i = 0; i < state.dim; i++)
        {
            bool imag_bit = bit_set_at(state.imaginary_part, i);
            
            std::size_t z_vector = 0;

            // TODO explain why this works
            for (std::size_t j = 0; j < state.dim; j++)
            {
                z_vector ^= pivot_vectors[j] * (state.quadratic_form[integral_pow_2(i) ^ integral_pow_2(j)] ^ (imag_bit * bit_set_at(state.imaginary_part, j) ));
            }

            bool sign_bit = bit_set_at(state.real_linear_part, i) ^ imag_bit ^ f2_dot_product(z_vector, state.shift);

            Pauli pauli(number_qubits, state.basis_vectors[i], z_vector, sign_bit, imag_bit);
            paulis.push_back(pauli);
            x_stabilisers.push_back(&pauli);
        }
    }

    std::vector<std::complex<float>> Check_Matrix::get_state_vector()
    {
        return Stabiliser_State(*this).get_state_vector();
    }

    void Check_Matrix::row_reduce()
    {
        if (row_reduced) {return;}

        row_reduce_x_stabilisers();
        row_reduce_z_only_stabilisers();

        row_reduced = true;
    }

    void Check_Matrix::row_reduce_x_stabilisers()
    {
        for (std::size_t i = 0; i < x_stabilisers.size(); i++)
        {
            Pauli *pauli = x_stabilisers[i];
            int pivot_index = integral_log_2((*pauli).x_vector);

            if (pivot_index == -1)
            {
                x_stabilisers.erase(x_stabilisers.begin() + i);
                z_only_stabilisers.push_back(pauli);
                i--;
            }
            else
            {
                auto u_pivot_index = (std::size_t) pivot_index;
                for (auto &other_pauli : x_stabilisers)
                {
                    if (other_pauli != pauli && bit_set_at((*other_pauli).x_vector, u_pivot_index))
                    {
                        (*other_pauli).multiply_by_pauli_on_right(*pauli);
                    }
                }
            }
        }
    }

    // TODO this repeats alot of code from reducing x_stabilisers, optimise?
    void Check_Matrix::row_reduce_z_only_stabilisers()
    {
        for (std::size_t i = 0; i < z_only_stabilisers.size(); i++)
        {
            Pauli *pauli = z_only_stabilisers[i];
            std::size_t pivot_index = integral_log_2((*pauli).z_vector);

            for (auto &other_pauli : z_only_stabilisers)
            {
                if (other_pauli != pauli && bit_set_at((*other_pauli).z_vector, pivot_index))
                {
                    (*other_pauli).multiply_by_pauli_on_right(*pauli);
                }
            }
        }
    }
}
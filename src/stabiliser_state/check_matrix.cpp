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

        std::unordered_map<std::size_t, bool> m_quadratic_form = stabiliser_state.get_quadratic_form_as_map();
        stabiliser_state.row_reduce_basis(m_quadratic_form);

        paulis.reserve(number_qubits);

        std::set<int> pivot_indices;

        for(const auto &vector : stabiliser_state.basis_vectors)
        {
            pivot_indices.insert(integral_log_2(vector));
        }

        add_z_only_stabilisers(pivot_indices);
        add_x_stabilisers(pivot_indices, m_quadratic_form);

        row_reduced = true;
    }

    void Check_Matrix::add_z_only_stabilisers(std::set<int> &pivot_indices)
    {
        for(std::size_t i = 0; i < number_qubits; i++)
        {
            if (!pivot_indices.contains(i))
            {
                std::size_t alpha = integral_pow_2(i);

                // make alpha perpendicular to the basis vectors
                for (std::size_t j = 0; j < pivot_indices.size(); j++)
                {
                    alpha |= integral_pow_2(j)*(basis_vectors[j], i);
                }
            }
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
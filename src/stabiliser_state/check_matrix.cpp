#include "check_matrix.h"
#include "f2_helper.h"

namespace fst
{
    Check_Matrix::Check_Matrix(std::vector<Pauli> paulis)
        : paulis(paulis)
    {
        number_qubits = paulis.size();
        categorise_paulis();
    }

    void Check_Matrix::categorise_paulis()
    {
        for (auto &pauli : paulis)
        {
            if (integral_log_2(pauli.x_vector) == -1)
            {
                z_only_stabilisers.push_back(&pauli);
            }
            else{
                x_stabilisers.push_back(&pauli);
            }
        }
    }

    Stabiliser_State Check_Matrix::get_stabiliser_state()
    {
        row_reduce();
        
        std::size_t dim = x_stabilisers.size();

        Stabiliser_State state(number_qubits, dim);

        set_support(state);
        set_linear_and_quadratic_forms(state);
        state.row_reduced = true;

        return state;
    }

    std::vector<std::complex<float>> Check_Matrix::get_state_vector()
    {
        return get_stabiliser_state().get_state_vector();
    }

    void Check_Matrix::set_support(Stabiliser_State &state) const
    {
        set_basis_vectors(state);
        set_shift(state);
    }

    void Check_Matrix::set_linear_and_quadratic_forms(Stabiliser_State &state) const
    {
        std::size_t imaginary_part = 0;
        std::size_t real_linear_part = 0;
        std::vector<std::size_t> quadratic_form;

        for (std::size_t j = 0; j < state.dim; j++)
        {
            std::size_t v_j = state.basis_vectors[j];
            Pauli *p_j = x_stabilisers[j];
            std::size_t beta_j = (*p_j).z_vector;
            std::size_t imag_bit = (*p_j).imag_bit;

            imaginary_part |= integral_pow_2(j) * imag_bit;
            real_linear_part |= integral_pow_2(j) * ((*p_j).sign_bit ^ f2_dot_product(beta_j, v_j ^ state.shift));

            for (std::size_t i = 0; i < j; i++)
            {
                std::size_t v_i = state.basis_vectors[i];
                std::size_t other_imag_bit = (*x_stabilisers[i]).imag_bit; // TODO we are accessing the imag_bits alot, optimise?

                if (f2_dot_product(beta_j, v_i) ^ imag_bit*other_imag_bit)
                {
                    quadratic_form.push_back(integral_pow_2(i) | integral_pow_2(j));
                }
            }
        }

        state.imaginary_part = imaginary_part;
        state.real_linear_part = real_linear_part;
        state.quadratic_form = quadratic_form;
    }

    void Check_Matrix::set_basis_vectors(fst::Stabiliser_State &state) const
    {
        std::vector<std::size_t> basis_vectors;
        basis_vectors.reserve(state.dim);

        for (const auto &pauli : x_stabilisers)
        {
            basis_vectors.push_back((*pauli).x_vector);
        }

        state.basis_vectors = std::move(basis_vectors);
    }

    void Check_Matrix::set_shift(fst::Stabiliser_State &state) const
    {
        std::size_t shift = 0;

        for (const auto &pauli : z_only_stabilisers)
        {
            std::size_t pivot_index = integral_log_2((*pauli).z_vector);
            shift |= (integral_pow_2(pivot_index) * (*pauli).sign_bit);
        }

        state.shift = shift;
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
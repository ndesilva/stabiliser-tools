#ifndef _FAST_STABILISER_CHECK_MATRIX_H
#define _FAST_STABILISER_CHECK_MATRIX_H

#include "pauli.h"
#include "stabiliser_state.h"

#include <vector>
#include <complex>

namespace fst
{
    /// TODO add documentation
    struct Check_Matrix
    {
        std::size_t number_qubits = 0;
        bool row_reduced = false;
        
        std::vector<Pauli> paulis;
        std::vector<Pauli *> z_only_stabilisers;
        std::vector<Pauli *> x_stabilisers;

        Check_Matrix(std::vector<Pauli> paulis);

        void categorise_paulis();
        
        Stabiliser_State get_stabiliser_state();
        std::vector<std::complex<float>> get_state_vector();

        void row_reduce();
        void row_reduce_x_stabilisers();
        void row_reduce_z_only_stabilisers();

        void set_support(Stabiliser_State &state) const;
        void set_basis_vectors(fst::Stabiliser_State &state) const;
        void set_shift(fst::Stabiliser_State &state) const;
        void set_linear_and_quadratic_forms(Stabiliser_State &state) const;
    };
}

#endif
#ifndef _FAST_STABILISER_CHECK_MATRIX_H
#define _FAST_STABILISER_CHECK_MATRIX_H

#include "pauli.h"

#include <vector>
#include <complex>
#include <set>

namespace fst
{
    struct Stabiliser_State;

    /// TODO add documentation
    struct Check_Matrix
    {
        std::size_t number_qubits = 0;
        
        std::vector<Pauli> paulis;
        std::vector<Pauli *> z_only_stabilisers;
        std::vector<Pauli *> x_stabilisers;
        
        bool row_reduced;

        explicit Check_Matrix(std::vector<Pauli> paulis, bool row_reduced = false);
        explicit Check_Matrix(Stabiliser_State &stabiliser_state);

        void categorise_paulis();
        
        void add_z_only_stabilisers(std::set<int> &pivot_indices);
		void add_x_stabilisers(std::set<int> &pivot_indices, std::map<std::size_t, bool> &m_quadratic_form);
        
        std::vector<std::complex<float>> get_state_vector();

        void row_reduce();
        void row_reduce_x_stabilisers();
        void row_reduce_z_only_stabilisers();
    };
}

#endif
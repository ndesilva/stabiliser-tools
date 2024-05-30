#ifndef _FAST_STABILISER_CHECK_MATRIX_H
#define _FAST_STABILISER_CHECK_MATRIX_H

#include "pauli.h"

#include <vector>
#include <complex>
#include <unordered_set>

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

        std::vector<std::complex<float>> get_state_vector();
        
        void row_reduce();

        private:
        
        void categorise_paulis();
        
        void add_z_only_stabilisers(std::vector<std::size_t> &pivot_vectors, std::unordered_set<std::size_t> pivot_indices_set, Stabiliser_State &state);
		void add_x_stabilisers(std::vector<std::size_t> &pivot_vectors, Stabiliser_State &state);  

        void row_reduce_x_stabilisers();
        void row_reduce_z_only_stabilisers();
    };
}

#endif
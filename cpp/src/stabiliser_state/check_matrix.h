#ifndef _FAST_STABILISER_CHECK_MATRIX_H
#define _FAST_STABILISER_CHECK_MATRIX_H

#include "pauli/pauli.h"

#include <vector>
#include <complex>
#include <unordered_set>

namespace fst
{
    struct Stabiliser_State;

    /// The class used to represent a list of n commuting paulis, an alternative representation
    /// of a stabiliser state
    struct Check_Matrix
    {
        std::size_t number_qubits = 0;
        
        std::vector<Pauli> paulis;
        
        /// Paulis are sorted into 2 types: "z_only", which have no X component, and "x_stabilisers",
        /// which may have both an x and z component
        std::vector<Pauli *> z_only_stabilisers;
        std::vector<Pauli *> x_stabilisers;
        
        bool row_reduced;

        explicit Check_Matrix(const std::vector<Pauli> &paulis, const bool row_reduced = false);
        explicit Check_Matrix(Stabiliser_State &stabiliser_state);

        /// Return the state vector of length 2^n stabilised by each of the Paulis in the check matrix
        std::vector<std::complex<float>> get_state_vector();
        
        /// Row reduce the check_matrix, giving a new set of paulis that generate the same stabiliser group.
        /// The new paulis have the x_vectors of the "x_stabiliser" paulis, and z_vectors of the "z_only" stabilisers
        /// in reduced row_echelon form. Note that the collection of all the pauli's z_vectors may NOT be in reduced row
        /// echelon form.
        void row_reduce();

        private:
        
        void categorise_paulis();
        
        void add_z_only_stabilisers(const std::vector<std::size_t> &pivot_vectors, const std::unordered_set<std::size_t> &pivot_indices_set, const Stabiliser_State &state);
		void add_x_stabilisers(const std::vector<std::size_t> &pivot_vectors, const Stabiliser_State &state);  

        void row_reduce_x_stabilisers();
        void row_reduce_z_only_stabilisers();
    };
}

#endif
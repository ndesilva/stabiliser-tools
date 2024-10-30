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
        
        // Get the list of Stabilisers
        const std::vector<Pauli>& get_paulis() const;
        // Set the list of Stabilisers
        // TODO : use std::forward to reduce overhead?
        void set_paulis(std::vector<Pauli> paulis_);
        
        /// Paulis are sorted into 2 types: "z_only", which have no X component, and "x_stabilisers",
        /// which may have both an x and z component
        const std::vector<Pauli *>& get_z_only_stabilisers() const;
        const std::vector<Pauli *>& get_x_stabilisers() const;

        /// IF THE CHECK MATRIX IS ROW REDUCED, then this returns a list of the pivot columns of the "z_only"
        /// stabilisers (correspdonding to the order of the z_only_stabiliser list). The pivot column of a "z_only"
        /// stabiliser should contain a 1, where all other "z_only" stabiliers have a zero there. Moreover, it should *not*
        /// be a pivot_column for the "x_stabilisers".
        const std::vector<size_t> & get_z_only_pivots() const;
        
        bool row_reduced;

        explicit Check_Matrix(const std::vector<Pauli> paulis, const bool row_reduced = false);
        explicit Check_Matrix(Stabiliser_State &stabiliser_state);

        /// Return the state vector of length 2^n stabilised by each of the Paulis in the check matrix
        std::vector<std::complex<float>> get_state_vector();
        
        /// Row reduce the check_matrix, giving a new set of paulis that generate the same stabiliser group.
        /// The new paulis have the x_vectors of the "x_stabiliser" paulis, and z_vectors of the "z_only" stabilisers
        /// in reduced row echelon form. Note that the collection of all the paulis' z_vectors may NOT be in reduced row
        /// echelon form.
        void row_reduce();

        private:

        std::vector<Pauli> paulis;
        std::vector<Pauli *> z_only_stabilisers;
        std::vector<Pauli *> x_stabilisers;
        
        std::vector<size_t> z_only_pivots;
                
        void categorise_paulis();
        
        void add_z_only_stabilisers(const std::vector<std::size_t> &pivot_vectors, const std::unordered_set<std::size_t> &pivot_indices_set, const Stabiliser_State &state);
		void add_x_stabilisers(const std::vector<std::size_t> &pivot_vectors, const Stabiliser_State &state);  

        void row_reduce_x_stabilisers();
        void row_reduce_z_only_stabilisers();

        void set_z_only_pivots();
    };
}

#endif
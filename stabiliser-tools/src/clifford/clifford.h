#ifndef _FAST_STABILISER_CLIFFORD_H
#define _FAST_STABILISER_CLIFFORD_H

#include "pauli/pauli.h"

#include <vector>
#include <complex>

namespace fst
{
    /// The class used to represent a Clifford operator U.
    /// Represented by its action on the Pauli basis:
    /// z_conjugates[i] = UZ_iU*, x_conjugates[i] = UX_iU*
    struct Clifford
    {
        std::size_t number_qubits = 0;

        std::vector<Pauli> z_conjugates;
        std::vector<Pauli> x_conjugates;

        std::complex<float> global_phase;

        Clifford(const std::vector<Pauli> z_conjugates, const std::vector<Pauli> x_conjugates, const std::complex<float> global_phase = 1.0f);

        /// Returns the matrix of the Clifford (with respect to the computational basis) 
        std::vector<std::vector<std::complex<float>>> get_matrix() const; 
    };
}

#endif
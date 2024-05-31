#ifndef _FAST_STABILISER_TEST_UTIL_H
#define _FAST_STABILISER_TEST_UTIL_H

#include <vector>
#include <complex>
#include <unordered_map>

#include "pauli.h"

using namespace fst;

namespace test {
    
    static constexpr std::complex<float> I = {0, 1};

    /// Returns the matrix-vector product of the given matrix and vector
    std::vector<std::complex<float>> matrix_vector_mult(const std::vector<std::vector<std::complex<float>>> &matrix, const std::vector<std::complex<float>> &vector);

    /// Given a vector of non-zero coefficient of the quadratic form, create the correspdonding unordered_map (quadratic_form[coefficient] = 1 if present, 0 if not)
    std::unordered_map<std::size_t, bool> get_quadratic_from_from_vector(const std::size_t dimension, const std::vector<std::size_t> &non_zero_coeffs);

    /// Checks whether the given vector of F2 vectors (represented as ints) is row reduced
    bool vectors_row_reduced(std::vector<std::size_t> &vectors);

    struct Pauli_Hasher
    {
        std::size_t operator() (const Pauli & pauli) const;
    };
}

#endif
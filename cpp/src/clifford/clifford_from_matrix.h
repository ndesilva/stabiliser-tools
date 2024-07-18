#ifndef _FAST_STABILISER_CLIFFORD_FROM_MATRIX_H
#define _FAST_STABILISER_CLIFFORD_FROM_MATRIX_H

#include <complex>

#include "clifford.h"

namespace fst
{
    /// Convert a 2^n by 2^n matrix with complex entries into a clifford object.
	///
	/// Assuming valid is faster, but will result in undefined behaviour if the matrix is not in fact a
	/// valid clifford operator
    Clifford clifford_from_matrix (const std::vector<std::vector<std::complex<float>>> &matrix, const bool assume_valid = false);

    /// Test wheter a matrix with complex entries corresponds to a clifford state.
    bool is_clifford_matrix(const std::vector<std::vector<std::complex<float>>> &matrix);
}

#endif
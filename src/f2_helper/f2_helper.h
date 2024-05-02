#ifndef _FAST_STABILISER_F2_HELPER_H
#define _FAST_STABILISER_F2_HELPER_H

#include <complex>
#include <vector>

namespace fst {

/// Gives the F_2 inner product between 2 F_2 vectors (represented
/// as integers)
int f2_dot_product(const size_t x, const size_t y);

/// Gives (-1)^(x.y), where . is the F_2 inner product between 2
/// F_2 vectors (represented as integers)
int sign_f2_dot_product(const size_t x, const size_t y);

/// Gives (i)^(x.y), where . is the F_2 inner product between 2
/// F_2 vectors (represented as integers)
std::complex<float> imag_f2_dot_product(const size_t x, const size_t y);

/// Given vector_index, the column vector of an element of the vector
/// space (represented as an integer) with respect to the vector basis, and
/// Q a quadratic form with respect to the same basis (represented as a list of 
/// coefficients, i.e. 101 corresponds to x_0 x_1), find the value of (-1)^Q(vector_index).
int evaluate_quadratic_form(const size_t vector_index, const std::vector<size_t> &quadratic_form);

/// Return the number of binary digits of number, i.e. the ceiling of log_2(number)
int integral_log_2(const size_t number);

}

#endif
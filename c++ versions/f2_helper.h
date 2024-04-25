#ifndef _FAST_STABILISER_F2_HELPER_H
#define _FAST_STABILISER_F2_HELPER_H

#include <complex>

namespace fst {

/// Gives the F_2 inner product between 2 F_2 vectors (represented
/// as integers)
int f2_dot_product(int x, int y);

/// Gives (-1)^(x.y), where . is the F_2 inner product between 2
/// F_2 vectors (represented as integers)
int sign_f2_dot_product(int x, int y);

/// Gives (i)^(x.y), where . is the F_2 inner product between 2
/// F_2 vectors (represented as integers)
std::complex<int> imag_f2_dot_product(int x, int y);

}

#endif
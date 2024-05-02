#include "f2_helper.h"

#include <bit>

using namespace fst;

int fst::f2_dot_product(const size_t x, const size_t y) {
    const auto product =  x & y;
    const int hamming_weight = std::popcount(product);
    return hamming_weight % 2;
}

int fst::sign_f2_dot_product(const size_t x, const size_t y) {
    return 1 - 2*f2_dot_product(x,y);
}

std::complex<float> fst::imag_f2_dot_product(const size_t x, const size_t y) {
    int dot_product = f2_dot_product(x,y);
    return std::complex<float> (1-dot_product, dot_product);
}

int fst::evaluate_quadratic_form(const size_t vector_index, const std::vector<size_t> &quadratic_form) {
    int mod2_result = 0;
    
    for(const auto &term : quadratic_form) {
        mod2_result ^= ((term & vector_index) == term);
    }

    return 1 - 2*mod2_result;
}

// TODO : test and benchmark? fastest way 
int fst::integral_log_2(const size_t number){
    return std::bit_width(number) - 1;
}
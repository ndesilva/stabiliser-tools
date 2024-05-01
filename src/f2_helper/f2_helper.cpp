#include "f2_helper.h"

using namespace fst;

int fst::f2_dot_product(int x, int y) {
    int product = x & y;
    int weight = 0;

    while(product) {
        product &= product - 1;
        weight ^= 1;
    } 

    return weight;
}

int fst::sign_f2_dot_product(int x, int y){
    return 1 - 2*f2_dot_product(x,y);
}

std::complex<float> fst::imag_f2_dot_product(int x, int y) {
    int dot_product = f2_dot_product(x,y);
    return std::complex<float> (1-dot_product, dot_product);
}

int fst::evaluate_quadratic_form(int vector_index, const std::vector<int> &quadratic_form) {
    int mod2_result = 0;
    
    for(const auto &term : quadratic_form) {
        mod2_result ^= ((term & vector_index) == term);
    }

    return 1 - 2*mod2_result;
}

int fst::integral_log_2(int number){
    return ceil(log2(number));
}
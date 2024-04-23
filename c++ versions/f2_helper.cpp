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

std::complex<int> fst::imag_f2_dot_product(int x, int y){
    int dot_product = f2_dot_product(x,y);
    return std::complex<int> (dot_product, 1-dot_product);
}

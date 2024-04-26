#include <catch2/catch_test_macros.hpp>

#include "f2_helper.h"

using namespace fst;

TEST_CASE("mod 2 product","[f2 helper]"){
    int x = 26; // 11010
    int y = 19; // 10011
    int z = 3;  // 00011

    REQUIRE(f2_dot_product(x,y) == 0);
    REQUIRE(f2_dot_product(x,z) == 1);
}

TEST_CASE("signed mod 2 product", "[f2 helper]") {
    int x = 51; // 110011
    int y = 3;  // 000011
    int z = 49; // 110101

    REQUIRE(sign_f2_dot_product(x,y) == 1);
    REQUIRE(sign_f2_dot_product(x,z) == -1);
}

TEST_CASE("imaginary mod 2 product", "[f2 helper]") {
    int x = 43; // 101011
    int y = 3;  // 000011
    int z = 45; // 101101

    REQUIRE(imag_f2_dot_product(x,y) == float(1));
    REQUIRE(imag_f2_dot_product(x,z) == std::complex<float> (0, 1));
}
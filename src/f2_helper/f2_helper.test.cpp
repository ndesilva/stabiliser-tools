#include <catch2/catch_test_macros.hpp>

#include "f2_helper.h"

using namespace fst;

TEST_CASE("mod 2 product","[f2 helper]"){
    size_t x = 26; // 11010
    size_t y = 19; // 10011
    size_t z = 3;  // 00011

    REQUIRE(f2_dot_product(x,y) == 0);
    REQUIRE(f2_dot_product(x,z) == 1);
}

TEST_CASE("signed mod 2 product", "[f2 helper]") {
    size_t x = 51; // 110011
    size_t y = 3;  // 000011
    size_t z = 49; // 110101

    REQUIRE(sign_f2_dot_product(x,y) == 1);
    REQUIRE(sign_f2_dot_product(x,z) == -1);
}

TEST_CASE("imaginary mod 2 product", "[f2 helper]") {
    size_t x = 43; // 101011
    size_t y = 3;  // 000011
    size_t z = 45; // 101101

    REQUIRE(imag_f2_dot_product(x,y) == float(1));
    REQUIRE(imag_f2_dot_product(x,z) == std::complex<float> (0, 1));
}

TEST_CASE("evaluate quadratic form", "[f2 helper]") {
    SECTION("dimension 3") {
        std::vector<size_t> quadratic_form = {3, 6}; // x_0 x_1 + x_1 x_2

        size_t point_1 = 3; // 011
        size_t point_2 = 7; // 111

        REQUIRE(evaluate_quadratic_form(point_1, quadratic_form) == -1);
        REQUIRE(evaluate_quadratic_form(point_2, quadratic_form) == 1);
    }

    SECTION("dimension 4") {
        std::vector<size_t> quadratic_form = {3, 6, 9}; // x_0 x_1 + x_1 x_2 + x_0 x_3
        size_t point_1 = 15; // 1111
        size_t point_2 = 11; // 1011  

        REQUIRE(evaluate_quadratic_form(point_1, quadratic_form) == -1);
        REQUIRE(evaluate_quadratic_form(point_2, quadratic_form) == 1);
    }
}
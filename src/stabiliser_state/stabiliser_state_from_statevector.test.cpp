#include <catch2/catch_test_macros.hpp>

#include "stabiliser_state_from_statevector.h"
#include <vector>
#include <iostream>

using namespace fst;

std::vector<std::complex<float>> get_three_qubit_stabiliser_statevector() {
    return std::vector<std::complex<float>> {0, 0, -0.5, {0, 0.5}, 0.5, {0, 0.5}, 0, 0};
}

TEST_CASE("testing correct stabiliser states", "[statevector -> stabiliser state]") {
    SECTION("dimension 3") {
        std::vector<std::complex<float>> statevector = get_three_qubit_stabiliser_statevector();

        // for (auto elt: statevector){
        //     std::cout << elt << ", ";
        // };

        Stabiliser_From_Vector_Convertor convertor (statevector);

        REQUIRE(convertor.is_stabiliser_state);
    }
}
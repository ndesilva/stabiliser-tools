#include <catch2/catch_test_macros.hpp>

#include "stabiliser_state_from_statevector.h"
#include <vector>
#include <iostream>

using namespace fst;

std::vector<std::complex<float>> get_three_qubit_stabiliser_statevector() {
    return std::vector<std::complex<float>> {0, 0, -0.5, {0, 0.5}, 0.5, {0, 0.5}, 0, 0};
}

std::vector<std::complex<float>> get_five_qubit_stabiliser_statevector() {
    std::vector<std::complex<float>> statevector (32, 0);

    statevector[1] = 1/sqrt(8);
    statevector[7] = {0, -1/sqrt(8)};
    statevector[8] = 1/sqrt(8);
    statevector[14] = {0, -1/sqrt(8)};
    statevector[17] = -1/sqrt(8);
    statevector[23] = {0, 1/sqrt(8)};
    statevector[24] = 1/sqrt(8);
    statevector[30] = {0, -1/sqrt(8)};
    
    return statevector;
}

TEST_CASE("testing correct stabiliser states", "[statevector -> stabiliser state]") {
    SECTION("3 qubits, dimension 2") {
        std::vector<std::complex<float>> statevector = get_three_qubit_stabiliser_statevector();

        Stabiliser_From_Vector_Convertor convertor (statevector);

        REQUIRE(convertor.is_stabiliser_state);
    }

    SECTION("5 qubits, dimension 3") {
        std::vector<std::complex<float>> statevector = get_five_qubit_stabiliser_statevector();

        for (int i = 0; i<32; i++){
            std::cout << i << ": " << statevector[i] << std::endl;
        };

        Stabiliser_From_Vector_Convertor convertor (statevector);

        REQUIRE(convertor.is_stabiliser_state);
    }
}
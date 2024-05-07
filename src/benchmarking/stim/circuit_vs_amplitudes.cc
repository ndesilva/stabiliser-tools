#include "circuit_vs_amplitudes.h"

#include "vector_simulator.h"
#include "twiddle.h"

#include <limits>
#include <cassert>

using namespace stim;

inline static size_t biggest_index(const std::vector<std::complex<float>> &state_vector) {
    size_t best_index = 0;
    float best_size = std::norm(state_vector[0]);
    for (size_t k = 1; k < state_vector.size(); k++) {
        float size = std::norm(state_vector[k]);
        if (size > best_size) {
            best_size = size;
            best_index = k;
        }
    }
    return best_index;
}

inline static size_t compute_occupation(const std::vector<std::complex<float>> &state_vector) {
    size_t c = 0;
    for (const auto &v : state_vector) {
        if (v != std::complex<float>{0, 0}) {
            c++;
        }
    }
    return c;
}

bool stim::stabilizer_state_vector_to_circuit(
    const std::vector<std::complex<float>> &state_vector) {
    if (!is_power_of_2(state_vector.size())) {
        return false;
    }

    uint8_t num_qubits = floor_lg2(state_vector.size());
    double weight = 0;
    for (const auto &c : state_vector) {
        weight += std::norm(c);
    }
    if (abs(weight - 1) > 0.125) {
        return false;
    }

    VectorSimulator sim(num_qubits);
    sim.state = state_vector;

    // Move biggest amplitude to start of state vector..
    size_t pivot = biggest_index(state_vector);
    for (size_t q = 0; q < num_qubits; q++) {
        if ((pivot >> q) & 1) {
            sim.apply_X(q);
        }
    }
    sim.smooth_stabilizer_state(sim.state[0]);
    size_t occupation = compute_occupation(sim.state);
    if (!is_power_of_2(occupation)) {
        return false;
    }

    // Repeatedly cancel amplitudes
    while (occupation > 1) {
        size_t k = 1;
        for (; k < state_vector.size(); k++) {
            if (sim.state[k].real() || sim.state[k].imag()) {
                break;
            }
        }
        if (k == state_vector.size()) {
            break;
        }

        size_t base_qubit = SIZE_MAX;
        for (size_t q = 0; q < num_qubits; q++) {
            if ((k >> q) & 1) {
                if (base_qubit == SIZE_MAX) {
                    base_qubit = q;
                } else {
                    sim.apply_CX(base_qubit, q);
                }
            }
        }

		auto s = sim.state[ size_t( 1 ) << base_qubit ];
        if (s == std::complex<float>{-1, 0}) {
            sim.apply_Z(base_qubit);
        } else if (s == std::complex<float>{0, 1}) {
            sim.apply_S_DAG(base_qubit);
        } else if (s == std::complex<float>{0, -1}) {
            sim.apply_S(base_qubit);
        }
        sim.apply_H(base_qubit);

        if (!sim.smooth_stabilizer_state(sim.state[0])){
            return false;
        }
        
        if (compute_occupation(sim.state) * 2 != occupation) {
            return false;
        }
        occupation >>= 1;
    }

    return true;
}
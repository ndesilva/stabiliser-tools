// Copyright 2021 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "vector_simulator.h"

#include <cassert>
#include <iostream>

using namespace stim;

static constexpr std::complex<float> s = 0.7071067811865475244f;
static constexpr std::complex<float> i = std::complex<float>(0, 1);

VectorSimulator::VectorSimulator(size_t num_qubits) {
    state.resize(size_t{1} << num_qubits, 0.0f);
    state[0] = 1;
}

inline std::vector<std::complex<float>> mat_vec_mul(
    const std::vector<std::vector<std::complex<float>>> &matrix, const std::vector<std::complex<float>> &vec) {
    std::vector<std::complex<float>> result;
    for (size_t row = 0; row < vec.size(); row++) {
        std::complex<float> v = 0;
        for (size_t col = 0; col < vec.size(); col++) {
            v += matrix[row][col] * vec[col];
        }
        result.push_back(v);
    }
    return result;
}

void VectorSimulator::apply(
    const std::vector<std::vector<std::complex<float>>> &matrix, const std::vector<size_t> &qubits) {
    size_t n = size_t{1} << qubits.size();
    assert(matrix.size() == n);
    std::vector<size_t> masks;
    for (size_t k = 0; k < n; k++) {
        size_t m = 0;
        for (size_t q = 0; q < qubits.size(); q++) {
            if ((k >> q) & 1) {
                m |= size_t{1} << qubits[q];
            }
        }
        masks.push_back(m);
    }
    assert(masks.back() < state.size());
    for (size_t base = 0; base < state.size(); base++) {
        if (base & masks.back()) {
            continue;
        }
        std::vector<std::complex<float>> in;
        in.reserve(masks.size());
        for (auto m : masks) {
            in.push_back(state[base | m]);
        }
        auto out = mat_vec_mul(matrix, in);
        for (size_t k = 0; k < masks.size(); k++) {
            state[base | masks[k]] = out[k];
        }
    }
}

void VectorSimulator::apply_X(size_t qubit) {
    apply({{0, 1}, {1, 0}}, {qubit});
}

void VectorSimulator::apply_Z(size_t qubit) {
    apply({{1, 0}, {0, -1}}, {qubit});
}

void VectorSimulator::apply_S(size_t qubit) {
    apply({{1, 0}, {0, i}}, {qubit});
}

void VectorSimulator::apply_S_DAG(size_t qubit) {
    apply({{1, 0}, {0, -i}}, {qubit});
}

void VectorSimulator::apply_H(size_t qubit) {
    apply({{s, s}, {s, -s}}, {qubit});
}

void VectorSimulator::apply_CX(size_t qubit1, size_t qubit2) {
    apply({{1, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}, {0, 1, 0, 0}}, {qubit1, qubit2});
}

bool VectorSimulator::smooth_stabilizer_state(std::complex<float> base_value) {
    std::vector<std::complex<float>> ratio_values{
        {0, 0},
        {1, 0},
        {-1, 0},
        {0, 1},
        {0, -1},
    };
    for (size_t k = 0; k < state.size(); k++) {
        auto ratio = state[k] / base_value;
        bool solved = false;
        for (const auto &r : ratio_values) {
            if (std::norm(ratio - r) < 0.125) {
                state[k] = r;
                solved = true;
            }
        }
        if (!solved) {
            return false;
        }
    }

    return true;
}
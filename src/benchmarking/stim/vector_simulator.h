/*
 * Copyright 2021 Google LLC
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef _STIM_SIMULATORS_VECTOR_SIMULATOR_H
#define _STIM_SIMULATORS_VECTOR_SIMULATOR_H

#include <complex>
#include <iostream>
#include <random>
#include <unordered_map>
#include "gates.h"

namespace stim {

/// A state vector quantum circuit simulator.
///
/// Not intended to be particularly performant. Mostly used as a reference when testing.
struct VectorSimulator {
    std::vector<std::complex<float>> state;

    /// Creates a state vector for the given number of qubits, initialized to the zero state.
    explicit VectorSimulator(size_t num_qubits);

    /// Applies a unitary operation to the given qubits, updating the state vector.
    void apply(const std::vector<std::vector<std::complex<float>>> &matrix, const std::vector<size_t> &qubits);

    /// Helper method for applying named single qubit gates.
    void apply(GateType gate, size_t qubit);
    /// Helper method for applying named two qubit gates.
    void apply(GateType gate, size_t qubit1, size_t qubit2);

    /// Modifies the state vector to be EXACTLY entries of 0, 1, -1, i, or -i.
    ///
    /// Each entry is a ratio relative to the given base value.
    /// If any entry has a ratio not near the desired set, an exception is raised.
    void smooth_stabilizer_state(std::complex<float> base_value);
};

}  // namespace stim

#endif

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

#ifndef _STIM_GATES_GATE_DATA_H
#define _STIM_GATES_GATE_DATA_H

#include <cassert>
#include <complex>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <vector>

#include "stim/mem/fixed_cap_vector.h"

namespace stim {

template <size_t W>
struct Tableau;

template <size_t W>
struct Flow;

template <size_t W>
struct PauliString;

constexpr inline uint16_t gate_name_to_hash(const char *v, size_t n) {
    // HACK: A collision is considered to be an error.
    // Just do *anything* that makes all the defined gates have different values.

    size_t result = n;
    if (n > 0) {
        auto c_first = v[0] | 0x20;
        auto c_last = v[n - 1] | 0x20;
        result += c_first ^ (c_last << 1);
    }
    if (n > 2) {
        auto c1 = v[1] | 0x20;
        auto c2 = v[2] | 0x20;
        result ^= c1;
        result += c2 * 11;
    }
    if (n > 5) {
        auto c3 = v[3] | 0x20;
        auto c5 = v[5] | 0x20;
        result ^= c3 * 61;
        result += c5 * 77;
    }
    return result & 0x1FF;
}

constexpr inline uint16_t gate_name_to_hash(const char *c) {
    return gate_name_to_hash(c, std::char_traits<char>::length(c));
}

constexpr const size_t NUM_DEFINED_GATES = 7;

enum class GateType : uint8_t {
    NOT_A_GATE = 0,
    X,
    CX,
    Z,
    S_DAG,
    S,
    H
};

struct Gate {
    /// The gate's type, such as stim::GateType::X or stim::GateType::MRZ.
    GateType id;
    
    /// A unitary matrix describing the gate. (Size 0 if the gate is not unitary.)
    FixedCapVector<FixedCapVector<std::complex<float>, 4>, 4> unitary_data;

    std::vector<std::vector<std::complex<float>>> unitary() const;
};

struct GateDataMapHashEntry {
    GateType id;
    const char *expected_name;
    size_t expected_name_len;
};

struct GateDataMap {
   private:
    void add_gate(bool &failed, const Gate &data);
    void add_gate_alias(bool &failed, const char *alt_name, const char *canon_name);

   public:
    std::array<GateDataMapHashEntry, 512> hashed_name_to_gate_type_table;
    std::array<Gate, NUM_DEFINED_GATES> items;
    GateDataMap();

    inline const Gate &operator[](GateType g) const {
        return items[(uint64_t)g];
    }

    inline const Gate &at(const char *text, size_t text_len) const {
        auto h = gate_name_to_hash(text, text_len);
        const auto &entry = hashed_name_to_gate_type_table[h];
        if (_case_insensitive_mismatch(text, text_len, entry.expected_name, entry.expected_name_len)) {
            throw std::out_of_range("Gate not found: '" + std::string(text, text_len) + "'");
        }
        // Canonicalize.
        return (*this)[entry.id];
    }

    inline const Gate &at(const char *text) const {
        return at(text, strlen(text));
    }

    inline const Gate &at(const std::string &text) const {
        return at(text.data(), text.size());
    }

    inline bool has(const std::string &text) const {
        auto h = gate_name_to_hash(text.data(), text.size());
        const auto &entry = hashed_name_to_gate_type_table[h];
        return !_case_insensitive_mismatch(text.data(), text.size(), entry.expected_name, entry.expected_name_len);
    }
};

extern const GateDataMap GATE_DATA;

}  // namespace stim

#endif

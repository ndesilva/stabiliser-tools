#ifndef _STIM_UTIL_TOP_CIRCUIT_VS_AMPLITUDES_H
#define _STIM_UTIL_TOP_CIRCUIT_VS_AMPLITUDES_H

#include <complex>
#include <vector>

namespace stim
{

    /// Synthesizes a circuit to generate the given state vector.
    ///
    /// Args:
    ///     stabilizer_state_vector: The vector of amplitudes to produce using a circuit.
    ///     little_endian: Whether the vector is using little endian or big endian ordering.
    ///     inverted_circuit: If false, returns a circuit that sends |000...0> to the state vector.
    ///         If true, returns a circuit that sends the state vector to |000...0> instead of a cir.
    ///
    /// Returns:
    ///     A circuit that outputs the given state vector (up to global phase).
    ///
    /// Throws:
    ///     std::invalid_argument: The given state vector cannot be produced by a stabilizer circuit.
    bool stabilizer_state_vector_to_circuit(
        const std::vector<std::complex<float>> &stabilizer_state_vector, const bool assume_valid = false
        );

} // namespace stim

#endif

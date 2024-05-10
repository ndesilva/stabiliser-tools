#ifndef _FAST_STABILISER_PAULI_H
#define _FAST_STABILISER_PAULI_H

#include <complex>
#include <vector>

namespace fst
{
    struct Pauli
    {
        int number_qubits = 0;
        std::size_t x_vector = 0;
        std::size_t z_vector = 0;

        bool sign_bit = 0;
        bool imag_bit = 0;
        std::complex<float> phase = 1;

        Pauli(const int number_qubits, const std::size_t x_vector, const std::size_t z_vector, const bool sign_bit, const bool imag_bit);

        void update_phase();

        bool is_hermitian() const;
        
        bool commutes_with(const Pauli &other_pauli) const;
        bool anticommutes_with(const Pauli &other_pauli) const;

        std::vector<std::vector<std::complex<float>>> get_matrix() const;
        std::vector<std::complex<float>> multiply_vector (const std::complex<float> &vector) const;

        void multiply_by_pauli_on_right(const Pauli &other_pauli);
    };
}

#endif
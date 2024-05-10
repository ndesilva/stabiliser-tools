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

        void update_phase();

        std::vector<std::vector<std::complex<float>>> get_matrix() const;

        void multiply_by_pauli_on_right(const Pauli &other_pauli);

        std::vector<std::complex<float>> multiply_vector (const std::complex<float> &vector) const;

        bool is_hermitian() const;
        
        bool commutes_with(const Pauli &other_pauli) const;
        bool anti_commutes_with(const Pauli &other_pauli) const;
    };
}

#endif
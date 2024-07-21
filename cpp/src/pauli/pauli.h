#ifndef _FAST_STABILISER_PAULI_H
#define _FAST_STABILISER_PAULI_H

#include <complex>
#include <vector>

//TODO: make sign_bit and imag_bit bools for memory efficiency. Update f2_dot_product etc. to also return bools
namespace fst
{
    /// The class used to represent a Pauli operator.
    /// Puali is (-1)^(sign_bit) * (-i)^(imag_bit) * X^(x_vector) * Z^(z_vector)
    /// The phase of the Pauli is (-1)^(sign_bit) * (-i)^(imag_bit)
    struct Pauli
    {
        std::size_t number_qubits = 0;
        std::size_t x_vector = 0;
        std::size_t z_vector = 0;

        unsigned int sign_bit = 0;
        unsigned int imag_bit = 0;

        public:
        Pauli() = default;
        Pauli(const std::size_t number_qubits, const std::size_t x_vector, const std::size_t z_vector, const bool sign_bit, const bool imag_bit);

        /// Returns whether the pauli operator is Hermitian
        bool is_hermitian() const;

        /// Given another Paulis, used to check whether it commutes/anticommutes
        /// with this Pauli
        bool commutes_with(const Pauli &other_pauli) const;
        bool anticommutes_with(const Pauli &other_pauli) const;

        /// Given a statevector x on the same number of qubits as the Pauli P, check
        /// whether or not Px = (-1)^(eig_sign) x, i.e. whether x is an eigenstate of P
        /// with eigenvalue (-1)^(eig_sign).
        bool has_eigenstate(const std::vector<std::complex<float>> &vector, const unsigned int eign_sign) const;

        /// Returns the matrix of the Pauli (with respect to the computational basis)
        std::vector<std::vector<std::complex<float>>> get_matrix() const;

        /// Given a vector x on the same number of qubits as the Pauli P, return Px
        std::vector<std::complex<float>> multiply_vector(const std::vector<std::complex<float>> &vector) const;

        /// Given another pauli Q, multiply this Pauli on the right by Q
        /// Note, the current instance is set to the result.
        void multiply_by_pauli_on_right(const Pauli &other_pauli);

        /// Gets the current phase of the pauli: (-1)^(sign_bit) * (-i)^(imag_bit)
        std::complex<float> get_phase() const;

        bool operator==(const Pauli &other) const = default;
    };
}

#endif
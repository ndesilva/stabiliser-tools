#include "pauli.h"
#include "util/f2_helper.h"

namespace fst
{
    Pauli::Pauli(const std::size_t number_qubits, const std::size_t x_vector, const std::size_t z_vector, const bool sign_bit, const bool imag_bit)
        : number_qubits(number_qubits), x_vector(x_vector), z_vector(z_vector), sign_bit(sign_bit), imag_bit(imag_bit)
    {}

    std::complex<float> Pauli::get_phase() const
    {
        return {f_min1_pow(sign_bit) * float_not(imag_bit), -f_min1_pow(sign_bit) * static_cast<float>(imag_bit)};
    }

    bool Pauli::is_hermitian() const
    {
        return imag_bit == f2_dot_product(x_vector, z_vector);
    }

    bool Pauli::anticommutes_with(const Pauli &other_pauli) const
    {
        return f2_dot_product(x_vector, other_pauli.z_vector) ^ f2_dot_product(z_vector, other_pauli.x_vector);
    }

    bool Pauli::commutes_with(const Pauli &other_pauli) const
    {
        return !anticommutes_with(other_pauli);
    }

    std::vector<std::vector<std::complex<float>>> Pauli::get_matrix() const
    {
        const std::size_t size = integral_pow_2(number_qubits);
        std::vector<std::vector<std::complex<float>>> matrix(size, std::vector<std::complex<float>>(size, 0));

        std::complex<float> phase = get_phase(); 

        for (size_t col_index = 0; col_index < size; col_index++)
        {
            matrix[col_index ^ x_vector][col_index] =
                phase * sign_f2_dot_product(col_index, z_vector);
        }

        return matrix;
    }

    std::vector<std::complex<float>> Pauli::multiply_vector(const std::vector<std::complex<float>> &vector) const
    {
        if (integral_pow_2(number_qubits) != vector.size())
        {
            throw std::invalid_argument("Invalid vector dimension for pauli-vector multiplication");
        }
        
        const size_t size = vector.size();
        std::vector<std::complex<float>> result(size, 0);
        std::complex<float> phase = get_phase();

        for (size_t index = 0; index < size; index++)
        {
            result[index ^ x_vector] = phase * sign_f2_dot_product(index, z_vector) * vector[index];
        }

        return result;
    }

    void Pauli::multiply_by_pauli_on_right(const Pauli &other_pauli)
    {
        if (number_qubits != other_pauli.number_qubits)
        {
            throw std::invalid_argument("Paulis act on a different number of qubits");
        }
        
        std::size_t sign_bit_update = f2_dot_product(z_vector, other_pauli.x_vector) ^ (imag_bit & other_pauli.imag_bit) ^ other_pauli.sign_bit;
        imag_bit ^= other_pauli.imag_bit;
        sign_bit ^= sign_bit_update;

        x_vector ^= other_pauli.x_vector;
        z_vector ^= other_pauli.z_vector;
    }

    bool Pauli::has_eigenstate(const std::vector<std::complex<float>> &vector, const unsigned int eig_sign) const
    {
        if (integral_pow_2(number_qubits) != vector.size())
        {
            throw std::invalid_argument("Invalid vector dimension");
        }
        
        std::complex<float> phase = get_phase();

        const std::size_t size = integral_pow_2(number_qubits);
        const std::complex<float> vector_phase = f_min1_pow(eig_sign) * phase;

        for (size_t index = 0; index < size; index++)
        {
            std::complex<float> expected_phase = vector_phase * sign_f2_dot_product(index, z_vector) * vector[index];
            if (vector[index ^ x_vector] != expected_phase)
            {
                return false;
            }
        }

        return true;
    }
}
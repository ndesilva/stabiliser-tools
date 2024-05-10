#include "pauli.h"
#include "f2_helper.h"

namespace fst
{
    Pauli::Pauli(const std::size_t number_qubits, const std::size_t x_vector, const std::size_t z_vector, const bool sign_bit, const bool imag_bit)
    : number_qubits(number_qubits)
    , x_vector(x_vector)
    , z_vector(z_vector)
    , sign_bit(sign_bit)
    , imag_bit(imag_bit)
    {
        update_phase();
    }

    void Pauli::update_phase()
    {
        phase = { f_min1_pow( sign_bit ) * float_not( imag_bit ), f_min1_pow( sign_bit ) * static_cast<float>( imag_bit ) };
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
        const std::size_t size = integral_pow_2( number_qubits );
        std::vector<std::vector<std::complex<float>>> matrix (size, std::vector<std::complex<float>> (size, 0));

        for ( size_t col_index = 0; col_index < size; col_index++ )
        {
            matrix.at( col_index ^ x_vector).at(col_index) = 
                phase * sign_f2_dot_product( col_index, z_vector );
        }
        
        return matrix;
    }

    std::vector<std::complex<float>> Pauli::multiply_vector (const std::vector<std::complex<float>> &vector) const
    {
        const size_t size = vector.size();
        std::vector<std::complex<float>> result(size, 0);

        for( size_t index = 0; index < size; index++)
        {
            result[index^x_vector] = phase * sign_f2_dot_product(index, z_vector) ;
        }

        return result;
    }

    void Pauli::multiply_by_pauli_on_right(const Pauli &other_pauli)
    {
        imag_bit ^= other_pauli.imag_bit;
        int sign_bit_update = f2_dot_product(z_vector, other_pauli.x_vector) ^ (imag_bit & other_pauli.imag_bit);
        sign_bit ^= sign_bit_update;
        update_phase();

        x_vector ^= other_pauli.x_vector;
        z_vector ^= other_pauli.z_vector;
    }

    bool Pauli::has_eigenstate(const std::vector<std::complex<float>> &vector, const unsigned int eig_sign) const
    {
        const std::size_t size = integral_pow_2(number_qubits);
        const std::complex<float> vector_phase = f_min1_pow( eig_sign ) * phase;

        for ( size_t index = 0; index < size; index++)
        {
            std::complex<float> expected_phase = vector_phase * sign_f2_dot_product(index, z_vector) * vector.at(index);
            if (vector.at(index ^ x_vector) != expected_phase) {
                return false;
            }
        }

        return true;
    }

}
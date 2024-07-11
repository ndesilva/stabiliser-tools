#include "clifford.h"
#include "util/f2_helper.h"
#include "stabiliser_state/check_matrix.h"
#include "stabiliser_state/stabiliser_state.h"

namespace fst
{
    Clifford::Clifford(const std::vector<Pauli> z_conjugates, const std::vector<Pauli> x_conjugates, const std::complex<float> global_phase )
        : z_conjugates(z_conjugates), x_conjugates(x_conjugates), global_phase(global_phase)
        {
            number_qubits = z_conjugates.size();
        }

    std::vector<std::vector<std::complex<float>>> Clifford::get_matrix() const
    {
        const std::size_t size = integral_pow_2(number_qubits);
        std::vector<std::vector<std::complex<float>>> transposed_matrix(size, std::vector<std::complex<float>> (size, 0) );

        Check_Matrix first_col_check_matrix(z_conjugates);
        Stabiliser_State state (first_col_check_matrix);
        state.global_phase = global_phase;

        transposed_matrix[0] = state.get_state_vector();

        std::size_t old_col_index = 0;

        for(std::size_t i = 1; i < size; i++)
        {
            // Iterate through the Gray code
            std::size_t new_col_index = i ^ (i >> 1);
            std::size_t bit_flipped = integral_log_2( old_col_index ^ new_col_index );

            transposed_matrix[new_col_index] = x_conjugates.at(bit_flipped).multiply_vector(transposed_matrix[old_col_index]);

            old_col_index = new_col_index;
        }

        std::vector<std::vector<std::complex<float>>> matrix(size, std::vector<std::complex<float>>(size, 0));
        
        for(std::size_t i = 0; i < size; i++)
        {
            for(std::size_t j = 0; j < size; j++)
            {
                matrix[i][j] = transposed_matrix[j][i];
            }
        }

        return matrix;
    }
}
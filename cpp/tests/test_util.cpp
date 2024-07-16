#include "test_util.h"
#include "util/f2_helper.h"

using namespace fst;

namespace test {

    std::vector<std::complex<float>> matrix_vector_mult(const std::vector<std::vector<std::complex<float>>> &matrix, const std::vector<std::complex<float>> &vector)
    {
        std::size_t size = vector.size();

        std::vector<std::complex<float>> result (size, 0);

        for(std::size_t i = 0; i < size; i ++)
        {
            for (std::size_t j = 0; j < size; j++) 
            {
                result[i] += matrix[i][j] * vector[j];
            }
        }

        return result;
    }

    std::unordered_map<std::size_t, bool> get_quadratic_from_from_vector(const std::size_t dimension, const std::vector<std::size_t> &non_zero_coeffs)
	{
		std::unordered_map<std::size_t, bool> quadratic_form;
		quadratic_form.reserve(dimension * (dimension + 1)/2);
		quadratic_form[0] = 0;

		for (std::size_t j = 0; j < dimension; j++)
		{
			for (std::size_t k = j + 1; k < dimension; k++)
			{
				quadratic_form[integral_pow_2(j) ^ integral_pow_2(k)] = 0;
			}
		}

		for (const auto coeff : non_zero_coeffs)
		{
			quadratic_form[coeff] = 1;
		}

		return quadratic_form;
	}

    bool vectors_row_reduced(std::vector<std::size_t> &vectors)
    {
        std::size_t length = vectors.size();
        
        for (std::size_t i = 0; i < length; i++)
        {
            std::size_t vector = vectors[i];
            
            if (vector == 0)
            {
                return false;
            }
            
            std::size_t pivot_index = integral_log_2(vector);

            for (std::size_t j = 0; j < length; j ++)
            {
                if (i != j && bit_set_at(vectors[j], pivot_index))
                {
                    return false;
                }
            }
        }

        return true;
    }

    std::size_t Pauli_Hasher::operator() (const Pauli & pauli) const
        {
            return (
                (std::hash<std::size_t>()(pauli.x_vector)
                ^ std::hash<std::size_t>()(pauli.z_vector << 1) << 2)
                ^ ((std::size_t) pauli.sign_bit << 1)
                ^ ((std::size_t) pauli.imag_bit)
            );
        }
}


#ifndef _FAST_STABILISER_STABILISER_STATE_FROM_VECTOR_H
#define _FAST_STABILISER_STABILISER_STATE_FROM_VECTOR_H

#include <vector>
#include <complex>
#include "stabiliser_state.h"

namespace fst{
    struct Stabiliser_From_Vector_Convertor {
        bool is_stabiliser_state = false;

        int number_qubits;
        int support_size;
        int dimension;
        std::vector<int> basis_vectors;
        int shift;

        int real_linear_part;
        int imaginary_part;
        std::vector<int> quadratic_form;
        std::complex<float> gloabal_factor;

        explicit Stabiliser_From_Vector_Convertor(std::vector<std::complex<float>> &state_vector, bool assume_stabiliser_state = 0);

        Stabiliser_State get_stabiliser_state() const;
    };
}

#endif
#include "stabiliser_state_from_vector.h"

#include <algorithm>
#include "f2_helper/f2_helper.h"

using namespace fst;

// TODO fix behaviour on zero vector
Stabiliser_From_Vector_Convertor::Stabiliser_From_Vector_Convertor(std::vector<std::complex<float>> &state_vector, bool assume_stabiliser_state = false) {
    int state_vector_size = state_vector.size();

    number_qubits = ceil(log2(state_vector_size));

    if (1 << number_qubits != state_vector_size){
        return;
    };

    int shift = 0;

    while (state_vector[shift] == float(0)) {
        shift++;
    }

    std::vector<int> vector_space_indices{0};

    for (int index = 0; index < state_vector_size; index++) {
        if (state_vector[index] != float(0)) {
            vector_space_indices.push_back(shift ^ index);
        }
    };

    support_size = vector_space_indices.size();
    dimension = ceil(log2(state_vector_size));

    if (1 << dimension != support_size){
        return;
    }

    float normalisation_factor = sqrt(support_size);
    std::complex<float> first_entry = state_vector[shift];
    gloabal_factor = normalisation_factor * first_entry;

    std::sort(vector_space_indices.begin(), vector_space_indices.end());

    for (int j = 0; j < dimension; j++) {
        int weight_one_string = 1 << j;

        int basis_vector = vector_space_indices[weight_one_string];
        basis_vectors.push_back(basis_vector);

        std::complex<float> phase = state_vector[basis_vector ^ shift]/first_entry;

        if (std::norm(phase - float(-1)) < 0.125) {
            real_linear_part ^= weight_one_string;
        }
        else if (std::norm(phase - std::complex<float>{0,1}) < 0.125) {
            imaginary_part ^= weight_one_string;
        }
        else if (std::norm(phase - std::complex<float>{0,-1}) < 0.125) {
            real_linear_part ^= weight_one_string;
            imaginary_part ^= weight_one_string;
        }
        else if (std::norm(phase - float(1)) >= 0.125){
            return;
        }

    for (int j=0; j < dimension; j++) {
        for (int i = j + 1; i < dimension; i++) {
            int vector_index = (1<<i) | (1<<j);

            float real_linear_eval = sign_f2_dot_product(vector_index, real_linear_part);
            std::complex<float> imag_linear_eval = imag_f2_dot_product(vector_index, imaginary_part);
            std::complex<float> linear_eval = real_linear_eval * imag_linear_eval;
            
            int total_index = vector_space_indices[vector_index]; ^ shift;

            std::complex<float> quadratic_form_eval = state_vector[total_index]/(first_entry * linear_eval);

            if (std::norm(quadratic_form_eval - float(-1)) < 0.125) {
                quadratic_form.push_back(vector_index);
            }
            else if (std::norm(quadratic_form_eval - float(-1)) >= 0.125) {
                return;
            }
        }
    }

    if (assume_stabiliser_state){
        is_stabiliser_state = true;
        return;
    }

    is_stabiliser_state = check_remaining_entries(state_vector);
}

bool Stabiliser_From_Vector_Convertor::check_remaining_entries(std::vector<std::complex<float>> &state_vector) const {
    for(int vector_index = 1; vector_index < support_size; vector_index++) {
        int total_index = vector_space_indices[vector_index] ^ shift;
        complex<float> actual_phase = state_vector[total_index];

        float real_linear_eval = sign_f2_dot_product(vector_index, real_linear_part);
        std::complex<float> imag_linear_eval = imag_f2_dot_product(vector_index, imaginary_part);
        
        // TODO change this
        std::complex<float> quadratic_eval = 1;
        std::complex<float> phase_eval = real_linear_eval * imag_linear_eval * quadratic_eval;

        if (std::norm(phase_eval * first_entry - actual_phase) >= 0.125) {
            return false;
        }
    }

    return true;
}

#include "stabiliser_state_from_statevector.h"

#include <algorithm>
#include "f2_helper/f2_helper.h"
#include <iostream>

using namespace fst;

// TODO fix behaviour on zero vector
Stabiliser_From_Vector_Convertor::Stabiliser_From_Vector_Convertor(std::vector<std::complex<float>> &statevector, bool assume_stabiliser_state) {
    std::cout << "\n calling converter\n";
    
    int state_vector_size = statevector.size();
    number_qubits = integral_log_2(state_vector_size);

    if (1 << number_qubits != state_vector_size){
        return;
    };

    std::cout << "passed power of 2 size check\n";

    shift = 0;

    while (statevector[shift] == float(0)) {
        shift++;
    }

    std::vector<int> vector_space_indices{0};

    for (int index = shift + 1; index < state_vector_size; index++) {
        if (statevector[index] != float(0)) {
            vector_space_indices.push_back(shift ^ index);
        }
    };

    support_size = vector_space_indices.size();
    dimension = integral_log_2(support_size);

    if (1 << dimension != support_size){
        return;
    }

    std::cout << "passed power of 2 support size check\n";

    float normalisation_factor = sqrt(support_size);
    first_entry = statevector[shift];
    gloabal_factor = normalisation_factor * first_entry;

    std::sort(vector_space_indices.begin(), vector_space_indices.end());

    for (int j = 0; j < dimension; j++) {
        int weight_one_string = 1 << j;

        int basis_vector = vector_space_indices[weight_one_string];
        basis_vectors.push_back(basis_vector);

        std::complex<float> phase = statevector[basis_vector ^ shift]/first_entry;

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
    }

    std::cout << "passed linear extraction\n";

    for (int j=0; j < dimension; j++) {
        for (int i = j + 1; i < dimension; i++) {
            int vector_index = (1<<i) | (1<<j);

            float real_linear_eval = sign_f2_dot_product(vector_index, real_linear_part);
            std::complex<float> imag_linear_eval = imag_f2_dot_product(vector_index, imaginary_part);
            std::complex<float> linear_eval = real_linear_eval * imag_linear_eval;
            
            int total_index = vector_space_indices[vector_index] ^ shift;

            std::complex<float> quadratic_form_eval = statevector[total_index]/(first_entry * linear_eval);

            if (std::norm(quadratic_form_eval - float(-1)) < 0.125) {
                quadratic_form.push_back(vector_index);
            }
            else if (std::norm(quadratic_form_eval - float(1)) >= 0.125) {
                return;
            }
        }
    }

    std::cout << "passed quadratic extraction\n";

    if (assume_stabiliser_state){
        is_stabiliser_state = true;
        return;
    }

    is_stabiliser_state = check_remaining_entries(statevector);

    std::cout << "remaining entries check gave " << is_stabiliser_state << std::endl;
};

bool Stabiliser_From_Vector_Convertor::check_remaining_entries(std::vector<std::complex<float>> &statevector) const {
    // std::cout << "checking remaining entires\n";
    
    int old_vector_index = 0;
    int total_index = shift;
    
    for(int i = 1; i < support_size; i++) {
        // iterate through the gray code
        int new_vector_index = i ^ (i >> 1);

        // std::cout << "checking vector with index ";
        // std::cout << new_vector_index << std::endl;

        int flipped_bit = integral_log_2(new_vector_index ^ old_vector_index);
        
        // std::cout << "shift is now ";
        // std::cout << shift << std::endl;

        total_index ^= basis_vectors[flipped_bit];

        // std::cout << "expanded vector is ";
        // std::cout << total_index << std::endl;

        std::complex<float> actual_phase = statevector[total_index];

        // std::cout << "done actual phase\n";
        float real_linear_eval = sign_f2_dot_product(new_vector_index, real_linear_part);
        
        // std::cout << "done real eval\n";

        std::complex<float> imag_linear_eval = imag_f2_dot_product(new_vector_index, imaginary_part);
        // std::cout << "done imag eval\n";
        
        std::complex<float> quadratic_eval = evaluate_quadratic_form(new_vector_index, quadratic_form);
        // std::cout << "done quadratic eval\n";
        
        std::complex<float> phase_eval = real_linear_eval * imag_linear_eval * quadratic_eval;

        // std::cout << "found phase\n";

        if (std::norm(phase_eval * first_entry - actual_phase) >= 0.125) {
            return false;
        }

        old_vector_index = new_vector_index;
    }

    return true;
}

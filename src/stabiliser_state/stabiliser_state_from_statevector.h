#ifndef _FAST_STABILISER_STABILISER_STATE_FROM_VECTOR_H
#define _FAST_STABILISER_STABILISER_STATE_FROM_VECTOR_H

#include <vector>
#include <complex>
#include "stabiliser_state.h"

namespace fst{
    /// Class used in constructing a stabiliser_state object from an explicit
    /// state vector of complex amplitudes and/or testing if a state vector 
    /// correspdonds to a stabiliser state. To convert/test, call the object constructor.
    struct Stabiliser_From_Vector_Convertor {
        bool is_stabiliser_state = false;

    private:
        int number_qubits;
        int support_size;
        int dimension;
        std::vector<int> basis_vectors;
        int shift;

        int real_linear_part = 0;
        int imaginary_part = 0;
        std::vector<int> quadratic_form;
        
        std::complex<float> global_phase;
        std::complex<float> first_entry;

    public:
        /// Convert a state vector of complex amplitudes into a stabiliser state object.
        /// if assume_stabiliser_state is set to true, the function runs much faster, but may
        /// have unexpected behaviour if the state vector is not in fact a stabiliser state.
        /// once the constructor has run, the flag is_stabiliser_state is set to true iff. the input
        /// vector was a stabiliser state. To retreive the stabiliser state object, call the get_stabiliser_state()
        /// method
        explicit Stabiliser_From_Vector_Convertor(std::vector<std::complex<float>> &statevector, bool assume_stabiliser_state = false);

        /// Extract the stabiliser_state object. Throws an error if is_stabiliser state is false  
        Stabiliser_State get_stabiliser_state() const;

    private:
        bool check_remaining_entries(std::vector<std::complex<float>> &statevector) const;
    };
}

#endif
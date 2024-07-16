#ifndef _STIM_CIRCUIT_VS_AMPLITUDES_PYBIND_H
#define _STIM_CIRCUIT_VS_AMPLITUDES_PYBIND_H

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>

#include "circuit_vs_amplitudes.h"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace stim;

// TODO: Try and export the operator == 
namespace stim_pybind {
    void init_circuit_vs_amplitudes(py::module_ &m)
    {
        m.def("stabilizer_state_vector_to_circuit", &stim::stabilizer_state_vector_to_circuit, "stabilizer_state_vector"_a);
    }
}

#endif
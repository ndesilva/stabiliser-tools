#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>

#include "circuit_vs_amplitudes.h"

namespace py = pybind11;
using namespace stim;

namespace stim_pybind {
    PYBIND11_MODULE(stim_mock, m)
    {
        m.def("circuit_from_statevector", &stabilizer_state_vector_to_circuit, py::arg("statevector"), py::arg("assume_valid") = false);
    }
}
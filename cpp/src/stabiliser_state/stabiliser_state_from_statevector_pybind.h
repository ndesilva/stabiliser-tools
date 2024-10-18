#ifndef _FAST_STABILISER_STABILISER_STATE_FROM_VECTOR_PYBIND_H
#define _FAST_STABILISER_STABILISER_STATE_FROM_VECTOR_PYBIND_H

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>

#include "stabiliser_state_from_statevector.h"

namespace py = pybind11;
using namespace fst;

// Eliminate copying by making types opaque, accept np arrays?
namespace fst_pybind
{
    void init_stabiliser_state_from_statevector(py::module_ &m)
    {
        m.def("stabiliser_state_from_statevector", &stabiliser_from_statevector, py::arg("statevector"), py::arg("assume_valid") = false, "Converts a state vector of complex amplitudes into a stabiliser state object. Assuming valid is faster, but will result in undefined behaviour if the state vector is not in fact a valid stabiliser state");
        m.def("is_stabiliser_state", &is_stabiliser_state, py::arg("statevector"), "Tests whether a state vector of complex amplitudes corresponds to a stabiliser state");
        m.def("stab_in_the_dark", &stab_in_the_dark, py::arg("statevector"), ";)");
    }
}

#endif
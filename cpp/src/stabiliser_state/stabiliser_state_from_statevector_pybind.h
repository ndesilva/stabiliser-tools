#ifndef _FAST_STABILISER_STABILISER_STATE_FROM_VECTOR_PYBIND_H
#define _FAST_STABILISER_STABILISER_STATE_FROM_VECTOR_PYBIND_H

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>

#include "stabiliser_state_from_statevector.h"

namespace py = pybind11;
using namespace fst;

// Eliminate copying by making types opaque?
namespace fst_pybind {
    void init_stabiliser_state_from_statevector(py::module_ &m)
    {
        m.def("stabiliser_state_from_statevector", &stabiliser_from_statevector, py::arg("statevector"), py::arg("assume_valid") = false);
        m.def("is_stabiliser_state", &is_stabiliser_state);
    }
}

#endif
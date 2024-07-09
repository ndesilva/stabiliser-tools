#ifndef _FAST_STABILISER_CHECK_MATRIX_PYBIND_H
#define _FAST_STABILISER_CHECK_MATRIX_PYBIND_H

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>

#include "check_matrix.h"

namespace py = pybind11;
using namespace fst;

// TODO: Return numpy arrays rather than arrays. Eliminate copying by making types opaque?
namespace fst_pybind {
    void init_check_matrix(py::module_ &m)
    {
        py::class_<Check_Matrix>(m, "Check_Matrix")
            .def_readwrite("number_qubits", &Check_Matrix::number_qubits);
    }
}

#endif
#ifndef _FAST_STABILISER_CLIFFORD_FROM_MATRIX_PYBIND_H
#define _FAST_STABILISER_CLIFFORD_FROM_MATRIX_PYBIND_H

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>

#include "clifford_from_matrix.h"

namespace py = pybind11;
using namespace fst;

// Eliminate copying by making types opaque, accept np arrays?

namespace fst_pybind
{
    void init_clifford_from_matrix(py::module_ &m)
    {
        m.def("clifford_from_matrix", &clifford_from_matrix, py::arg("matrix"), py::arg("assume_valid") = false, "Converts a 2^n by 2^n matrix with complex entries into a Clifford object. Assuming valid is faster, but will result in undefined behaviour if the matrix is not in fact a valid Clifford operator");
        m.def("is_clifford_matrix", &is_clifford_matrix, py::arg("matrix"), "Tests whether a matrix with complex entries corresponds to a Clifford");
    }
}

#endif
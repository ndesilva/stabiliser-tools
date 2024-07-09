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
            .def_readwrite("number_qubits", &Check_Matrix::number_qubits)
            .def_readwrite("paulis", &Check_Matrix::paulis)
            .def_readwrite("row_reduced", &Check_Matrix::row_reduced)
            .def(py::init<const std::vector<Pauli>, const bool>(), py::arg("paulis"), py::arg("row_reduced") = false)
            .def(py::init<Stabiliser_State&>())
            .def("get_state_vector", &Check_Matrix::get_state_vector)
            .def("row_reduce", &Check_Matrix::row_reduce);
    }
}

#endif
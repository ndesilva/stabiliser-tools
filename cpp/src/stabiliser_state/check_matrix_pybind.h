#ifndef _FAST_STABILISER_CHECK_MATRIX_PYBIND_H
#define _FAST_STABILISER_CHECK_MATRIX_PYBIND_H

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>

#include "check_matrix.h"

namespace py = pybind11;
using namespace fst;

// TODO: Return numpy arrays rather than arrays. Eliminate copying by making types opaque?
namespace fst_pybind
{
    void init_check_matrix(py::module_ &m)
    {
        py::class_<Check_Matrix>(m, "Check_Matrix")
            .def_readwrite("number_qubits", &Check_Matrix::number_qubits, "int\t\tThe number of qubits")
            .def_readwrite("row_reduced", &Check_Matrix::row_reduced, "bool")
            .def("set_paulis", &Check_Matrix::set_paulis, py::arg("paulis"), "Sets the list of stabilisers for the stabiliser state")
            .def("get_paulis", &Check_Matrix::get_paulis, "Gets the list[Pauli] of stabilisers for the stabiliser state")
            .def(py::init<const std::vector<Pauli>, const bool>(), py::arg("paulis"), py::arg("row_reduced") = false)
            .def(py::init<Stabiliser_State &>(), py::arg("stabiliser_state"))
            .def("get_state_vector", &Check_Matrix::get_state_vector, "Returns the state vector of length 2^n stabilised by each of the Paulis in the check matrix")
            .def("row_reduce", &Check_Matrix::row_reduce, "Row reduces the check matrix, giving a new set of Paulis that generates the same stabiliser group.\n\nPaulis are sorted into 2 types: \"z_only\", which have no X component, and \"x_stabilisers\", which may have both an x and z component. After performing this function, the x_vectors of the new \"x_stabiliser\" Paulis and the z_vectors of the new \"z_only\" stabilisers are in reduced row echelon form. Note that the collection of all the Paulis' z_vectors may NOT be in reduced row echelon form")
            .doc() = "The class used to represent a list of n commuting Paulis, an alternative representation of a stabiliser state";
    }
}

#endif
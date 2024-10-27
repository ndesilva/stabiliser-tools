#ifndef _FAST_STABILISER_STABILISER_STATE_PYBIND_H
#define _FAST_STABILISER_STABILISER_STATE_PYBIND_H

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>

#include "stabiliser_state.h"

namespace py = pybind11;
using namespace fst;
using namespace pybind11::literals;

// TODO: Try and export the operator ==
namespace fst_pybind
{
    void init_stabiliser_state(py::module_ &m)
    {
        py::class_<Stabiliser_State>(m, "Stabiliser_State")
            .def_readwrite("number_qubits", &Stabiliser_State::number_qubits, "int\t\tThe number of qubits")
            .def_readwrite("basis_vectors", &Stabiliser_State::basis_vectors, "list[int]\tBasis vectors for the vector space")
            .def_readwrite("dim", &Stabiliser_State::dim, "int\t\tThe dimension of the vector space")
            .def_readwrite("shift", &Stabiliser_State::shift, "int\t\tA constant vector that shifts the vector space to the affine space")
            .def_readwrite("real_linear_part", &Stabiliser_State::real_linear_part, "int\t\tThe diagonal part of the quadratic form")
            .def_readwrite("imaginary_part", &Stabiliser_State::imaginary_part, "int\t\tThe linear form")
            .def_readwrite("quadratic_form", &Stabiliser_State::quadratic_form, "dict[int]\tThe (rest of the) quadratic form. The quadratic form is stored as a map. It should always have quadratic_form[0] = 0. Q(e_i, e_j) is stored as quadratic_form[2^i ^ 2^j]")
            .def_readwrite("global_phase", &Stabiliser_State::global_phase, "complex\t\tThe global phase")
            .def_readwrite("row_reduced", &Stabiliser_State::row_reduced, "bool\t\tWhether the matrix of basis vectors is row reduced")
            .def(py::init<const std::size_t>(), "number_qubits"_a) // TODO: Do we want this?
            .def(py::init<Check_Matrix &>(), "check_matrix"_a)
            .def("get_state_vector", &Stabiliser_State::get_state_vector, "Returns the state vector of length 2^n of the stabiliser state (with respect to the computational basis), as type list[complex]")
            .def("row_reduce_basis", &Stabiliser_State::row_reduce_basis, "Row reduces the basis to reduced row-echelon form. Note that the quadratic form and the real and imaginary linear parts are also updated, so the instance represents the same stabiliser state")
            .doc() = "The class used to represent a stabiliser state. The state is stored using the ideas of Dehaene & De Moore, as an affine space, and a quadratic and linear form over that space. More precisely, it is stored as a list of basis vectors for a vector space, a constant vector that is added to every element of the vector space to reach, the affine space, and a quadratic and linear form defined on the vector space";
    }
}

#endif
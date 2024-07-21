#ifndef _FAST_STABILISER_STABILISER_STATE_PYBIND_H
#define _FAST_STABILISER_STABILISER_STATE_PYBIND_H

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>

#include "stabiliser_state.h"

namespace py = pybind11;
using namespace fst;


// TODO: Try and export the operator == 
namespace fst_pybind
{
    void init_stabiliser_state(py::module_ &m)
    {
        py::class_<Stabiliser_State>(m, "Stabiliser_State")
            .def_readwrite("number_qubits", &Stabiliser_State::number_qubits)
            .def_readwrite("basis_vectors", &Stabiliser_State::basis_vectors)
            .def_readwrite("dim", &Stabiliser_State::dim)
            .def_readwrite("shift", &Stabiliser_State::shift)
            .def_readwrite("real_linear_part", &Stabiliser_State::real_linear_part)
            .def_readwrite("imaginary_part", &Stabiliser_State::imaginary_part)
            .def_readwrite("quadratic_form", &Stabiliser_State::quadratic_form)
            .def_readwrite("global_phase", &Stabiliser_State::global_phase)
            .def_readwrite("row_reduced", &Stabiliser_State::row_reduced)
            .def(py::init<const std::size_t, const std::size_t>())
            .def(py::init<const std::size_t>())
            .def(py::init<Check_Matrix&>())
            .def("get_state_vector", &Stabiliser_State::get_state_vector)
            .def("row_reduce_basis", &Stabiliser_State::row_reduce_basis);
    }
}

#endif
#ifndef _FAST_STABILISER_PAULI_PYBIND_H
#define _FAST_STABILISER_PAULI_PYBIND_H

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>

#include "pauli.h"

namespace py = pybind11;
using namespace fst;

// TODO: Export the operator ==
// TODO: Return numpy arrays rather than arrays. Eliminate copying by making types opaque?
namespace fst_pybind {
    void init_pauli(py::module_ &m)
    {
        py::class_<Pauli>(m, "Pauli")
            .def_readwrite("number_qubits", &Pauli::number_qubits)
            .def_readwrite("x_vector", &Pauli::x_vector)
            .def_readwrite("z_vector", &Pauli::z_vector)
            .def_readwrite("sign_bit", &Pauli::sign_bit)
            .def_readwrite("imag_bit", &Pauli::imag_bit)
            .def(py::init<const std::size_t, const std::size_t, const std::size_t, const bool, const bool>())
            .def("is_hermitian", &Pauli::is_hermitian)
            .def("commutes_with", &Pauli::commutes_with)
            .def("anticommutes_with", &Pauli::anticommutes_with)
            .def("has_eigenstate", &Pauli::has_eigenstate)
            .def("get_matrix", &Pauli::get_matrix)
            .def("multiply_vector", &Pauli::multiply_vector)
            .def("multiply_by_pauli_on_right", &Pauli::multiply_by_pauli_on_right)
            .def("get_phase", &Pauli::get_phase);
    }
}

#endif
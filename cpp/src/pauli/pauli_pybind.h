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
namespace fst_pybind
{
    void init_pauli(py::module_ &m)
    {
        py::class_<Pauli>(m, "Pauli")
            .def_readwrite("number_qubits", &Pauli::number_qubits, "int\t\tThe number of qubits")
            .def_readwrite("x_vector", &Pauli::x_vector, "int")
            .def_readwrite("z_vector", &Pauli::z_vector, "int")
            .def_readwrite("sign_bit", &Pauli::sign_bit, "int")
            .def_readwrite("imag_bit", &Pauli::imag_bit, "int")
            .def(py::init<const std::size_t, const std::size_t, const std::size_t, const bool, const bool>(), py::arg("number_qubits"), py::arg("x_vector"), py::arg("z_vector"), py::arg("sign_bit"), py::arg("imag_bit"))
            .def("is_hermitian", &Pauli::is_hermitian, "Returns whether the pauli operator is Hermitian")
            .def("commutes_with", &Pauli::commutes_with, py::arg("other_pauli"), "Given another Pauli, used to check whether it commutes with this Pauli")
            .def("anticommutes_with", &Pauli::anticommutes_with, py::arg("other_pauli"), "Given another Pauli, used to check whether it anticommutes with this Pauli")
            .def("has_eigenstate", &Pauli::has_eigenstate, py::arg("vector"), py::arg("eig_sign"), "Given a statevector x on the same number of qubits as the Pauli P, checks whether or not Px = (-1)^(eig_sign) x, i.e. whether x is an eigenstate of P with eigenvalue (-1)^(eig_sign)")
            .def("get_matrix", &Pauli::get_matrix, "Returns the matrix of the Pauli (with respect to the computational basis)")
            .def("multiply_vector", &Pauli::multiply_vector, py::arg("vector"), "Given a vector x on the same number of qubits as the Pauli P, returns Px")
            .def("multiply_by_pauli_on_right", &Pauli::multiply_by_pauli_on_right, py::arg("other_pauli"), "Given another pauli Q, multiplies this Pauli on the right by Q. Note, the current instance is set to the result")
            .def("get_phase", &Pauli::get_phase, "Gets the current phase of the pauli: (-1)^(sign_bit) * (-i)^(imag_bit)")
            .doc() = "The class used to represent a Pauli operator. A Pauli is (-1)^(sign_bit) * (-i)^(imag_bit) * X^(x_vector) * Z^(z_vector). The phase of the Pauli is (-1)^(sign_bit) * (-i)^(imag_bit)";
    }
}

#endif
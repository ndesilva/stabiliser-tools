#ifndef _FAST_STABILISER_CLIFFORD_PYBIND_H
#define _FAST_STABILISER_CLIFFORD_PYBIND_H

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>

#include "clifford.h"

namespace py = pybind11;
using namespace fst;

namespace fst_pybind
{
    void init_clifford(py::module_ &m)
    {
        py::class_<Clifford>(m, "Clifford")
            .def_readwrite("number_qubits", &Clifford::number_qubits, "int\t\tThe number of qubits")
            .def_readwrite("z_conjugates", &Clifford::z_conjugates, "list[Pauli]")
            .def_readwrite("x_conjugates", &Clifford::x_conjugates, "list[Pauli]")
            .def_readwrite("global_phase", &Clifford::global_phase, "complex")
            .def(py::init<const std::vector<Pauli>, const std::vector<Pauli>, const std::complex<float>>(), py::arg("z_conjugates"), py::arg("x_conjugates"), py::arg("global_phase") = 1.0f)
            .def("get_matrix", &Clifford::get_matrix, "Returns the matrix of the Clifford (with respect to the computational basis)")
            .doc() = "The class used to represent a Clifford operator U. Represented by its action on the Pauli basis: z_conjugates[i] = UZ_iU*, x_conjugates[i] = UX_iU*";
    }
}

#endif
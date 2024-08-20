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
            .def_readwrite("number_qubits", &Clifford::number_qubits)
            .def_readwrite("z_conjugates", &Clifford::z_conjugates)
            .def_readwrite("x_conjugates", &Clifford::x_conjugates)
            .def_readwrite("global_phase", &Clifford::global_phase)
            .def(py::init<const std::vector<Pauli>, const std::vector<Pauli>, const std::complex<float>>(), py::arg("z_conjugates"), py::arg("x_conjugates"), py::arg("global_phase") = 1.0f )
            .def("get_matrix", &Clifford::get_matrix);
    }
}

#endif
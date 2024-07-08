#include <pybind11/pybind11.h>

#include "pauli.h"

namespace py = pybind11;
using namespace fst;

PYBIND11_MODULE(fast, m)
{
    py::class_<Pauli>(m, "Pauli")
        .def(py::init<const std::size_t, const std::size_t, const std::size_t, const bool, const bool>())
        .def("is_hermitian", &Pauli::is_hermitian);
}
#include <pybind11/pybind11.h>

#include "pauli/pauli_pybind.h"
#include "stabiliser_state/check_matrix_pybind.h"
#include "stabiliser_state/stabiliser_state_pybind.h"
#include "stabiliser_state/stabiliser_state_from_statevector_pybind.h"
#include "clifford/clifford_pybind.h"
#include "clifford/clifford_from_matrix_pybind.h"

namespace py = pybind11;
using namespace fst;

namespace fst_pybind {

    void init_pauli(py::module_ &);
    void init_check_matrix(py::module_ &);
    void init_stabiliser_state(py::module_ &);
    void init_stabiliser_state_from_statevector(py::module_ &);
    void init_clifford(py::module_ &);
    void init_clifford_from_matrix(py::module_ &);
    
    PYBIND11_MODULE(stab_tools, m)
    {
        init_pauli(m);
        init_check_matrix(m);
        init_stabiliser_state(m);
        init_stabiliser_state_from_statevector(m);
        init_clifford(m);
        init_clifford_from_matrix(m);
    }
}
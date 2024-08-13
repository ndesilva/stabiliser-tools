import generators as gs
import qiskit.quantum_info as qi
import stim

import sys
PATH_TO_LIBRARY = './build/ninja-multi-vcpkg/cpp/src/Release'
PATH_TO_STIM_MOCK = './build/ninja-multi-vcpkg/cpp/benchmarking/stim/Release'
sys.path.extend([PATH_TO_LIBRARY, PATH_TO_STIM_MOCK])
import fast as fst
import stim_mock as sm

def qiskit_C1_converter(matrix, assume_valid = True):
    return qi.Clifford.from_matrix(matrix)

def stim_C1_convertor(matrix, assume_valid = True):
    return stim.Tableau.from_unitary_matrix(matrix, endian = "little")

def qiskit_C1_test(matrix):
    try:
        qi.Clifford.from_matrix(matrix)
    except:
        pass

def stim_C1_test(matrix):
    try:
        stim.Tableau.from_unitary_matrix(matrix, endian = "little")
    except:
        pass


configs = [
    {
        "pre_string": "converting S1 to efficient rep",
        "functions_to_time": [
            sm.circuit_from_statevector,
            fst.stabiliser_state_from_statevector
        ],
        "function_strings": [
            "stim",
            "our method"
        ],
        "generation_types": [
            gs.random_stab_state_with_assump,
            gs.random_stab_state,
            gs.computational_zero,
            gs.random_full_support_stab_state,
        ],
        "generation_strings": [
            "random stab state with assump",
            "random stab state without assump",
            "computational zero",
            "random full support stabiliser state"
        ],
        "min_qubit_number" : 1,
        "max_qubit_number" : 11,
        "reps" : int(1e3)
    },

    {
        "pre_string": "testing S1",
        "functions_to_time": [
            sm.circuit_from_statevector,
            fst.is_stabiliser_state
        ],
        "function_strings": [
            "stim",
            "our method"
        ],
        "generation_types": [
            gs.random_stab_state,
            gs.random_full_support_almost_stab_state
        ],
        "generation_strings": [
            "random stabiliser state",
            "almost stab state"
        ],
        "min_qubit_number" : 1,
        "max_qubit_number" : 11,
        "reps" : int(1e3)
    },

    {
        "pre_string": "converting C1 to efficient rep",
        "functions_to_time": [
            fst.clifford_from_matrix,
            qiskit_C1_converter,
            stim_C1_convertor
        ],
        "function_strings": [
            "our method",
            "Qiskit",
            "stim"
        ],
        "generation_types": [
            gs.random_clifford,
            gs.random_clifford_with_assumption,
            gs.get_identity_matrix,
            gs.get_Hadamard_matrix,
            gs.get_anti_identiy_matrix,
        ],
        "generation_strings": [
            "random clifford without assump",
            "random clifford with assump",
            "identity matrix",
            "Hadamard matrix",
            "anti-identity matrix"
        ],
        "min_qubit_number" : 1,
        "max_qubit_number" : 9,
        "reps" : int(1e3) 
    },

    {
        "pre_string": "testing C1",
        "functions_to_time":[
            fst.is_clifford_matrix,
            qiskit_C1_test,
            stim_C1_test
        ],
        "function_strings": [
            "our method",
            "Qiskit",
            "stim"
        ],
        "generation_types": [
            gs.random_clifford,
            gs.random_almost_clifford,
        ],
        "generation_strings": [
            "random clifford",
            "almost clifford"
        ],
        "min_qubit_number" : 1,
        "max_qubit_number" : 9,
        "reps" : int(1e3) 
    }
]

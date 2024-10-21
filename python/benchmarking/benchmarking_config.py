"""We have the following benchmarks:

1. Testing S_V
2. S_V -> succinct representation
3. Succinct representation -> S_V
4. S_P -> succinct representation
5. Succinct representation -> S_P
6. S_V -> S_P
7. S_P -> S_V

8. Testing C_U
9. C_U -> succinct representation
10. Succinct representation -> C_U
"""

import generators as gs
import qiskit.quantum_info as qi
import stim

import sys
PATH_TO_LIBRARY = './build/ninja-multi-vcpkg/cpp/src/Release'
# PATH_TO_STIM_MOCK = './build/ninja-multi-vcpkg/cpp/benchmarking/stim/Release'
sys.path.extend([PATH_TO_LIBRARY])
import fast as fst # type: ignore
# import stim_mock as sm

def stim_S_V_test(matrix):
    try:
        stim.Tableau.from_state_vector(matrix, endian='big')
    except:
        pass

def our_check_matrix_to_succinct():
    pass

def stim_check_matrix_to_succinct():
    pass

def our_succinct_to_check_matrix():
    pass

def stim_succinct_to_check_matrix():
    pass

def our_statevector_to_check_matrix():
    pass

def stim_statevector_to_check_matrix():
    pass

def our_check_matrix_to_statevector():
    pass

def stim_check_matrix_to_statevector():
    pass

def qiskit_C1_converter(matrix, assume_valid = True):
    return qi.Clifford.from_matrix(matrix)

def stim_C1_convertor(matrix, assume_valid = True):
    return stim.Tableau.from_unitary_matrix(matrix, endian = "big")

def qiskit_C1_test(matrix):
    try:
        qi.Clifford.from_matrix(matrix)
    except:
        pass

def stim_C1_test(matrix):
    try:
        stim.Tableau.from_unitary_matrix(matrix, endian = "big")
    except:
        pass


configs = [
    {
        "pre_string": "Testing S_V",
        "functions_to_time": [
            stim_S_V_test,
            fst.is_stabiliser_state
        ],
        "function_strings": [
            "stim",
            "our method"
        ],
        "generation_types": [
            gs.rand_s_v_function
        ],
        "generation_strings": [
            "random stabiliser state or almost stab state"
        ],
        "min_qubit_number" : 1,
        "max_qubit_number" : 12,
        "reps" : int(1e3)
    },
    
    {
        "pre_string": "S_V -> succinct representation",
        "functions_to_time": [
            stim_S_V_test,
            fst.stabiliser_state_from_statevector
        ],
        "function_strings": [
            "stim",
            "our method"
        ],
        "generation_types": [
            gs.rand_s_v_to_succinct_function
        ],
        "generation_strings": [
            "rand_s_v_to_succinct_function"
        ],
        "min_qubit_number" : 1,
        "max_qubit_number" : 12,
        "reps" : int(1e3)
    },

    {
        "pre_string": "Succinct representation -> S_V",
        "functions_to_time": [
            stim_S_V_test,
            fst.stabiliser_state_from_statevector
        ],
        "function_strings": [
            "stim",
            "our method"
        ],
        "generation_types": [ # TODO
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
        "max_qubit_number" : 12,
        "reps" : int(1e3)
    },

    {
        "pre_string": "S_P -> succinct representation",
        "functions_to_time": [
            stim_check_matrix_to_succinct,
            our_check_matrix_to_succinct
        ],
        "function_strings": [
            "stim",
            "our method"
        ],
        "generation_types": [ # TODO
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
        "max_qubit_number" : 12,
        "reps" : int(1e3)
    },

    {
        "pre_string": "Succinct representation -> S_P",
        "functions_to_time": [
            stim_succinct_to_check_matrix,
            our_succinct_to_check_matrix
        ],
        "function_strings": [
            "stim",
            "our method"
        ],
        "generation_types": [ # TODO
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
        "max_qubit_number" : 12,
        "reps" : int(1e3)
    },

    {
        "pre_string": "S_V -> S_P",
        "functions_to_time": [
            stim_statevector_to_check_matrix,
            our_statevector_to_check_matrix
        ],
        "function_strings": [
            "stim",
            "our method"
        ],
        "generation_types": [ # TODO
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
        "max_qubit_number" : 12,
        "reps" : int(1e3)
    },

    {
        "pre_string": "S_P -> S_V",
        "functions_to_time": [
            stim_check_matrix_to_statevector,
            our_check_matrix_to_statevector
        ],
        "function_strings": [
            "stim",
            "our method"
        ],
        "generation_types": [ # TODO
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
        "max_qubit_number" : 12,
        "reps" : int(1e3)
    },

    {
        "pre_string": "Testing C_U",
        "functions_to_time":[
            fst.is_clifford_matrix,
        ],
        "function_strings": [
            "our method",
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
    },

    {
        "pre_string": "C_U -> succinct representation",
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
        "pre_string": "Succinct representation -> C_U",
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
        "generation_types": [ # TODO
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
]
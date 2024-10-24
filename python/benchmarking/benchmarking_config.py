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
import pickle
import stim

import sys
PATH_TO_LIBRARY = './build/ninja-multi-vcpkg/cpp/src/Release'
# PATH_TO_STIM_MOCK = './build/ninja-multi-vcpkg/cpp/benchmarking/stim/Release'
sys.path.extend([PATH_TO_LIBRARY])
import fast as fst # type: ignore
# import stim_mock as sm

def stim_S_V_test(statevector):
    try:
        stim.Tableau.from_state_vector(statevector, endian='big')
    except:
        pass

def stim_S_V_to_succinct(statevector, assume_valid=True):
    stim.Tableau.from_state_vector(statevector, endian='big')

def our_S_V_to_succinct(statevector, assume_valid=True):
    fst.stabiliser_state_from_statevector(statevector, assume_valid)

def stim_succinct_to_S_V(our_succinct: fst.Stabiliser_State, stim_succinct: stim.Tableau):
    stim_succinct.to_state_vector()

def our_succinct_to_S_V(our_succinct: fst.Stabiliser_State, stim_succinct: stim.Tableau):
    our_succinct.get_state_vector()

def our_check_matrix_to_succinct(check_matrix: fst.Check_Matrix):
    fst.Stabiliser_State(check_matrix)

def our_succinct_to_check_matrix(our_succinct: fst.Stabiliser_State):
    fst.Check_Matrix(our_succinct)

def stim_S_V_to_check_matrix(statevector):
    stim.Tableau.from_state_vector(statevector, endian='big').to_stabilizers()

def our_S_V_to_check_matrix(statevector):
    fst.Check_Matrix(fst.stabiliser_state_from_statevector(statevector))

def stim_check_matrix_to_statevector(our_check_matrix: fst.Check_Matrix, stim_check_matrix: stim.Tableau):
    stim_check_matrix.to_state_vector()

def our_check_matrix_to_statevector(our_check_matrix: fst.Check_Matrix, stim_check_matrix: stim.Tableau):
    fst.Stabiliser_State(our_check_matrix)

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

def our_succinct_to_C_U(our_clifford: fst.Clifford, qiskit_clifford: qi.Clifford, stim_clifford: stim.Tableau):
    our_clifford.get_matrix()

def qiskit_succinct_to_C_U(our_clifford: fst.Clifford, qiskit_clifford: qi.Clifford, stim_clifford: stim.Tableau):
    qiskit_clifford.to_matrix()

def stim_succinct_to_C_U(our_clifford: fst.Clifford, qiskit_clifford: qi.Clifford, stim_clifford: stim.Tableau):
    stim_clifford.to_unitary_matrix(endian='big')


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
            gs.rand_s_v
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
            stim_S_V_to_succinct,
            our_S_V_to_succinct
        ],
        "function_strings": [
            "stim",
            "our method"
        ],
        "generation_types": [
            gs.rand_s_v_to_succinct
        ],
        "generation_strings": [
            "rand_s_v_to_succinct"
        ],
        "min_qubit_number" : 1,
        "max_qubit_number" : 12,
        "reps" : int(1e3)
    },

    {
        "pre_string": "Succinct representation -> S_V",
        "functions_to_time": [
            stim_succinct_to_S_V,
            our_succinct_to_S_V
        ],
        "function_strings": [
            "stim",
            "our method"
        ],
        "generation_types": [
            gs.rand_succinct
        ],
        "generation_strings": [
            "rand_succinct"
        ],
        "min_qubit_number" : 1,
        "max_qubit_number" : 12,
        "reps" : int(1e3)
    },

    {
        "pre_string": "S_P -> succinct representation",
        "functions_to_time": [
            our_check_matrix_to_succinct
        ],
        "function_strings": [
            "our method"
        ],
        "generation_types": [
            gs.rand_our_check_matrix
        ],
        "generation_strings": [
            "rand_our_check_matrix"
        ],
        "min_qubit_number" : 1,
        "max_qubit_number" : 12,
        "reps" : int(1e3)
    },

    {
        "pre_string": "Succinct representation -> S_P",
        "functions_to_time": [
            our_succinct_to_check_matrix
        ],
        "function_strings": [
            "our method"
        ],
        "generation_types": [
            gs.rand_our_succinct
        ],
        "generation_strings": [
            "rand_our_succinct"
        ],
        "min_qubit_number" : 1,
        "max_qubit_number" : 12,
        "reps" : int(1e3)
    },

    {
        "pre_string": "S_V -> S_P",
        "functions_to_time": [
            stim_S_V_to_check_matrix,
            our_S_V_to_check_matrix
        ],
        "function_strings": [
            "stim",
            "our method"
        ],
        "generation_types": [
            gs.random_stab_state
        ],
        "generation_strings": [
            "random_stab_state"
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
        "generation_types": [
            gs.rand_check_matrix
        ],
        "generation_strings": [
            "rand_check_matrix"
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
            gs.rand_clifford_test
        ],
        "generation_strings": [
            "rand_clifford_test"
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
            gs.rand_clifford_matrix
        ],
        "generation_strings": [
            "rand_clifford_matrix"
        ],
        "min_qubit_number" : 1,
        "max_qubit_number" : 9,
        "reps" : int(1e3) 
    },

    {
        "pre_string": "Succinct representation -> C_U",
        "functions_to_time": [
            our_succinct_to_C_U,
            qiskit_succinct_to_C_U,
            stim_succinct_to_C_U
        ],
        "function_strings": [
            "our method",
            "Qiskit",
            "stim"
        ],
        "generation_types": [
            gs.rand_clifford_succinct
        ],
        "generation_strings": [
            "rand_clifford_succinct"
        ],
        "min_qubit_number" : 1,
        "max_qubit_number" : 9,
        "reps" : int(1e3) 
    },
]

# configs = configs[-1:]
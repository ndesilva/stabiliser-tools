"""We have the following benchmarks:

1. Testing S_V
2. S_V to succinct representation
3. Succinct representation to S_V
4. S_P to succinct representation
5. Succinct representation to S_P
6. S_V to S_P
7. S_P to S_V

8. Testing C_U
9. C_U to succinct representation
10. Succinct representation to C_U
"""

import generators as gs
import qiskit.quantum_info as qi
import pickle
import stim
# import stab_tools as fst

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
    fst.Check_Matrix(fst.stabiliser_state_from_statevector(statevector, assume_valid=True))

def stim_check_matrix_to_statevector(our_check_matrix: fst.Check_Matrix, stim_check_matrix: stim.Tableau):
    stim_check_matrix.to_state_vector()

def our_check_matrix_to_statevector(our_check_matrix: fst.Check_Matrix, stim_check_matrix: stim.Tableau):
    fst.Stabiliser_State(our_check_matrix)

def our_C_U_converter(matrix, assume_valid=True):
    return fst.clifford_from_matrix(matrix, assume_valid)

def qiskit_C_U_converter(matrix, assume_valid = True):
    return qi.Clifford.from_matrix(matrix)

def stim_C_U_converter(matrix, assume_valid = True):
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
        "title": r"Testing $[S_V]$",
        "functions_to_time": [
            fst.is_stabiliser_state,
            stim_S_V_test
        ],
        "function_strings": [
            "our method",
            "stim"
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
        "pre_string": "S_V to succinct representation",
        "title": r"$[S_V]$ to succinct rep",
        "functions_to_time": [
            our_S_V_to_succinct,
            stim_S_V_to_succinct
        ],
        "function_strings": [
            "our method",
            "stim"
        ],
        "generation_types": [
            gs.rand_s_v_to_succinct
        ],
        "generation_strings": [
            "rand_s_v_to_succinct"
        ],
        "min_qubit_number" : 3,
        "max_qubit_number" : 12,
        "reps" : int(1e3)
    },

    {
        "pre_string": "Succinct representation to S_V",
        "title": r"Succinct rep to $[S_V]$",
        "functions_to_time": [
            our_succinct_to_S_V,
            stim_succinct_to_S_V
        ],
        "function_strings": [
            "our method",
            "stim"
        ],
        "generation_types": [
            gs.rand_succinct
        ],
        "generation_strings": [
            "rand_succinct"
        ],
        "min_qubit_number" : 3,
        "max_qubit_number" : 12,
        "reps" : int(1e3)
    },

    {
        "pre_string": "S_P to succinct representation",
        "title": r"$[S_P]$ to succinct rep",
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
        "pre_string": "Succinct representation to S_P",
        "title": r"Succinct rep to $[S_P]$",
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
        "pre_string": "S_V to S_P",
        "title": r"$[S_V]$ to $[S_P]$",
        "functions_to_time": [
            our_S_V_to_check_matrix,
            stim_S_V_to_check_matrix
        ],
        "function_strings": [
            "our method",
            "stim"
        ],
        "generation_types": [
            gs.random_stab_state
        ],
        "generation_strings": [
            "random_stab_state"
        ],
        "min_qubit_number" : 3,
        "max_qubit_number" : 12,
        "reps" : int(1e3)
    },

    {
        "pre_string": "S_P to S_V",
        "title": r"$[S_P]$ to $[S_V]$",
        "functions_to_time": [
            our_check_matrix_to_statevector,
            stim_check_matrix_to_statevector
        ],
        "function_strings": [
            "our method",
            "stim"
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
        "title": r"Testing $[C_U]$",
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
            gs.rand_clifford_test
        ],
        "generation_strings": [
            "rand_clifford_test"
        ],
        "min_qubit_number" : 1,
        "max_qubit_number" : 10,
        "reps" : int(1e3) 
    },

    {
        "pre_string": "C_U to succinct representation",
        "title": r"$[C_U]$ to succinct rep",
        "functions_to_time": [
            our_C_U_converter,
            qiskit_C_U_converter,
            stim_C_U_converter
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
        "min_qubit_number" : 3,
        "max_qubit_number" : 10,
        "reps" : int(1e3) 
    },

    {
        "pre_string": "Succinct representation to C_U",
        "title": r"Succinct rep to $[C_U]$",
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
        "max_qubit_number" : 10,
        "reps" : int(1e3) 
    },
]

# configs = configs[0:3]
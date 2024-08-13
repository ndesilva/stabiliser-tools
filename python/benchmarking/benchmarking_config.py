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
]
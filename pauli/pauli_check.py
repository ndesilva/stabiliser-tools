import numpy as np
import numba
import F2_helper.F2_helper as f2

# Assuming matrix of size 2^n by 2^n, returns whether matrix is in the Pauli group.
def is_pauli(matrix : np.ndarray, allow_global_factor = False) -> bool:
    size = matrix.shape[0]
    n = f2.fast_log2(size)

    if np.count_nonzero(matrix) != size:
        #print('reject due to total non-zero')
        return False
    
    first_col_nonzero = np.nonzero(matrix[:,0])[0]

    if len(first_col_nonzero) != 1:
        #print('reject due to first col non-zero')
        return False

    phase = matrix[first_col_nonzero[0], 0]

    if not (allow_global_factor or is_valid_pauli_entry(phase)):
        #print('reject due to first col invalid entry')
        return False
    
    q = first_col_nonzero[0]
    p = 0

    for j in range(n):
        col = 1 << j

        entry = matrix[col ^ q, col]/phase

        match entry:
            case 1:
                pass
            case -1:
                p |= col #add a Z operator corresponding to this column
            case _:
                #print('reject due to q col invalid entry')
                return False
    
    for col in range(1, size): # we are repeating the q columns here, speed up?
        entry = matrix[col ^ q, col]
        value = phase*f2.sign_mod2product(p, col)

        if entry != value: # if on every loop is maybe not ideal; how to speed this up?
            #print('reject due to a remaining entry invalid')
            return False
        
    return True

def is_valid_pauli_entry(entry : float) -> bool:
    match entry:
        case 1:
            return True
        case -1:
            return True
        case 1j:
            return True
        case -1j:
            return True
        case _:
            return False

def is_pm_one(entry : float) -> bool:
    match entry:
        case 1:
            return True
        case -1:
            return True
        case _:
            return False
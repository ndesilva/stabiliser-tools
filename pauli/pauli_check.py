import numpy as np
import numba

# Assuming matrix of size 2^n by 2^n, returns whether matrix is in the Pauli group.
def is_pauli(matrix : np.ndarray, allow_global_factor = False) -> bool:
    size = matrix.shape[0]
    n = int(np.log2(size))

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
    
    p = first_col_nonzero[0]
    q = 0

    for j in range(n):
        col = 1 << j

        entry = matrix[col ^ p, col]/phase

        match entry:
            case 1:
                pass
            case -1:
                q |= col #add a Z operator corresponding to this column
            case _:
                #print('reject due to q col invalid entry')
                return False
    
    for col in range(1, size): # we are repeating the q columns here, speed up?
        entry = matrix[col ^ p, col]
        value = phase*phase_mod2product(q, col)

        if entry != value: # if on every loop is maybe not ideal; how to speed this up?
            #print('reject due to a remaining entry invalid')
            return False
        
    return True

# generates the n-qubit pauli matrix sX^pZ^q
def generate_pauli(n : int, s : complex, p : int, q : int) -> np.ndarray:
    size = 2**n
    matrix = np.zeros((size, size), dtype=complex)

    for j in range(size):
        matrix[j^p, j] = s*phase_mod2product(j, q)

    return matrix

# returns (-1) ** the mod 2 inner product of the binary representations of x, y
@numba.njit()
def phase_mod2product(x : int, y : int) -> int:
    product = x & y

    pairity = 0
    
    while product:
        pairity ^= 1
        product &= product -1

    return 1-2*pairity

@numba.njit()
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

@numba.njit()
def is_pm_one(entry : float) -> bool:
    match entry:
        case 1:
            return True
        case -1:
            return True
        case _:
            return False
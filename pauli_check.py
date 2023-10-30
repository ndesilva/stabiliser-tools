import numpy as np

# assuming matrix of size 2^n by 2^n
def is_square_matrix_pauli(matrix : np.ndarray) -> bool:
    size = matrix.shape[0]
    n = np.log2(size)

    if not correct_nonzero_entries(matrix, size):
        return False
    
    first_col_nonzero = np.nonzero(matrix[:,0])[0]

    if len(first_col_nonzero) != 1:
        return False

    phase = matrix[first_col_nonzero[0], 0]

    if not is_valid_pauli_entry(phase):
        return False
    
    p = first_col_nonzero

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
                return False
    
    for col in range(1, size): # we are repeating the q columns here, speed up?
        entry = matrix[col ^ p, col]
        value = phase*mod2product(q, col)

        if entry != value: # if on every loop is maybe not ideal; how to speed this up?
            return False
        
    return True

# Find the mod 2 inner product of the binary representations of x.y
def mod2product(x : int, y : int) -> int:
    prod = x & y

    wt = 0
    
    while prod:
        wt ^= 1
        prod &= prod -1

    return wt
        
def correct_nonzero_entries(matrix : np.ndarray, size : int) -> bool:
    return len(np.flatnonzero(matrix)[0]) == size

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
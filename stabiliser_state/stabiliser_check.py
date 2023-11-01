import numpy as np
import functools

# Assuming vector of length 2^n, returns whether vector is a stabiliser state. Currently assumes all entries are +-1, generalise to complex entries
def is_pauli(state_vector : np.ndarray, allow_global_factor = False) -> bool:    
    nonzero_indices = np.nonzero(state_vector)[0]
    support_size = len(nonzero_indices)
    
    k = np.log2(support_size)

    if not k.is_integer():
        return False
    
    k = int(k)
    n = int(np.log2(state_vector))

    support_indicies_bitstrings = np.array([[int(bit) for bit in f'{index:0{n}b}']
                    for index in nonzero_indices])
    
    shift = support_indicies_bitstrings[0]
    shifted_bitstrings = support_indicies_bitstrings ^ shift

    row_reduced = gaussian_elim_F2(shifted_bitstrings)
    non_zero_rows = row_reduced[~np.all(row_reduced == 0, axis=1)] # could this be made faster by just having gaussian_elim_F2 chop off the zero rows?

    if non_zero_rows.shape[0] != k:
        return False
    
    non_zero_coeffs = state_vector[nonzero_indices]
    phase = non_zero_coeffs[0]

    if not (allow_global_factor or is_valid_stabiliser_entry(phase)):
        #print('reject due to first col invalid entry')
        return False

    linear_real_part = []
    imag_part = []

    # get linear terms
    for j in range(k):
        index = 1 << j
        match state_vector[index]/phase:
            case 1:
                pass
            case -1:
                linear_real_part.append(index)
            case 1j:
                imag_part.append(index)
            case -1j:
                linear_real_part.append(index)
                imag_part.append(index)
            case _:
                return False
    
    quadratic_real_part = []

    # get quadratic terms:
    for a in range(k-1):
        for b in range(a+1, k):
            index = (1 << a) | (1 << b)
            
            linear_part = (1-2*evaluate_poly(linear_real_part, index))*(1+(1j-1)*evaluate_poly(imag_part, index))

            value = state_vector[index]/(phase*linear_part)

            match value:
                case 1:
                    pass
                case -1:
                    quadratic_real_part.append(index)
                case _:
                    return False
                
    real_part = linear_real_part + quadratic_real_part

    for index in range(2**k): # We are repeating columns of Hamming weight 1,2 - fast way to not do this?
        value = phase * (1-2*evaluate_poly(real_part, index)) * (1+(1j-1)*evaluate_poly(imag_part, index))

        if state_vector[index] != value:
            return False
        
    return True

def evaluate_poly(non_zero_coeffs : list[int], integer : int) -> int:
    b = [ integer & index == index for index in non_zero_coeffs]
    return int(functools.reduce(lambda x,y : x^y, b))


def is_valid_stabiliser_entry(entry : float) -> bool:
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


def gaussian_elim_F2(M : np.ndarray) -> np.ndarray:
    num_rows, num_cols = M.shape
    M = np.copy(M)

    row_to_comp = 0
    for j in range(num_cols):
        col = M[row_to_comp:, j]
        if np.count_nonzero(col) > 0:
            i = np.nonzero(col)[0][0] + row_to_comp
            M[(row_to_comp, i), :] = M[(i, row_to_comp), :]

            for ii in range(num_rows):
                if ii != row_to_comp and M[ii, j] != 0:
                    M[ii, :] ^= M[row_to_comp, :]

            row_to_comp += 1
            if row_to_comp == num_rows:
                break

    return M
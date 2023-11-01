import numpy as np
import functools
from operator import itemgetter

# Assuming vector of length 2^n, returns whether vector is a stabiliser state. Currently assumes all entries are +-1, generalise to complex entries
def is_pauli(state_vector : np.ndarray, allow_global_factor = False) -> bool:    
    
    nonzero_indices = np.nonzero(state_vector)[0]
    support_size = len(nonzero_indices)
    
    k = np.log2(support_size)

    # check support is a power of 2
    if not k.is_integer():
        return False
    
    k = int(k)
    n = int(np.log2(state_vector))

    shift = nonzero_indices[0]

    # shift affine space to a vector space, remembering the pairing of the indicies and entries in the state vector
    vector_space_value_pairs = [(index ^ shift, state_vector[index]) for index in nonzero_indices]
    
    # sort the vector space in increasing F_2 order; can then apply lemma on sorted vector spaces
    vector_space_value_pairs.sort(key = itemgetter(0)) # itemgetter is faster than declaring a lambda function lambda x : x[0] according to the internet

    vector_space_indicies = [pair[0] for pair in vector_space_value_pairs]
    
    weight_one_bitstrings = [1<<j for j in range(k)]
    basis_vectors = [vector_space_indicies[index] for index in weight_one_bitstrings]

    # Csing lemma, check that the indicies form an F2 vector space
    for j in range(2**k): # we check the basis vectors again here, fast way to not do that?
        vectors = [(j & weight_one_bitstrings[l] == weight_one_bitstrings[l])*basis_vectors[l] for l in range(k)]
        value = functools.reduce(lambda x,y : x^y, vectors)

        if vector_space_indicies[j] != value:
            return False
    
    non_zero_coeffs = [pair[1] for pair in vector_space_value_pairs]
    phase = non_zero_coeffs[0]

    if not (allow_global_factor or is_valid_stabiliser_entry(phase)):
        #print('reject due to first non-zero entry invalid entry')
        return False

    linear_real_part = []
    imag_part = []

    # get linear terms
    for index in weight_one_bitstrings:
        match non_zero_coeffs[index]/phase:
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

            value = non_zero_coeffs[index]/(phase*linear_part)

            match value:
                case 1:
                    pass
                case -1:
                    quadratic_real_part.append(index)
                case _:
                    return False
    
    # Combine real and quadratic terms in a quadratic form
    real_part = linear_real_part + quadratic_real_part

    for index in range(2**k): # We are repeating columns of Hamming weight 1,2 - fast way to not do this?
        value = phase * (1-2*evaluate_poly(real_part, index)) * (1+(1j-1)*evaluate_poly(imag_part, index))

        if non_zero_coeffs[index] != value:
            return False
        
    return True

def evaluate_poly(non_zero_coeffs : list[int], integer : int) -> int:
    terms = [ integer & index == index for index in non_zero_coeffs]
    return int(functools.reduce(lambda x,y : x^y, terms))

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
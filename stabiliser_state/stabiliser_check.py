import numpy as np
import functools
from operator import itemgetter
from F2_helper.F2_helper import sign_mod2product, imag_mod2product, sign_evaluate_poly
from stabiliser_state.Stabiliser_State import Stabiliser_State

# Assuming vector of length 2^n, returns whether vector is a stabiliser state. Currently assumes all entries are +-1, generalise to complex entries
def is_stabiliser_state(state_vector : np.ndarray, allow_global_factor = False) -> bool | Stabiliser_State:    
    nonzero_indices = np.nonzero(state_vector)[0]
    support_size = len(nonzero_indices)
    
    dimension = np.log2(support_size)

    # check support is a power of 2
    if not dimension.is_integer():
        #print('support not power of 2')
        return False
    dimension = int(dimension)

    shift = nonzero_indices[0]

    # shift affine space to a vector space, remembering the pairing of the indicies and entries in the state vector
    vector_space_value_pairs = [(index ^ shift, state_vector[index]) for index in nonzero_indices]
    
    # sort the vector space in increasing F_2 order; can then apply lemma on sorted vector spaces
    vector_space_value_pairs.sort(key = itemgetter(0)) # itemgetter is faster than declaring a lambda function lambda x : x[0] according to the internet

    vector_space_indicies = [pair[0] for pair in vector_space_value_pairs]
    
    weight_one_bitstrings = [1<<j for j in range(dimension)]
    basis_vectors = [vector_space_indicies[index] for index in weight_one_bitstrings]

    # Csing lemma, check that the indicies form an F2 vector space
    for j in range(1<<dimension): # we check the basis vectors again here, fast way to not do that?
        vectors = [(j & weight_one_bitstrings[l] == weight_one_bitstrings[l])*basis_vectors[l] for l in range(dimension)]
        value = functools.reduce(lambda x,y : x^y, vectors, 0)

        if vector_space_indicies[j] != value:
            #print('Support not affine space')
            return False
    
    non_zero_coeffs = [pair[1] for pair in vector_space_value_pairs]
    first_entry = non_zero_coeffs[0]

    if not (allow_global_factor or is_valid_stabiliser_entry(first_entry*(np.sqrt(support_size)))):
        #print('invalid first entry')
        return False

    linear_real_part = 0
    imag_part = 0

    # get linear terms
    for index in weight_one_bitstrings:
        match non_zero_coeffs[index]/first_entry:
            case 1:
                pass
            case -1:
                linear_real_part |= index
            case 1j:
                imag_part |= index
            case -1j:
                linear_real_part |= index
                imag_part |= index
            case _:
                #print('invalid linear term')
                return False
    
    quadratic_real_part = []

    # get quadratic terms:
    for a in range(dimension-1):
        for b in range(a+1, dimension):
            index = (1 << a) | (1 << b)
            
            linear_part = sign_mod2product(index, linear_real_part)*imag_mod2product(index, imag_part)

            value = non_zero_coeffs[index]/(first_entry*linear_part)

            match value:
                case 1:
                    pass
                case -1:
                    quadratic_real_part.append(index)
                case _:
                    #print('invalid quadratic term')
                    return False

    for index in range(1<<dimension): # We are repeating columns of Hamming weight 1,2 - fast way to not do this?
        value = first_entry*imag_mod2product(index, imag_part)*sign_mod2product(index, linear_real_part)*sign_evaluate_poly(quadratic_real_part, index)

        if non_zero_coeffs[index] != value:
            #print('inconsistent remainder')
            return False
        
    #print('State accepted \n')
    return True

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
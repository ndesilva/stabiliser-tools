import functools
import numba

def evaluate_poly(non_zero_coeffs : list[int], integer : int) -> int:
    terms = [ integer & index == index for index in non_zero_coeffs]
    return int(functools.reduce(lambda x,y : x^y, terms, 0))

def sign_evaluate_poly(non_zero_coeffs : list[int], integer : int) -> int:
    return 1-2*evaluate_poly(non_zero_coeffs, integer)

@numba.njit()
def mod2product(x : int, y : int) -> int:
    product = x & y

    wt = 0

    while product:
        product &= product -1
        wt ^= 1

    return wt

@numba.njit()
def sign_mod2product(x : int, y : int) -> int:
    return 1-2*mod2product(x,y)

@numba.njit()
def imag_mod2product(x : int, y : int) -> complex:
    return 1 + (1j-1)*mod2product(x,y)

def fast_log2(x : int) -> int:
    return x.bit_length() - 1

def get_bit_at(x : int, bit_location: int) -> int:
    power_of_two = 1 << bit_location

    return (x & power_of_two) == power_of_two

def get_vector_expansion(dimension : int, basis_vectors : list[int], coefficients : int) -> int:
    terms = [get_bit_at(coefficients, j)*basis_vectors[j] for j in range(dimension)]
    return functools.reduce(lambda x,y : x ^ y, terms, 0)
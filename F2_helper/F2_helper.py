import functools
import numpy as np

def evaluate_poly(non_zero_coeffs : list[int], integer : int) -> int:
    terms = [ integer & index == index for index in non_zero_coeffs]
    return int(functools.reduce(lambda x,y : x^y, terms, 0))

def sign_evaluate_poly(non_zero_coeffs : list[int], integer : int) -> int:
    return 1-2*evaluate_poly(non_zero_coeffs, integer)

def mod2product(x : int, y : int) -> int:
    product = x & y

    pairity = 0
    
    while product:
        pairity ^= 1
        product &= product -1

    return pairity

def sign_mod2product(x : int, y : int) -> int:
    return 1-2*mod2product(x,y)

def imag_mod2product(x : int, y : int) -> complex:
    return 1 + (1j-1)*mod2product(x,y)

def num_bin_digits(x : int) -> int:
    return len(bin(x)[2:])

def fast_log2(x : int) -> int:
    return num_bin_digits(x) - 1

def get_bit_at(x : int, bit_location: int) -> int:
    power_of_two = 1 << bit_location

    return (x & power_of_two) == power_of_two
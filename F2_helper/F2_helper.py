import functools

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
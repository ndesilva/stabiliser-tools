import functools

def evaluate_poly(non_zero_coeffs : list[int], integer : int) -> int:
    terms = [ integer & index == index for index in non_zero_coeffs]
    return int(functools.reduce(lambda x,y : x^y, terms))

def phase_mod2product(x : int, y : int) -> int:
        product = x & y

        pairity = 0
        
        while product:
            pairity ^= 1
            product &= product -1

        return 1-2*pairity
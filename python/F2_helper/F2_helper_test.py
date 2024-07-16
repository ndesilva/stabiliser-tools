import unittest
from F2_helper.F2_helper import *

class Test_Pauli_Class(unittest.TestCase): 
    def test_mod2product(self):
        x = 26 #binary 11010
        y = 19 #binary 10011
        z = 3  #binary 00011

        self.assertEqual(mod2product(x,y), 0)
        self.assertEqual(mod2product(x,z), 1)

    def test_sign_mod2product(self):
        x = 51 #binary 110011
        y = 3  #binary 000011
        z = 49 #binary 110101

        self.assertEqual(sign_mod2product(x,y), 1)
        self.assertEqual(sign_mod2product(x,z), -1)

    def test_imag_mod2product(self):
        x = 43 #binary 101011
        y = 3  #binary 000011
        z = 45 #binary 101101

        self.assertEqual(imag_mod2product(x,y), 1)
        self.assertEqual(imag_mod2product(x,z), 1j)

    def test_evaluate_poly_case_one(self):
        poly = [1,2,7] # x_1 + x_2 + x_1 x_2 x_3

        point_1 = 3 # binary 011
        point_2 = 7 # binary 111

        self.assertEqual(evaluate_poly(poly, point_1), 0)
        self.assertEqual(evaluate_poly(poly, point_2), 1)

    def test_evaluate_poly_case_two(self):
        poly = [3, 12] # x_1 x_2 + x_3 x_4

        point_1 = 15 # binary 1111
        point_2 = 7  # binary 0111

        self.assertEqual(evaluate_poly(poly, point_1), 0)
        self.assertEqual(evaluate_poly(poly, point_2), 1)

    def test_sign_evaluate_poly(self):
        poly = [3] #x_1 x_2
        
        point_1 = 1 # 01
        point_2 = 3 # 11

        self.assertEqual(sign_evaluate_poly(poly, point_1), 1)
        self.assertEqual(sign_evaluate_poly(poly, point_2), -1)

    def test_fast_log2_non_power_2(self):
        x = 9
        self.assertEqual(fast_log2(x), 3)

    def test_fast_log2_power_2(self):
        x = 16
        self.assertEqual(fast_log2(x), 4)

    def test_fast_log2_on_0(self):
        x = 0
        self.assertEqual(fast_log2(x), -1)

    def test_get_bit_at_returns_0(self):
        x = 54 # 110110
        
        self.assertEqual(get_bit_at(x,0), False)
        self.assertEqual(get_bit_at(x, 3), False)
    
    def test_get_bit_at_returns_1(self):
        x = 5 # 101

        self.assertEqual(get_bit_at(x, 0), True)
        self.assertEqual(get_bit_at(x, 2), True)

    def tets_get_bit_at_returns_0_when_index_larger_than_bit_rep(self):
        x = 7 #111

        self.assertEqual(get_bit_at(x, 3), False)

    def test_get_vector_expansion(self):
        dim = 3
        basis_vectors = [1, 2, 4] #001, 010, 100
        coeffs = 5 #101

        expected_sum = 5 # 101

        self.assertEqual(get_vector_expansion(dim, basis_vectors, coeffs), expected_sum)

    def test_get_vector_expansion_case_two(self):
        dim = 3
        basis_vectors = [17, 12, 26] #10001, 01100, 11010 
        coeffs = 7 #111
        
        expected_sum = 7 #00111


        self.assertEqual(get_vector_expansion(dim, basis_vectors, coeffs), expected_sum)

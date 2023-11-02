import unittest
from .F2_helper import phase_mod2product, evaluate_poly

class Test_Pauli_Class(unittest.TestCase): 
    def test_phase_mod2product_returns_1(self):
        x = 26 #binary 11010
        y = 19 #binary 10011

        self.assertEqual(phase_mod2product(x,y), 1)

    def test_phase_mod2product_returns_minus_1(self):
        x = 51 #binary 110011
        y = 49 #binary 110101

        self.assertEqual(phase_mod2product(x,y), -1)

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
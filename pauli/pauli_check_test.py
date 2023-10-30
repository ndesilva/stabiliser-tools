import pauli_check
import unittest
import numpy as np

class Test_Pauli_Check(unittest.TestCase):

    def test_is_pauli_accepts_single_qubit_pauli(self):
        n = 1
        s = -1j
        p = 1
        q = 1

        pauli = pauli_check.generate_pauli(n, s, p, q)

        self.assertTrue(pauli_check.is_pauli(pauli))

    def test_is_pauli_accepts_five_qubit_pauli(self):
        n = 5
        s = 1j
        p = 28 # binary 10110
        q = 3  # binary 00011

        pauli = pauli_check.generate_pauli(n,s,p,q)
        self.assertTrue(pauli_check.is_pauli(pauli))

    def test_is_pauli_rejects_pauli_with_extra_entry(self):
        n = 3
        s = 1
        p = 4 # binary 100
        q = 1 # binary 001

        corrupt_pauli = pauli_check.generate_pauli(n,s,p,q)

        corrupt_pauli[1,1] = 1 # add an extra 1

        self.assertFalse(pauli_check.is_pauli(corrupt_pauli))

    def test_is_pauli_rejects_pauli_with_first_column_all_zeros(self):
        n = 3
        s = 1
        p = 4 # binary 100
        q = 1 # binary 001

        corrupt_pauli = pauli_check.generate_pauli(n,s,p,q)

        corrupt_pauli[1,1] = 1 # add an extra 1
        corrupt_pauli[4,0] = 0 # remove only element in first column

        self.assertFalse(pauli_check.is_pauli(corrupt_pauli))

    def test_is_pauli_rejects_pauli_with_first_column_invalid_entry(self):
        n = 3
        s = 1
        p = 4 # binary 100
        q = 1 # binary 001

        corrupt_pauli = pauli_check.generate_pauli(n,s,p,q)

        corrupt_pauli[4,0] = 1 + 1j # modify element in first column to be invalid 

        self.assertFalse(pauli_check.is_pauli(corrupt_pauli))

    def test_is_pauli_rejects_pauli_with_q_column_invalid_entry(self):
        n = 3
        s = 1
        p = 4 # binary 100
        q = 1 # binary 001

        corrupt_pauli = pauli_check.generate_pauli(n,s,p,q)

        corrupt_pauli[5,1] = 1 + 1j # modify element in second column (q column) to be invalid

        self.assertFalse(pauli_check.is_pauli(corrupt_pauli))

    def test_is_pauli_rejects_pauli_with_remaining_column_invalid_entry(self):
        n = 3
        s = -1
        p = 4 # binary 100
        q = 1 # binary 001

        corrupt_pauli = pauli_check.generate_pauli(n,s,p,q)

        corrupt_pauli[1,5] = -1 # modify element in 5th column (remaining column) to be incorrect

        self.assertFalse(pauli_check.is_pauli(corrupt_pauli))

    def test_generate_pauli_single_qubit(self):
        n = 1
        s = 1j
        p = 1
        q = 1

        Y = np.array([[0, -1j],[1j, 0]]) # type: ignore

        self.assertTrue(np.array_equal(pauli_check.generate_pauli(n,s,p,q), Y))

    def test_generate_pauli_three_qubits(self):
        n = 3
        s = -1
        p = 6 # binary 110
        q = 5 # binary 101

        pauli = np.array([
            [0  ,0  ,0  ,0  ,0  ,0  ,1  ,0],
            [0  ,0  ,0  ,0  ,0  ,0  ,0  ,-1],
            [0  ,0  ,0  ,0  ,1  ,0  ,0  ,0],
            [0  ,0  ,0  ,0  ,0  ,-1 ,0  ,0],
            [0  ,0  ,-1 ,0  ,0  ,0  ,0  ,0],
            [0  ,0  ,0  ,1  ,0  ,0  ,0  ,0],
            [-1 ,0  ,0  ,0  ,0  ,0  ,0  ,0],
            [0  ,1  ,0  ,0  ,0  ,0  ,0  ,0]
        ])

        self.assertTrue(np.array_equal(pauli_check.generate_pauli(n,s,p,q), pauli))  
    
    def test_phase_mod2product_returns_1(self):
        x = 26 #binary 11010
        y = 19 #binary 10011

        self.assertEqual(pauli_check.phase_mod2product(x,y), 1)

    def test_phase_mod2product_returns_minus_1(self):
        x = 51 #binary 110011
        y = 49 #binary 110101

        self.assertEqual(pauli_check.phase_mod2product(x,y), -1)

if __name__ == '__main__':
    unittest.main()        
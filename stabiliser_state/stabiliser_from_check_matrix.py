import stabiliser_state.Check_Matrix as cm
import stabiliser_state.Stabiliser_State as ss
import F2_helper.F2_helper as f2

class Stabiliser_from_check_matrix:
    def __init__(self, check_matrix : cm.Check_Matrix):        
        self.number_qubits = check_matrix.number_qubits

        self.vector_basis = [pauli.x_vector for pauli in check_matrix.non_zero_x]
        self.dimension = len(self.vector_basis)

        self.__get_shift_vector(check_matrix)

        self.imag_part = 0
        self.linear_real_part = 0
        self.quadratic_form = []

        for j in range(self.dimension):
            v_j = self.vector_basis[j]
            beta_j = check_matrix.non_zero_x[j].z_vector
            imag_bit = check_matrix.non_zero_x[j].i_bit

            self.imag_part |= (1 << j) * imag_bit
            self.linear_real_part |= (1 << j) * ( check_matrix.non_zero_x[j].sign_bit ^ f2.mod2product(beta_j, v_j ^ self.shift) )

            for i in range(j):
                v_i = self.vector_basis[i]
                other_imag_bit = check_matrix.non_zero_x[i].sign_bit

                if f2.mod2product(beta_j, v_i) ^ other_imag_bit*imag_bit:
                    self.quadratic_form.append( 1 << i | 1 << j)

    def get_stab_state(self) -> ss.Stabiliser_State:
        return ss.Stabiliser_State(self.number_qubits, self.quadratic_form, self.linear_real_part, self.imag_part, self.vector_basis, self.shift, row_reduced = True)

    def __get_shift_vector(self, check_matrix):
        self.shift = 0

        for z_pauli in check_matrix.zero_x:
            pivot_index = f2.fast_log2(z_pauli.z_vector)
            self.shift |= (1<<pivot_index) * z_pauli.sign_bit
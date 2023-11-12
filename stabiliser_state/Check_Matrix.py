from __future__ import annotations
import stabiliser_state.Stabiliser_State as ss
import stabiliser_state.stabiliser_from_state_vector as ssv

import numpy as np
import pauli.Pauli as p
import F2_helper.F2_helper as f2

class Check_Matrix():

    @staticmethod
    def from_stabiliser_state(stab_state : ss.Stabiliser_State) -> Check_Matrix: # TODO test
        return stab_state.get_check_matrix()
    
    @staticmethod
    def from_statevector(state_vector : np.ndarray, assume_stab_state : bool = False) -> Check_Matrix: # TODO test
        state = ssv.Stabiliser_From_State_Vector(state_vector, allow_global_factor = True, assume_stab_state = assume_stab_state)

        if not state.is_stab_state:
            raise ValueError('State vector does not describe a stabiliser state')

        return Check_Matrix.from_stabiliser_state(state.get_stab_state())

    def __init__(self, paulis : list[p.Pauli], reduced_form : bool = False):
        self.paulis = paulis
        self.reduced_form = reduced_form
        self.number_qubits = len(paulis)

        self.non_zero_x : list[p.Pauli] = []
        self.zero_x : list[p.Pauli] = []
        
        self.__extract_zero_x_paulis()

    def get_stabiliser_state(self) -> ss.Stabiliser_State: # TODO test
        self.__put_into_reduced_form()
    
        vector_basis = [pauli.x_vector for pauli in self.non_zero_x]
        dimension = len(vector_basis)

        shift_vector = self.__get_shift_vector()

        imag_part = 0
        linear_real_part = 0
        quadratic_form = []

        imag_part, linear_real_part = self.__set_linear_and_quadratic_forms(vector_basis, dimension, shift_vector, imag_part, linear_real_part, quadratic_form)

        return ss.Stabiliser_State(self.number_qubits, quadratic_form, linear_real_part, imag_part, vector_basis, shift_vector, row_reduced = True)

    def __put_into_reduced_form(self) -> None:
        if self.reduced_form:
            return

        self.__row_reduce_non_zero_x()
        self.__row_reduce_zero_x()

        self.reduced_form = True

    def __set_linear_and_quadratic_forms(self, vector_basis : list[int], dimension : int, shift_vector : int, imag_part : int, linear_real_part : int, quadratic_form : list[int]) -> tuple[int, int]: # TODO test
        for j in range(dimension):
            v_j = vector_basis[j]
            beta_j = self.non_zero_x[j].z_vector
            imag_bit = self.non_zero_x[j].i_bit

            imag_part |= (1 << j) * imag_bit
            linear_real_part |= (1 << j) * ( self.non_zero_x[j].sign_bit ^ f2.mod2product(beta_j, v_j ^ shift_vector) )

            for i in range(j):
                v_i = vector_basis[i]
                other_imag_bit = self.non_zero_x[i].i_bit

                if f2.mod2product(beta_j, v_i) ^ other_imag_bit*imag_bit:
                    quadratic_form.append( 1 << i | 1 << j )

        return imag_part, linear_real_part
    
    def __get_shift_vector(self) -> int: # TODO test
        shift = 0

        for z_pauli in self.zero_x:
            pivot_index = f2.fast_log2(z_pauli.z_vector)
            shift |= (1 << pivot_index) * z_pauli.sign_bit

        return shift

    def __row_reduce_zero_x(self): # TODO test
        for pauli in self.zero_x:
            pivot_index = f2.fast_log2(pauli.z_vector)

            for other_pauli in self.zero_x:
                if other_pauli != pauli and f2.get_bit_at(other_pauli.z_vector, pivot_index):
                    other_pauli.multiply_by_pauli_on_right(pauli)

    def __row_reduce_non_zero_x(self): # TODO test
        to_remove = []
        
        for pauli in self.non_zero_x:
            pivot_index = f2.fast_log2(pauli.x_vector)

            if pivot_index == -1:
                self.zero_x.append(pauli)
                to_remove.append(pauli)
            
            else:
                for other_pauli in self.non_zero_x:
                    if other_pauli != pauli and f2.get_bit_at(other_pauli.x_vector, pivot_index):
                        other_pauli.multiply_by_pauli_on_right(pauli)
        
        for pauli in to_remove:
            self.non_zero_x.remove(pauli)

    def __extract_zero_x_paulis(self) -> None:
        for pauli in self.paulis:
            if f2.fast_log2(pauli.x_vector) == -1:
                self.zero_x.append(pauli)
            else:
                self.non_zero_x.append(pauli)
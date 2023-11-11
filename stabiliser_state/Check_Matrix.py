import pauli.Pauli as p
import F2_helper.F2_helper as f2

class Check_Matrix():
    def __init__(self, paulis : list[p.Pauli], assume_valid : bool = False, reduced_form : bool = False):
        self.paulis = paulis
        self.number_qubits = len(paulis)

        self.non_zero_x : list[p.Pauli] = []
        self.zero_x : list[p.Pauli] = []
        
        self.__extract_zero_x_paulis()

        if not reduced_form:
            self.__put_into_reduced_form()
        
    def __put_into_reduced_form(self) -> None: # TODO test
        if self.reduced_form:
            return

        self.__row_reduce_non_zero_x()
        self.__row_reduce_zero_x()

        self.reduced_form = True

    def __row_reduce_zero_x(self):
        for pauli in self.zero_x:
            pivot_index = f2.fast_log2(pauli.z_vector)

            for other_pauli in self.zero_x:
                if other_pauli != pauli and f2.get_bit_at(other_pauli.z_vector, pivot_index):
                    other_pauli.multiply_by_pauli_on_right(pauli)

    def __row_reduce_non_zero_x(self):
        for pauli in self.non_zero_x:
            pivot_index = f2.fast_log2(pauli.x_vector)

            if pivot_index == -1:
                self.zero_x.append(pauli)
                self.non_zero_x.remove(pauli)
            
            else:
                for other_pauli in self.non_zero_x:
                    if other_pauli != pauli and f2.get_bit_at(other_pauli.x_vector, pivot_index):
                        other_pauli.multiply_by_pauli_on_right(pauli)

    def __extract_zero_x_paulis(self) -> None: # TODO test
        for pauli in self.paulis:
            if f2.fast_log2(pauli.x_vector) == -1:
                self.zero_x.append(pauli)
            else:
                self.non_zero_x.append(pauli)
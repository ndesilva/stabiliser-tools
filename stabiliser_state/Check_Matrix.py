import pauli.Pauli as p

class Check_Matrix():
    def __init__(self, paulis : list[p.Pauli], assume_valid = False, reduced_form = False):
        self.paulis = paulis

        if not assume_valid:
            #self.__check_valid()
            pass
        
        self.reduced_form = reduced_form

    def put_into_reduced_form(self):
        if self.reduced_form:
            return
        
        pass

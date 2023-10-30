import time
import random

pauli_entries = [1, -1, 1j, -1j]

def list_check(float):
    return float in pauli_entries

def match_check(float):
    match float:
        case 1:
            return True
        case -1:
            return True
        case 1j:
            return True
        case _:
            return False
        
def or_check(float):
    if float == 1 or float == -1 or float == 1j or float == -1j:
        return True
    else:
        return False
    
functions_to_time = [list_check, match_check, or_check]
reps = int(1e7)

for function in functions_to_time:

    timer = 0

    for i in range(reps):
        r = random.choice([1, -1, 1j, -1j, 0.9, -0.9, 0.77, 0.005, 0, 0.1, -0.1j])
       
        st = time.time()       
        
        function(r)
        
        et = time.time()
        timer += et - st

    print(timer/reps)
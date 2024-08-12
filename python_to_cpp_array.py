import numpy as np
array = np.array([[ 0. +0.j ,  0. -0.5j,  0. +0.j ,  0. +0.5j,  0.5+0.j ,  0. +0.j ,
        -0.5+0.j ,  0. +0.j ],
       [ 0. +0.j ,  0. -0.5j,  0. +0.j ,  0. +0.5j, -0.5+0.j ,  0. +0.j ,
         0.5+0.j ,  0. +0.j ],
       [ 0.5+0.j ,  0. +0.j , -0.5+0.j ,  0. +0.j ,  0. +0.j ,  0. -0.5j,
         0. +0.j ,  0. +0.5j],
       [-0.5+0.j ,  0. +0.j ,  0.5+0.j ,  0. +0.j ,  0. +0.j ,  0. -0.5j,
         0. +0.j ,  0. +0.5j],
       [ 0. +0.j ,  0.5+0.j ,  0. +0.j ,  0.5+0.j ,  0. +0.5j,  0. +0.j ,
         0. +0.5j,  0. +0.j ],
       [ 0. +0.j , -0.5+0.j ,  0. +0.j , -0.5+0.j ,  0. +0.5j,  0. +0.j ,
         0. +0.5j,  0. +0.j ],
       [ 0. +0.5j,  0. +0.j ,  0. +0.5j,  0. +0.j ,  0. +0.j ,  0.5+0.j ,
         0. +0.j ,  0.5+0.j ],
       [ 0. +0.5j,  0. +0.j ,  0. +0.5j,  0. +0.j ,  0. +0.j , -0.5+0.j ,
         0. +0.j , -0.5+0.j ]])

array_string = '{'

for row in array:
    array_string += '{'

    for elt in row:
        if elt == 0:
            array_string += '.0f, '
        elif elt.imag == 0:
            array_string += f'{elt.real}f, '
        else:
            array_string += f'{elt.imag}f*I, '
    
    array_string = array_string[:-2]
    array_string += '},\n'

array_string = array_string[:-2]
array_string += '};'

print(array_string)
import sys

PATH_TO_LIBRARY = './build/ninja-multi-vcpkg/cpp/src/Release'

sys.path.append(PATH_TO_LIBRARY)

import fast

X = fast.Pauli(1,1,0,0,0)
print(X.x_vector)
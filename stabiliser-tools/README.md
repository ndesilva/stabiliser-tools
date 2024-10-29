# stabiliser-tools: Fast checks and conversions for stabiliser states, Cliffords and Paulis
**[GitHub project page](https://github.com/ndesilva/stabiliser-tools)**

**[Paper](https://arxiv.org/abs/2311.10357)**

`stabiliser-tools` supports the following specifications:
- **Stabiliser states:**
    - **[S<sub>V</sub>]**: a complex vector of amplitudes;
    - **[S<sub>P</sub>]**: a quadratic form, a linear map, and an affine subspace of $\mathbb{Z}_2^n$;
    - **[S<sub>Q</sub>]**: a check matrix, i.e. a compact list of Pauli gate generators for the stabiliser group.
- **Cliffords:**
    - **[C<sub>U</sub>]**: a unitary matrix;
    - **[C<sub>T</sub>]**: a list of $2n$ Pauli gates representing the images of basic Pauli gates under conjugation, i.e. a *tableau*.

After installing the package, its classes and functions can be accessed via
```
import stab_tools
```
Note that on Linux and other Unix-like operating systems, you may receive the following output upon importing:
```
WARNING: You may be about to receive an ImportError. (If you don't, you can safely ignore this message.) To resolve this, try running the following command in your terminal:

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<path-to-python-environment>/lib/python3.10/site-packages/stab_tools
```
This is due to the way shared libraries are linked in Linux. In order to stop this error from occurring in every new terminal session, consider updating `LD_LIBRARY_PATH` in your `.bashrc` or `.bashprofile` file.

After successfully importing `stab_tools`, run `help(stab_tools)` to view the list of functions and classes.

Available classes:
```
Clifford
Pauli
Check_Matrix
Stabiliser_State
```

Example code:
```
>>> import stab_tools
>>> import numpy as np
>>> v = np.array([0, 1, 0, 0, 0, 0, 1, 0]) / np.sqrt(2)
>>> stab_tools.is_stabiliser_state(v)
True
>>> s = stab_tools.stabiliser_state_from_statevector(v)
>>> paulis = stab_tools.Check_Matrix(s).get_paulis()
>>> paulis[0].get_matrix()
[[0j, 0j, 0j, 0j, 0j, 0j, 0j, (1-0j)], [0j, 0j, 0j, 0j, 0j, 0j, (1-0j), 0j], [0j, 0j, 0j, 0j, 0j, (1-0j), 0j, 0j], [0j, 0j, 0j, 0j, (1-0j), 0j, 0j, 0j], [0j, 0j, 0j, (1-0j), 0j, 0j, 0j, 0j], [0j, 0j, (1-0j), 0j, 0j, 0j, 0j, 0j], [0j, (1-0j), 0j, 0j, 0j, 0j, 0j, 0j], [(1-0j), 0j, 0j, 0j, 0j, 0j, 0j, 0j]]
>>> pauli_xs = [stab_tools.Pauli(3, 2**n, 0, 0, 0) for n in range(3)]
>>> pauli_xs[0].get_matrix()
[[0j, (1-0j), 0j, 0j, 0j, 0j, 0j, 0j], [(1-0j), 0j, 0j, 0j, 0j, 0j, 0j, 0j], [0j, 0j, 0j, (1-0j), 0j, 0j, 0j, 0j], [0j, 0j, (1-0j), 0j, 0j, 0j, 0j, 0j], [0j, 0j, 0j, 0j, 0j, (1-0j), 0j, 0j], [0j, 0j, 0j, 0j, (1-0j), 0j, 0j, 0j], [0j, 0j, 0j, 0j, 0j, 0j, 0j, (1-0j)], [0j, 0j, 0j, 0j, 0j, 0j, (1-0j), 0j]]
>>> H = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
>>> H2 = np.kron(H, H)
>>> H2
array([[ 0.5,  0.5,  0.5,  0.5],
       [ 0.5, -0.5,  0.5, -0.5],
       [ 0.5,  0.5, -0.5, -0.5],
       [ 0.5, -0.5, -0.5,  0.5]])
>>> stab_tools.is_clifford_matrix(H2)
True
```

## Compatibility
We have compiled and built the underlying C++ source code for a variety of different platforms and CPU architectures. However, if your machine is incompatible with any of the Python wheels, the full source code can be found in the source distribution, as well as on [GitHub](https://github.com/ndesilva/stabiliser-tools). You will need the libraries [Catch2](https://github.com/catchorg/Catch2) and [pybind11](https://github.com/pybind/pybind11).

If you are unable to build the code yourself, or if you have any other questions, contact the developers at ming_yin_2@sfu.ca.
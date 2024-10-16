# `stabiliser-tools`: Fast checks and conversions for stabiliser states, Cliffords and Paulis
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
Note that on Linux/UNIX, you may receive the following output upon importing:
```
WARNING: You may be about to receive an ImportError. (If you don't, you can safely ignore this message.) To resolve this, try running the following command in your terminal:

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<path-to-python-environment>/lib/python3.10/site-packages/stab_tools
```
This is due to the way shared libraries are linked in Linux/UNIX. In order to stop this error from occurring in every new terminal session, consider updating `LD_LIBRARY_PATH` in your `.bashrc` or `.bashprofile` file.

After successfully importing `stab_tools`, run `help(stab_tools)` to view the list of functions and classes.
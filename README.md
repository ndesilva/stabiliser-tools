# `stabiliser-tools`: Fast checks for stabiliser states, Cliffords and Paulis

## Installing `stabiliser-tools`
[...]

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

"""Fast checks and conversions for stabiliser states, Cliffords and Paulis.

Available classes (type 'help(stab_tools.<class>)' to see accessible functions and variables):
    Clifford
    Pauli
    Check_Matrix
    Stabiliser_State

Examples
--------
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
"""

import os
import sys

PATH_TO_LIBRARY = os.path.abspath(os.path.dirname(__file__))
ENV = ''
# print(f'Adding {PATH_TO_LIBRARY} to path if necessary')
if PATH_TO_LIBRARY not in sys.path:
    sys.path.extend([PATH_TO_LIBRARY])

if sys.platform in ['linux', 'darwin']:  # TODO Check
    try:
        match sys.platform:
            case 'linux':
                ENV = os.environ['LD_LIBRARY_PATH']
            case 'darwin':
                ENV = os.environ['LD_LIBRARY_PATH'] + os.environ['DYLD_LIBRARY_PATH'] + os.environ['LD_LIBRARY_PATH']
            case _:
                ENV = ''
        if PATH_TO_LIBRARY not in ENV:
            raise ImportError
    except (KeyError, ImportError):
        print('\nWARNING: You may be about to receive an ImportError. (If you don\'t, you can safely ignore this message.) '
              'To resolve this, try running the following command in your terminal:\n'
              f'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:{PATH_TO_LIBRARY}\n\n'
              'In order to stop this error from occurring in every new terminal session, consider updating LD_LIBRARY_PATH '
              'in your .bashrc or .bashprofile file.', file=sys.stderr)
        if sys.platform == 'darwin':
            print('(On Mac OS, the appropriate environment variable may be DYLD_LIBRARY_PATH or DYLD_FALLBACK_LIBRARY_PATH)',
                  file=sys.stderr)

from _stab_tools import *
del os, sys, PATH_TO_LIBRARY, ENV

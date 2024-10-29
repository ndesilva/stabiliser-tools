"""Fast checks and conversions for stabiliser states, Cliffords and Paulis.

Available classes (type 'help(stab_tools.<class>)' to see accessible functions and variables):
    Clifford
    Pauli
    Check_Matrix
    Stabiliser_State
"""

import os
import sys

PATH_TO_LIBRARY = os.path.abspath(os.path.dirname(__file__))
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
        if PATH_TO_LIBRARY not in ENV:
            raise ImportError
    except (KeyError, ImportError):
        print('\nWARNING: You may be about to receive an ImportError. (If you don\'t, you can safely ignore this message.) '
              'To resolve this, try running the following command in your terminal:\n'
              f'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:{PATH_TO_LIBRARY}\n', file=sys.stderr)
        if sys.platform == 'darwin':
            print('(On Mac OS, the appropriate environment variable may be DYLD_LIBRARY_PATH or DYLD_FALLBACK_LIBRARY_PATH)',
                  file=sys.stderr)

from _stab_tools import *
del os, sys, PATH_TO_LIBRARY

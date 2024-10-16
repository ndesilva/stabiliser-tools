# TODO Add documentation here?

import os
import sys

PATH_TO_LIBRARY = os.path.abspath(os.path.dirname(__file__))
# print(f'Adding {PATH_TO_LIBRARY} to path if necessary')
if PATH_TO_LIBRARY not in sys.path:
    sys.path.extend([PATH_TO_LIBRARY])

if sys.platform == 'linux':
    try:
        if PATH_TO_LIBRARY not in os.environ['LD_LIBRARY_PATH']:
            raise ImportError
    except (KeyError, ImportError):
        print('\nWARNING: You may be about to receive an ImportError. (If you don\'t, you can safely ignore this message.) '
            'To resolve this, try running the following command in your terminal:\n'
            f'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:{PATH_TO_LIBRARY}\n', file=sys.stderr)

from fast import *
del os, sys, PATH_TO_LIBRARY

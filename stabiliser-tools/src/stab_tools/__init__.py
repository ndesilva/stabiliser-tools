import os
import sys
PATH_TO_LIBRARY = f'{os.path.abspath(os.path.dirname(__file__))}'
if PATH_TO_LIBRARY not in sys.path:
    # print(f'Adding {PATH_TO_LIBRARY} to path')
    sys.path.extend([PATH_TO_LIBRARY])
from fast import *
del os, sys, PATH_TO_LIBRARY

"""
SymPy is a Python library for symbolic mathematics. It aims to become a
full-featured computer algebra system (CAS) while keeping the code as
simple as possible in order to be comprehensible and easily extensible.
SymPy is written entirely in Python and does not require any external
libraries, except optionally for plotting support.

See the webpage for more information and documentation:

    http://code.google.com/p/sympy/"""

__version__ = "0.7.2"

# Try to determine if 2to3 has been run. To do this, we look at long.__name__.
# If 2to3 has been run, it should convert long to int.
import sys
if sys.version_info[0] == 3:
    try:
        HAS_2TO3 = long.__name__ == "int"
    except NameError: # it tries to see long but long doesn't exist in Python 3
        HAS_2TO3 = False
else:
    HAS_2TO3 = long.__name__ == "int"

if sys.version_info[0] == 2:
    if HAS_2TO3:
        raise ImportError("It appears 2to3 has been run on the codebase. Use "
                          "Python 3 or get the original source code.")
    else:
        if sys.version_info[1] < 5:
            raise ImportError("Python Version 2.5 or above is required for SymPy.")
else: # Python 3
    if not HAS_2TO3:
        raise ImportError("This is the Python 2 version of SymPy. To use SymPy "
    "with Python 3, please obtain a Python 3 version from http://sympy.org, "
    "or use the bin/use2to3 script if you are using the git version.")
    # Here we can also check for specific Python 3 versions, if needed

del sys
del HAS_2TO3


def __sympy_debug():
    # helper function so we don't import os globally
    import os
    return eval(os.getenv('SYMPY_DEBUG', 'False'))
SYMPY_DEBUG = __sympy_debug()

from sympy.core import *
from logic import *
from assumptions import *
from polys import *
from series import *
from functions import *
from ntheory import *
from concrete import *
from simplify import *
from sets import *
from solvers import *
from matrices import *
from geometry import *
from utilities import *
from integrals import *
from tensor import *
from parsing import *
# Adds about .04-.05 seconds of import time
# from combinatorics import *
# This module is slow to import:
#from physics import units
from plotting import plot, Plot, textplot, plot_backends, plot_implicit
from printing import pretty, pretty_print, pprint, pprint_use_unicode, \
    pprint_try_use_unicode, print_gtk, print_tree, pager_print, TableForm
from printing import ccode, fcode, jscode, latex, preview
from printing import python, print_python, srepr, sstr, sstrrepr
from interactive import init_session, init_printing

evalf._create_evalf_table()

# This is slow to import:
#import abc

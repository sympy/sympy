"""SymPy is a Python library for symbolic mathematics. It aims to become a
full-featured computer algebra system (CAS) while keeping the code as
simple as possible in order to be comprehensible and easily extensible.
SymPy is written entirely in Python and does not require any external
libraries, except optionally for plotting support.

See the webpage for more information and documentation:

    http://sympy.org
"""

from __future__ import absolute_import, print_function
del absolute_import, print_function

try:
    import mpmath
except ImportError:
    raise ImportError("SymPy now depends on mpmath as an external library. "
    "See http://docs.sympy.org/latest/install.html#mpmath for more information.")

from sympy.release import __version__

import sys
if sys.version_info[0] == 2 and sys.version_info[1] < 6:
    raise ImportError("Python Version 2.6 or above is required for SymPy.")
else:  # Python 3
    pass
    # Here we can also check for specific Python 3 versions, if needed

sys.modules['sympy.mpmath'] = mpmath

del sys


def __sympy_debug():
    # helper function so we don't import os globally
    import os
    debug_str = os.getenv('SYMPY_DEBUG', 'False')
    if debug_str in ('True', 'False'):
        return eval(debug_str)
    else:
        raise RuntimeError("unrecognized value for SYMPY_DEBUG: %s" %
                           debug_str)
SYMPY_DEBUG = __sympy_debug()

from .core import *
from .logic import *
from .assumptions import *
from .polys import *
from .series import *
from .functions import *
from .ntheory import *
from .concrete import *
from .simplify import *
from .sets import *
from .solvers import *
from .matrices import *
from .geometry import *
from .utilities import *
from .integrals import *
from .tensor import *
from .parsing import *
from .calculus import *
# Adds about .04-.05 seconds of import time
# from combinatorics import *
# This module is slow to import:
#from physics import units
from .plotting import plot, textplot, plot_backends, plot_implicit
from .printing import pretty, pretty_print, pprint, pprint_use_unicode, \
    pprint_try_use_unicode, print_gtk, print_tree, pager_print, TableForm
from .printing import ccode, fcode, jscode, mathematica_code, octave_code, \
    latex, preview
from .printing import python, print_python, srepr, sstr, sstrrepr
from .interactive import init_session, init_printing

evalf._create_evalf_table()

# This is slow to import:
#import abc

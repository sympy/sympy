"""
SymPy is a Python library for symbolic mathematics. It aims to become a
full-featured computer algebra system (CAS) while keeping the code as simple as
possible in order to be comprehensible and easily extensible. SymPy is written
entirely in Python and does not require any external libraries, except
optionally for plotting support.

See the webpage for more information and documentation:

    http://code.google.com/p/sympy/
"""

__version__ = "0.5.4-svn"

from sympy.core import *

from series import *
from concrete import *
from functions import *
from simplify import *
from solvers import *
from matrices import *
from geometry import *
from polynomials import *
from utilities import *
from integrals import *
from plotting import Plot, textplot
from printing import pretty, pretty_print, pprint, pprint_use_unicode

import abc

#for _n, _cls in Basic.singleton.items():
#    exec '%s = _cls()' % (_n)

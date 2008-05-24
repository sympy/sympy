"""
SymPy is a Python library for symbolic mathematics. It aims to become a
full-featured computer algebra system (CAS) while keeping the code as simple as
possible in order to be comprehensible and easily extensible. SymPy is written
entirely in Python and does not require any external libraries, except
optionally for plotting support.

See the webpage for more information and documentation:

    http://code.google.com/p/sympy/
"""

__version__ = "0.5.15"

#put path to pyglet into the search path, so that it can be imported
import os.path
import sys
sys.path.insert(0, os.path.join(os.path.abspath(os.path.dirname(__file__)), "thirdparty", \
        "pyglet"))

import symbol as stdlib_symbol
from sympy.core import *

from series import *
from functions import *
from ntheory import *
from concrete import *
from simplify import *
from solvers import *
from matrices import *
from geometry import *
from polynomials import *
from utilities import *
from integrals import *
# This module is slow to import:
#from physics import units
from plotting import Plot, textplot
from printing import pretty, pretty_print, pprint, pprint_use_unicode, \
    pprint_try_use_unicode, print_gtk
from printing import latex, preview, view, pngview, pdfview, dviview

import abc

from polys import *

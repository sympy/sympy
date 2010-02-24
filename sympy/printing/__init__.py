"""Printing subsystem"""

from pretty import *
from latex import latex, print_latex
from mathml import mathml, print_mathml
from python import python, print_python
from ccode import ccode, print_ccode
from fcode import fcode, print_fcode
from gtk import *

from preview import preview

from str import StrPrinter, sstr, sstrrepr
_StrPrinter = StrPrinter()


from repr import srepr

# /cyclic/
from sympy.core import basic
from sympy.matrices import matrices
basic.StrPrinter = _StrPrinter
matrices.StrPrinter = _StrPrinter
del basic, matrices

from tree import print_tree

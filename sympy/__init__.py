__version__ = "0.5.2"

from sympy.core import *

from concrete import *
from functions import *
from series import *
from simplify import *
from solvers import *
from matrices import *
from geometry import *
from polynomials import *
from utilities import *
from integrals import *
from plotting import Plot
from printing.pretty import pprint

#for _n, _cls in Basic.singleton.items():
#    exec '%s = _cls()' % (_n)
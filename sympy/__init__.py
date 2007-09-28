__version__ = "0.5.3-svn"

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

#for _n, _cls in Basic.singleton.items():
#    exec '%s = _cls()' % (_n)

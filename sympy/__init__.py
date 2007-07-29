
__version__ = "0.5.0"

from sympy.core import *

from series import *
from simplify import *
from solvers import *
from matrices import *
from geometry import *
from polynomials import *
from utilities import *
#from specfun import *
from integrals import *

try:
    from plotting import Plot
except ImportError, e:
    class Plot(object):
        def __init__(*args, **kwargs):
            raise e

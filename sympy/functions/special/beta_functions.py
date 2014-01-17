from __future__ import print_function, division

from sympy.core import Add, S, C, sympify, oo, pi
from sympy.core.function import Function, ArgumentIndexError
from .zeta_functions import zeta
from .error_functions import erf
from sympy.core import Dummy, Rational
from sympy.functions.elementary.exponential import log
from sympy.functions.elementary.integers import floor
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.combinatorial.numbers import bernoulli
from sympy.functions.combinatorial.factorials import rf
from sympy.functions.combinatorial.numbers import harmonic
from sympy.core.compatibility import xrange

###############################################################################
############################ COMPLETE BETA  FUNCTION ##########################
###############################################################################
class beta(Function):
    """Beta function or Euler's first integral is closely associated with gamma function.

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Beta_function
    """
    nargs = 2
    unbranched = True

    def fdiff(self, argindex=2):
        x, y = self.args
        if argindex == 1:
            # Diff wrt x
            return beta(x, y)*(digamma(x) - digamma(x+y)) 
        elif argindex == 2:
            # Diff wrt y
            return beta(x, y)*(digamma(y) - digamma(x+y)) 
        else:
            raise ArgumentIndexError(self, argindex)

	@classmethod
    def eval(cls, x, y):
        return gamma(x)*gamma(y) / gamma(x + y)
         
        
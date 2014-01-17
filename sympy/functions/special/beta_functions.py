from __future__ import print_function, division

from sympy.core import Add, S, C, sympify, oo, pi
from sympy.core.function import Function, ArgumentIndexError
from sympy.functions.special.gamma_functions import gamma, digamma

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

    def fdiff(self, argindex):
        x, y = self.args 
        if argindex == 1:
            # Diff wrt x
            return beta(x, y)*(digamma(x) - digamma(x + y)) 
        elif argindex == 2:
            # Diff wrt y
            return beta(x, y)*(digamma(y) - digamma(x + y)) 
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, x, y):
        return gamma(x)*gamma(y) / gamma(x + y)

    def _eval_expand_func(self, **hints):
        x, y = self.args
        return gamma(x)*gamma(y) / gamma(x + y) 
         
        
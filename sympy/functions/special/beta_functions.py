from __future__ import print_function, division

from sympy.core.function import Function, ArgumentIndexError
from sympy.functions.special.gamma_functions import gamma, digamma

###############################################################################
############################ COMPLETE BETA  FUNCTION ##########################
###############################################################################

class beta(Function):
    r"""
    The beta function
    
    .. math::
        \Beta(x,y) := \int^{1}_{0} t^{x-1} (1-t)^{y-1} \mathrm{d}t.
        
                
    Beta function or Euler's first integral is closely associated with gamma function.
    The Beta function often used in probability theory and mathematical statistics. It satisfies
    properties like:
    B(a,1) = 1/a;
    B(a,b) = B(b,a) etc.

    Examples
    ========

    >>> from sympy import I, pi
    >>> from sympy.abc import x,y

    The Beta function obeys the mirror symmetry:
    >>> from sympy import conjugate
    >>> conjugate(beta(x,y))
    beta(conjugate(x), conjugate(y))

    Differentiation with respect to both x and y is supported:

    >>> from sympy import diff
    >>> diff(beta(x,y), x)
    (polygamma(0, x) - polygamma(0, x + y))*beta(x, y)

    >>> from sympy import diff
    >>> diff(beta(x,y), y)
    (polygamma(0, y) - polygamma(0, x + y))*beta(x, y)

    We can numerically evaluate the gamma function to arbitrary precision
    on the whole complex plane:

    >>> beta(pi,pi).evalf(40)
    0.02671848900111377452242355235388489324562

    >>> beta(1+I,1+I).evalf(20)
    -0.2112723729365330143 - 0.7655283165378005676*I


    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Beta_function
    .. [2] http://mathworld.wolfram.com/BetaFunction.html
    .. [3] http://dlmf.nist.gov/5.12
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
        pass

    def _eval_expand_func(self, **hints):
        x, y = self.args
        return gamma(x)*gamma(y) / gamma(x + y)

    def _eval_is_real(self):
        return self.args[0].is_real and self.args[1].is_real

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate(), self.args[1].conjugate())

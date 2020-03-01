from __future__ import print_function, division
from sympy.core.function import Function, ArgumentIndexError
from sympy.functions.special.gamma_functions import gamma, digamma
from sympy.functions.special.hyper import hyper
from sympy.functions.elementary.exponential import log
from sympy import Symbol, S

###############################################################################
############################ COMPLETE BETA  FUNCTION ##########################
###############################################################################

class beta(Function):
    r"""
    The beta integral is called the Eulerian integral of the first kind by
    Legendre:

    .. math::
        \mathrm{B}(x,y) := \int^{1}_{0} t^{x-1} (1-t)^{y-1} \mathrm{d}t.

    Explanation
    ===========

    The Beta function or Euler's first integral is closely associated
    with the gamma function. The Beta function is often used in probability
    theory and mathematical statistics. It satisfies properties like:

    .. math::
        \mathrm{B}(a,1) = \frac{1}{a} \\
        \mathrm{B}(a,b) = \mathrm{B}(b,a)  \\
        \mathrm{B}(a,b) = \frac{\Gamma(a) \Gamma(b)}{\Gamma(a+b)}

    Therefore for integral values of $a$ and $b$:

    .. math::
        \mathrm{B} = \frac{(a-1)! (b-1)!}{(a+b-1)!}

    Examples
    ========

    >>> from sympy import I, pi
    >>> from sympy.abc import x, y

    The Beta function obeys the mirror symmetry:

    >>> from sympy import beta
    >>> from sympy import conjugate
    >>> conjugate(beta(x, y))
    beta(conjugate(x), conjugate(y))

    Differentiation with respect to both $x$ and $y$ is supported:

    >>> from sympy import beta
    >>> from sympy import diff
    >>> diff(beta(x, y), x)
    (polygamma(0, x) - polygamma(0, x + y))*beta(x, y)

    >>> from sympy import beta
    >>> from sympy import diff
    >>> diff(beta(x, y), y)
    (polygamma(0, y) - polygamma(0, x + y))*beta(x, y)

    We can numerically evaluate the gamma function to arbitrary precision
    on the whole complex plane:

    >>> from sympy import beta
    >>> beta(pi, pi).evalf(40)
    0.02671848900111377452242355235388489324562

    >>> beta(1 + I, 1 + I).evalf(20)
    -0.2112723729365330143 - 0.7655283165378005676*I

    See Also
    ========

    gamma: Gamma function.
    uppergamma: Upper incomplete gamma function.
    lowergamma: Lower incomplete gamma function.
    polygamma: Polygamma function.
    loggamma: Log Gamma function.
    digamma: Digamma function.
    trigamma: Trigamma function.

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Beta_function
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
        if y is S.One:
            return 1/x
        if x is S.One:
            return 1/y

    def _eval_expand_func(self, **hints):
        x, y = self.args
        return gamma(x)*gamma(y) / gamma(x + y)

    def _eval_is_real(self):
        return self.args[0].is_real and self.args[1].is_real

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate(), self.args[1].conjugate())

    def _eval_rewrite_as_gamma(self, x, y, **kwargs):
        return self._eval_expand_func(**kwargs)


###############################################################################
########################### INCOMPLETE BETA  FUNCTION #########################
###############################################################################



class betainc(Function):
    r"""
    The incomplete beta function is a generalization of the beta function,

    .. math::
        \mathrm{B}_{z}(a,b) := \int^{z}_{0} t^{a-1} (1-t)^{b-1} \mathrm{d}t.

    Incomplete beta function satisfies properties like:

    .. math::
        \mathrm{B}_{1}(a,b) = \mathrm{B}(a,b) \\

    See Also
    ========
    Gamma function, Upper incomplete gamma function, Lower incomplete gamma function,
    Polygamma function, Log Gamma function, Digamma function, Trigamma function.

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Beta_function#Incomplete_beta_function
    .. [2] http://mathworld.wolfram.com/IncompleteBetaFunction.html
    """
    nargs = 3

    def fdiff(self, argindex):
        z, a, b = self.args
        if argindex not in (1,2,3):
            raise ArgumentIndexError(self,argindex)
        elif argindex == 1:
            # Diff wrt z
            return (1 - z)**(b - 1)*z**(a - 1)
        elif argindex == 2:
            # Diff wrt a
            term1 = betainc(z, a, b)*log(z)
            term2 = z**a*gamma(a)**2*hyper([a, a, 1-b], [a+1, a+1], z)
            return  term1 - term2
        elif argindex == 3:
            term1 = gamma(b)**2*(1 - z)**b*hyper([b, b, 1-a], [b+1, b+1], 1-z)
            term2 = log(1 - z)*self.func(1 - z, a, b)
            term3 = (digamma(b) - digamma(a + b))*beta(a, b)
            return term1 - term2 + term3

    @classmethod
    def eval(cls, z, a, b):
        #Case z = 0 , incomplete beta is zero
        if z == 0:
            return S.Zero
        #Beta function is special case of Incomplete beta function for z=1
        if z == 1:
            return beta(a,b)

    def _eval_rewrite_as_Integral(self,*args, **kwargs):
        z, a, b = self.args
        from sympy import Integral
        t = Symbol('t')
        integrand = t**(a - 1)*(1 - t)**(b - 1)
        return Integral(integrand,(t, 0, z))

    def _eval_rewrite_as_hyper(self,*args, **kwargs):
        z, a, b = self.args
        return gamma(a)*z**a*hyper([a, 1-b], [a+1], z)

    def _eval_is_real(self):
        #return true if a and b > 0 and 1 >= z >= 0 else None
        z, a, b = self.args
        return (a.is_positive and b.is_positive and
                z.is_nonnegative and (1-z).is_nonnegative)


    def _eval_conjugate(self):
        z, a, b = self.args
        return self.func(z.conjugate(), a.conjugate(), b.conjugate())

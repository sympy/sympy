from sympy.core import S
from sympy.core.function import Function, ArgumentIndexError
from sympy.core.symbol import Dummy
from sympy.functions.special.gamma_functions import gamma, digamma

###############################################################################
############################ COMPLETE BETA  FUNCTION ##########################
###############################################################################

class beta(Function):
    r"""
    The beta integral is called the Eulerian integral of the first kind by
    Legendre:

    .. math::
        \mathrm{B}(x,y)  \int^{1}_{0} t^{x-1} (1-t)^{y-1} \mathrm{d}t.

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

    A special case of the Beta function when `x = y` is the
    Central Beta function. It satisfies properties like:

    .. math::
        \mathrm{B}(x) = 2^{1 - 2x}\mathrm{B}(x, \frac{1}{2})
        \mathrm{B}(x) = 2^{1 - 2x} cos(\pi x) \mathrm{B}(\frac{1}{2} - x, x)
        \mathrm{B}(x) = \int_{0}^{1} \frac{t^x}{(1 + t)^{2x}} dt
        \mathrm{B}(x) = \frac{2}{x} \prod_{n = 1}^{\infty} \frac{n(n + 2x)}{(n + x)^2}

    Examples
    ========

    >>> from sympy import I, pi
    >>> from sympy.abc import x, y

    The Beta function obeys the mirror symmetry:

    >>> from sympy import beta, conjugate
    >>> conjugate(beta(x, y))
    beta(conjugate(x), conjugate(y))

    Differentiation with respect to both $x$ and $y$ is supported:

    >>> from sympy import beta, diff
    >>> diff(beta(x, y), x)
    (polygamma(0, x) - polygamma(0, x + y))*beta(x, y)

    >>> diff(beta(x, y), y)
    (polygamma(0, y) - polygamma(0, x + y))*beta(x, y)

    >>> diff(beta(x), x)
    2*(polygamma(0, x) - polygamma(0, 2*x))*beta(x, x)

    We can numerically evaluate the gamma function to arbitrary precision
    on the whole complex plane:

    >>> from sympy import beta
    >>> beta(pi).evalf(40)
    0.02671848900111377452242355235388489324562

    >>> beta(1 + I).evalf(20)
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

    def __new__(cls, x, y=None):
        if y == None:
            y = x
        return super(cls, cls).__new__(cls, x, y)

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

    def _eval_rewrite_as_gamma(self, x, y, piecewise=True, **kwargs):
        return self._eval_expand_func(**kwargs)

    def _eval_rewrite_as_Integral(self, x, y, **kwargs):
        from sympy.integrals.integrals import Integral
        t = Dummy('t')
        return Integral(t**(x - 1)*(1 - t)**(y - 1), (t, 0, 1))

###############################################################################
########################## INCOMPLETE BETA FUNCTION ###########################
###############################################################################

class betainc(Function):
    r"""
    The Generalized Incomplete Beta function is defined as

    .. math::
        \mathrm{B}_{(x_1, x_2)}(a, b) = \int_{x_1}^{x_2} t^{a - 1} (1 - t)^{b - 1} dt

    The Beta function is a special case of the
    Generalized Incomplete Beta function :

    .. math:: \mathrm{B}(a, b) = \mathrm{B}_{(0, 1)}(a, b)

    The Generalized Regularized Incomplete Beta function is given by

    .. math::
        \mathrm{I}_{(x_1, x_2)}(a, b) = \frac{\mathrm{B}_{(x_1, x_2)}(a, b)}{\mathrm{B}(a, b)}

    The Regularized Incomplete Beta function is the cumulative distribution
    function of the beta distribution.

    Examples
    ========

    >>> from sympy import I, pi, E, symbols
    >>> from sympy import betainc, conjugate
    >>> a, b, x, x1, x2 = symbols('a b x x1 x2')

    The Generalized Incomplete Beta function is
    given by

    >>> betainc(a, b, x1, x2)
    betainc(a, b, x1, x2, 0)

    The Generalized Regularized Incomplete Beta function
    is given by

    >>> betainc(a, b, x1, x2, regularized=True)
    betainc(a, b, x1, x2, 1)

    The Incomplete Beta function is a special case
    of the Generalized Incomplete Beta Function.
    It can be obtained as follows:

    >>> betainc(a, b, 0, x)
    betainc(a, b, 0, x, 0)

    The Regularized Incomplete Beta function is a special case
    of the Generalized Regularized Incomplete Beta Function.
    It can be obtained as follows:

    >>> betainc(a, b, 0, x, regularized=True)
    betainc(a, b, 0, x, 1)

    The Beta function obeys the mirror symmetry:

    >>> conjugate(betainc(a, b, x1, x2))
    betainc(conjugate(a), conjugate(b), conjugate(x1), conjugate(x2), 0)

    We can numerically evaluate the Incomplete Beta function
    to arbitrary precision on the whole complex plane:

    >>> from sympy import betainc
    >>> betainc(2, 3, 4, 5).evalf(10)
    56.08333333
    >>> betainc(0.75, 1 - 4*I, 0, 2 + 3*I).evalf(25)
    0.2241657956955709603655887 + 0.3619619242700451992411724*I
    >>> betainc(pi, E, 0, 1, regularized=True).evalf(5)
    1.00000

    See Also
    ========

    beta: Beta function
    hyper: Generalized Hypergeometric function

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Beta_function#Incomplete_beta_function
    .. [2] https://dlmf.nist.gov/8.17
    .. [3] https://functions.wolfram.com/GammaBetaErf/Beta4/
    .. [4] https://functions.wolfram.com/GammaBetaErf/BetaRegularized4/02/

    """
    nargs = 5
    unbranched = True

    def __new__(cls, a, b, x1, x2, regularized=False):
        return Function.__new__(cls, a, b, x1, x2, int(regularized == True))

    def _eval_is_real(self):
        return all(arg.is_real for arg in self.args)

    def _eval_conjugate(self):
        a, b, x1, x2, reg = [arg.conjugate() for arg in self.args]
        return self.func(a, b, x1, x2, reg == 1)

    def _eval_rewrite_as_Integral(self, a, b, x1, x2, reg, **kwargs):
        from sympy.integrals.integrals import Integral
        t = Dummy('t')
        expr = Integral(t**(a - 1)*(1 - t)**(b - 1), (t, x1, x2))
        if reg == 1: #Regularized
            return expr / Integral(t**(a - 1)*(1 - t)**(b - 1), (t, 0, 1))
        else:
            return expr

    def _eval_rewrite_as_hyper(self, a, b, x1, x2, reg, **kwargs):
        from sympy.functions.special.hyper import hyper
        expr = (x2**a * hyper((a, 1 - b), (a + 1,), x2) - x1**a * hyper((a, 1 - b), (a + 1,), x1)) / a
        if reg == 1: #Regularized
            return expr / beta(a, b)
        else:
            return expr

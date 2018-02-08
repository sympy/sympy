from __future__ import print_function, division

from sympy.stats.drv import SingleDiscreteDistribution, SingleDiscretePSpace
from sympy import factorial, exp, S, sympify, floor, sqrt, log, uppergamma, sign
from sympy.stats.rv import _value_check
from sympy.sets.sets import Interval
from mpmath import pi
import random

__all__ = ['Geometric', 'Poisson']


def rv(symbol, cls, *args):
    args = list(map(sympify, args))
    dist = cls(*args)
    dist.check(*args)
    return SingleDiscretePSpace(symbol, dist).value


class PoissonDistribution(SingleDiscreteDistribution):
    _argnames = ('lamda',)

    set = S.Naturals0

    @staticmethod
    def check(lamda):
        _value_check(lamda > 0, "Lambda must be positive")

    def pdf(self, k):
        return self.lamda**k / factorial(k) * exp(-self.lamda)

    def _cdf(self, k):
        return uppergamma(floor(k+1), self.lamda) / factorial(floor(k))

    def sample(self):
        if self.lamda < 30:
            """
            Poisson Random Variate Generation
            C. D. Kemp and Adrienne W. Kemp
            Journal of the Royal Statistical Society. Series C (Applied Statistics)
            Vol. 40, No. 1 (1991), pp. 143-158
            Simple sequential search (Algorithm 3 in the above paper) - O(lamda) so only for use in low lamda cases
            """
            p_0 = exp(-self.lamda)
            p = p_0
            F = p
            x = 0
            u = random.uniform(0,1)
            while True:
                if u > F:
                    x = x+1
                    p = self.lamda*p/x
                    F = F + p
                else:
                    return x
        else:
            """
            https://pdfs.semanticscholar.org/00f1/8468dd53e4e3f0c28a5d504f8cd10376a673.pdf
            Algorithm implements PTRD described in the above paper
            """
            log_factorial_dict = {}
            for k in range(10):
                log_factorial_dict[k]=log(factorial(k))
            smu = sqrt(self.lamda)
            b = 0.931 + 2.53*smu
            a = -0.059 + 0.02483*b
            alpha = 1.1239 + 1.1328 / (b - 3.4) #alpha here is defined as 1/alpha in the paper
            v_r = 0.9277 - 3.6224 / (b - 2)
            while True:
                V = random.uniform(0, 1)
                if V <= 0.86 * v_r:
                    U = V / v_r - 0.43
                    return floor((2 * a/ (0.5 - abs(U)) + b) * U + self.lamda + 0.445)
                elif V >= v_r:
                    U = random.uniform(-0.5, 0.5)
                else:
                    U = V / v_r - 0.93
                    U = sign(U) * 0.5 - U
                    V = random.uniform(0, v_r)

                us = 0.5 - abs(U)
                if us < 0.013 and V > us:
                    break

                k = floor((2 * a / us + b) * U + self.lamda + 0.445)
                V = V * alpha / (a / us**2 + b)

                if k >= 10 and log(V * smu) <= (k + 0.5) * log(self.lamda / k) - self.lamda - log(sqrt(2 * pi)) + k -(1/12 - 1/(360 * k**2)) / k:
                    return k
                elif k>= 0 and k <= 9 and log(V) <= k * log(self.lamda) - self.lamda - log_factorial_dict[k]:
                    return k
                else:
                    continue

def Poisson(name, lamda):
    r"""
    Create a discrete random variable with a Poisson distribution.

    The density of the Poisson distribution is given by

    .. math::
        f(k) := \frac{\lambda^{k} e^{- \lambda}}{k!}

    Parameters
    ==========

    lamda: Positive number, a rate

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Poisson, density, E, variance
    >>> from sympy import Symbol, simplify

    >>> rate = Symbol("lambda", positive=True)
    >>> z = Symbol("z")

    >>> X = Poisson("x", rate)

    >>> density(X)(z)
    lambda**z*exp(-lambda)/factorial(z)

    >>> E(X)
    lambda

    >>> simplify(variance(X))
    lambda

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Poisson_distribution
    [2] http://mathworld.wolfram.com/PoissonDistribution.html
    """
    return rv(name, PoissonDistribution, lamda)


class GeometricDistribution(SingleDiscreteDistribution):
    _argnames = ('p',)
    set = S.Naturals

    @staticmethod
    def check(p):
        _value_check(0 < p and p <= 1, "p must be between 0 and 1")

    def pdf(self, k):
        return (1 - self.p)**(k - 1) * self.p


def Geometric(name, p):
    r"""
    Create a discrete random variable with a Geometric distribution.

    The density of the Geometric distribution is given by

    .. math::
        f(k) := p (1 - p)^{k - 1}

    Parameters
    ==========

    p: A probability between 0 and 1

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Geometric, density, E, variance
    >>> from sympy import Symbol, S

    >>> p = S.One / 5
    >>> z = Symbol("z")

    >>> X = Geometric("x", p)

    >>> density(X)(z)
    (4/5)**(z - 1)/5

    >>> E(X)
    5

    >>> variance(X)
    20

    References
    ==========

    [1] http://en.wikipedia.org/wiki/Geometric_distribution
    [2] http://mathworld.wolfram.com/GeometricDistribution.html
    """
    return rv(name, GeometricDistribution, p)

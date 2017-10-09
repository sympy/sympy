from __future__ import print_function, division

from sympy.stats.drv import SingleDiscreteDistribution, SingleDiscretePSpace
from sympy import factorial, exp, S, sympify
from sympy.stats.rv import _value_check

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

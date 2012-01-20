from rv import (P, E, Density, Where, Given, pspace, CDF, Sample, sample_iter, random_symbols, independent, dependent)
from sympy import sqrt

__all__ = ['P', 'E', 'Density', 'Where', 'Given', 'Sample', 'CDF', 'pspace',
        'sample_iter', 'Var', 'Std', 'Skewness', 'Covar', 'dependent',
        'independent', 'random_symbols']

def variance(X, given=None, **kwargs):
    """
    Variance of a random expression

    Expectation of (X-E(X))**2

    Examples
    ========

    >>> from sympy.stats import Die, E, Bernoulli, Var
    >>> from sympy import simplify, Symbol

    >>> X = Die(6)
    >>> p = Symbol('p')
    >>> B = Bernoulli(p, 1, 0)

    >>> Var(2*X)
    35/3

    >>> simplify(Var(B))
    p*(-p + 1)
    """
    return E(X**2, given, **kwargs) - E(X, given, **kwargs)**2
Var = variance


def standard_deviation(X, given=None, **kwargs):
    """
    Standard Deviation of a random expression

    Square root of the Expectation of (X-E(X))**2

    Examples
    ========

    >>> from sympy.stats import Bernoulli, Std
    >>> from sympy import Symbol

    >>> p = Symbol('p')
    >>> B = Bernoulli(p, 1, 0)

    >>> Std(B)
    sqrt(-p**2 + p)
    """
    return sqrt(variance(X, given, **kwargs))
Std = standard_deviation

def covariance(X, Y, given=None, **kwargs):
    """
    Covariance of two random expressions

    The expectation that the two variables will rise and fall together

    Covariance(X,Y) = E( (X-E(X)) * (Y-E(Y)) )

    Examples
    ========

    >>> from sympy.stats import Exponential, Covar
    >>> from sympy import Symbol

    >>> rate = Symbol('lambda', positive=True, real=True, bounded = True)
    >>> X = Exponential(rate)
    >>> Y = Exponential(rate)

    >>> Covar(X, X)
    lambda**(-2)
    >>> Covar(X, Y)
    0
    >>> Covar(X, Y + rate*X)
    1/lambda
    """

    return E( (X-E(X, given, **kwargs)) * (Y-E(Y, given, **kwargs)),
            given, **kwargs)
Covar = covariance

def skewness(X, given=None, **kwargs):

    mu = E(X, given, **kwargs)
    sigma = Std(X, given, **kwargs)
    return E( ((X-mu)/sigma) ** 3 , given, **kwargs)
Skewness = skewness


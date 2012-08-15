from rv import (probability, expectation, density, where, given, pspace, cdf,
        sample, sample_iter, random_symbols, independent, dependent)
from sympy import sqrt

__all__ = ['P', 'E', 'density', 'where', 'given', 'sample', 'cdf', 'pspace',
        'sample_iter', 'variance', 'std', 'skewness', 'covariance', 'dependent',
        'independent', 'random_symbols']

def variance(X, condition=None, **kwargs):
    """
    Variance of a random expression

    Expectation of (X-E(X))**2

    Examples
    ========

    >>> from sympy.stats import Die, E, Bernoulli, variance
    >>> from sympy import simplify, Symbol

    >>> X = Die('X', 6)
    >>> p = Symbol('p')
    >>> B = Bernoulli('B', p, 1, 0)

    >>> variance(2*X)
    35/3

    >>> simplify(variance(B))
    p*(-p + 1)
    """
    return (expectation(X**2, condition, **kwargs) -
            expectation(X, condition, **kwargs)**2)

def standard_deviation(X, condition=None, **kwargs):
    """
    Standard Deviation of a random expression

    Square root of the Expectation of (X-E(X))**2

    Examples
    ========

    >>> from sympy.stats import Bernoulli, std
    >>> from sympy import Symbol

    >>> p = Symbol('p')
    >>> B = Bernoulli('B', p, 1, 0)

    >>> std(B)
    sqrt(-p**2 + p)
    """
    return sqrt(variance(X, condition, **kwargs))
std = standard_deviation

def covariance(X, Y, condition=None, **kwargs):
    """
    Covariance of two random expressions

    The expectation that the two variables will rise and fall together

    Covariance(X,Y) = E( (X-E(X)) * (Y-E(Y)) )

    Examples
    ========

    >>> from sympy.stats import Exponential, covariance
    >>> from sympy import Symbol

    >>> rate = Symbol('lambda', positive=True, real=True, bounded = True)
    >>> X = Exponential('X', rate)
    >>> Y = Exponential('Y', rate)

    >>> covariance(X, X)
    lambda**(-2)
    >>> covariance(X, Y)
    0
    >>> covariance(X, Y + rate*X)
    1/lambda
    """

    return expectation(
            (X - expectation(X, condition, **kwargs)) *
            (Y - expectation(Y, condition, **kwargs)),
                      condition, **kwargs)

def skewness(X, condition=None, **kwargs):

    mu = expectation(X, condition, **kwargs)
    sigma = std(X, condition, **kwargs)
    return expectation( ((X-mu)/sigma) ** 3 , condition, **kwargs)

P = probability
E = expectation

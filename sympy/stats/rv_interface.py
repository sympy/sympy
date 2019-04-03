from __future__ import print_function, division

from .rv import (probability, expectation, density, where, given, pspace, cdf,
        characteristic_function, sample, sample_iter, random_symbols, independent, dependent,
        sampling_density, moment_generating_function, _value_check)
from sympy import Piecewise, sqrt, solveset, Symbol, S

__all__ = ['P', 'E', 'density', 'where', 'given', 'sample', 'cdf', 'characteristic_function', 'pspace',
        'sample_iter', 'variance', 'std', 'skewness', 'covariance',
        'dependent', 'independent', 'random_symbols', 'correlation',
        'moment', 'cmoment', 'sampling_density', 'moment_generating_function',
        'quantile']



def moment(X, n, c=0, condition=None, **kwargs):
    """
    Return the nth moment of a random expression about c i.e. E((X-c)**n)
    Default value of c is 0.

    Examples
    ========

    >>> from sympy.stats import Die, moment, E
    >>> X = Die('X', 6)
    >>> moment(X, 1, 6)
    -5/2
    >>> moment(X, 2)
    91/6
    >>> moment(X, 1) == E(X)
    True
    """
    return expectation((X - c)**n, condition, **kwargs)


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
    p*(1 - p)
    """
    return cmoment(X, 2, condition, **kwargs)


def standard_deviation(X, condition=None, **kwargs):
    """
    Standard Deviation of a random expression

    Square root of the Expectation of (X-E(X))**2

    Examples
    ========

    >>> from sympy.stats import Bernoulli, std
    >>> from sympy import Symbol, simplify

    >>> p = Symbol('p')
    >>> B = Bernoulli('B', p, 1, 0)

    >>> simplify(std(B))
    sqrt(p*(1 - p))
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

    >>> rate = Symbol('lambda', positive=True, real=True, finite=True)
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


def correlation(X, Y, condition=None, **kwargs):
    """
    Correlation of two random expressions, also known as correlation
    coefficient or Pearson's correlation

    The normalized expectation that the two variables will rise
    and fall together

    Correlation(X,Y) = E( (X-E(X)) * (Y-E(Y)) / (sigma(X) * sigma(Y)) )

    Examples
    ========

    >>> from sympy.stats import Exponential, correlation
    >>> from sympy import Symbol

    >>> rate = Symbol('lambda', positive=True, real=True, finite=True)
    >>> X = Exponential('X', rate)
    >>> Y = Exponential('Y', rate)

    >>> correlation(X, X)
    1
    >>> correlation(X, Y)
    0
    >>> correlation(X, Y + rate*X)
    1/sqrt(1 + lambda**(-2))
    """
    return covariance(X, Y, condition, **kwargs)/(std(X, condition, **kwargs)
     * std(Y, condition, **kwargs))


def cmoment(X, n, condition=None, **kwargs):
    """
    Return the nth central moment of a random expression about its mean
    i.e. E((X - E(X))**n)

    Examples
    ========

    >>> from sympy.stats import Die, cmoment, variance
    >>> X = Die('X', 6)
    >>> cmoment(X, 3)
    0
    >>> cmoment(X, 2)
    35/12
    >>> cmoment(X, 2) == variance(X)
    True
    """
    mu = expectation(X, condition, **kwargs)
    return moment(X, n, mu, condition, **kwargs)


def smoment(X, n, condition=None, **kwargs):
    """
    Return the nth Standardized moment of a random expression i.e.
    E( ((X - mu)/sigma(X))**n )

    Examples
    ========

    >>> from sympy.stats import skewness, Exponential, smoment
    >>> from sympy import Symbol
    >>> rate = Symbol('lambda', positive=True, real=True, finite=True)
    >>> Y = Exponential('Y', rate)
    >>> smoment(Y, 4)
    9
    >>> smoment(Y, 4) == smoment(3*Y, 4)
    True
    >>> smoment(Y, 3) == skewness(Y)
    True
    """
    sigma = std(X, condition, **kwargs)
    return (1/sigma)**n*cmoment(X, n, condition, **kwargs)

def skewness(X, condition=None, **kwargs):
    """
    Measure of the asymmetry of the probability distribution

    Positive skew indicates that most of the values lie to the right of
    the mean

    skewness(X) = E( ((X - E(X))/sigma)**3 )

    Examples
    ========

    >>> from sympy.stats import skewness, Exponential, Normal
    >>> from sympy import Symbol
    >>> X = Normal('X', 0, 1)
    >>> skewness(X)
    0
    >>> rate = Symbol('lambda', positive=True, real=True, finite=True)
    >>> Y = Exponential('Y', rate)
    >>> skewness(Y)
    2
    """
    return smoment(X, 3, condition, **kwargs)

def quantile(X, p):
    r"""
    Return the :math:`p^{th}` order quantile of a probability distribution.

    Quantile is defined as the value at which the probability of the random
    variable is less than or equal to the given probability.

    ..math::
        Q(p) = inf{x \in (-\infty, \infty) such that p < F(x)}

    Examples
    ========

    >>> from sympy.stats import quantile, Die, Exponential
    >>> from sympy import Symbol, pprint
    >>> p = Symbol("p", positive=True)

    >>> l = Symbol("lambda", positive=True)
    >>> X = Exponential("x", l)
    >>> pprint(quantile(X, p), use_unicode=False)
               -log(1 - p)
    [0, oo) n {------------}
                  lambda

    >>> D = Die("d", 6)
    >>> pprint(quantile(D, p), use_unicode=False)
    /1  for p <= 1/6
    |
    |2  for p <= 1/3
    |
    |3  for p <= 1/2
    <
    |4  for p <= 2/3
    |
    |5  for p <= 5/6
    |
    \6   for p < 1
    """
    _value_check(p > 0, "The order p must be positive.")
    _value_check(p < 1, "The order p must be less than 1.")

    if pspace(X).is_Continuous or pspace(X).is_Discrete:
        x = Symbol("x")
        return solveset(cdf(X)(x) - p, x, S.Reals)
    else:
        set = tuple()
        for key, value in cdf(X).items():
            if value == S.One:
                set = set + ((key, p < value), )
            else:
                set = set + ((key, p <= value), )
        return Piecewise(*set)


P = probability
E = expectation

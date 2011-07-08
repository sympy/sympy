from rv import P, E, Density, Where, Given#, sample
from sympy import sqrt

def variance(X, given=None):
    """Variance of a random expression.

    Expectation of (X-E(X))**2

    >>> from sympy.statistics import Die, E, Bernoulli, Symbol, var, simplify
    >>> X = Die(6).value
    >>> p = Symbol('p')
    >>> B = Bernoulli(p, 1, 0).value

    >>> var(2*X)
    35/3

    >>> simplify(var(B))
    p*(-p + 1)

    """
    return E(X**2, given) - E(X, given)**2
var = variance


def standard_deviation(X, given=None):
    """Standard Deviation of a random expression.

    Square root of the Expectation of (X-E(X))**2

    >>> from sympy.statistics import Bernoulli, std, Symbol
    >>> p = Symbol('p')
    >>> B = Bernoulli(p, 1, 0).value

    >>> std(B)
    (-p**2 + p)**(1/2)

    """
    return sqrt(variance(X, given))
std = standard_deviation

def covariance(X, Y, given=None):
    """Covariance of two random expressions.

    The expectation that the two variables will rise and fall together

    Covariance(X,Y) = E( (X-E(X)) * (Y-E(Y)) )


    >>> from sympy.statistics import ExponentialProbabilitySpace, covar, Symbol
    >>> rate = Symbol('lambda', positive=True, real=True, bounded = True)
    >>> X = ExponentialProbabilitySpace(rate).value
    >>> Y = ExponentialProbabilitySpace(rate).value

    >>> covar(X, X)
    lambda**(-2)
    >>> covar(X, Y)
    0
    >>> covar(X, Y + rate*X)
    1/lambda
    """

    return E( (X-E(X, given)) * (Y-E(Y, given)), given )
covar = covariance

def skewness(X, given=None):

    mu = E(X, given)
    sigma = std(X, given)
    return E( ((X-mu)/sigma) ** 3 , given)


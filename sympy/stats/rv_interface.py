from rv import P, E, Density, Where, Given, pspace, CDF, Sample, sample_iter
from sympy import sqrt

def variance(X, given=None, **kwargs):
    """Variance of a random expression.

    Expectation of (X-E(X))**2

    >>> from sympy.stats import Die, E, Bernoulli, var
    >>> from sympy import simplify, Symbol

    >>> X = Die(6)
    >>> p = Symbol('p')
    >>> B = Bernoulli(p, 1, 0)

    >>> var(2*X)
    35/3

    >>> simplify(var(B))
    p*(-p + 1)

    """
    return E(X**2, given, **kwargs) - E(X, given, **kwargs)**2
var = variance


def standard_deviation(X, given=None, **kwargs):
    """Standard Deviation of a random expression.

    Square root of the Expectation of (X-E(X))**2

    >>> from sympy.stats import Bernoulli, std
    >>> from sympy import Symbol

    >>> p = Symbol('p')
    >>> B = Bernoulli(p, 1, 0)

    >>> std(B)
    (-p**2 + p)**(1/2)

    """
    return sqrt(variance(X, given, **kwargs))
std = standard_deviation

def covariance(X, Y, given=None, **kwargs):
    """Covariance of two random expressions.

    The expectation that the two variables will rise and fall together

    Covariance(X,Y) = E( (X-E(X)) * (Y-E(Y)) )


    >>> from sympy.stats import Exponential, covar
    >>> from sympy import Symbol

    >>> rate = Symbol('lambda', positive=True, real=True, bounded = True)
    >>> X = Exponential(rate)
    >>> Y = Exponential(rate)

    >>> covar(X, X)
    lambda**(-2)
    >>> covar(X, Y)
    0
    >>> covar(X, Y + rate*X)
    1/lambda
    """

    return E( (X-E(X, given, **kwargs)) * (Y-E(Y, given, **kwargs)),
            given, **kwargs)
covar = covariance

def skewness(X, given=None, **kwargs):

    mu = E(X, given, **kwargs)
    sigma = std(X, given, **kwargs)
    return E( ((X-mu)/sigma) ** 3 , given, **kwargs)

def dependent(a, b):
    """Dependence of two random expressions

    Two expressions are independent if knowledge of one does not change
    computations on the other

    >>> from sympy.stats import Die, dependent, Given
    >>> from sympy import Tuple

    >>> X, Y = Die(6), Die(6)
    >>> dependent(X, Y)
    False
    >>> dependent(2*X + Y, -Y)
    True
    >>> X, Y = Given(Tuple(X, Y), X>Y)
    >>> dependent(X, Y)
    True

    See Also:
        independent
    """
    a_symbols = pspace(b).symbols
    b_symbols = pspace(a).symbols
    return len(a_symbols.intersect(b_symbols)) != 0

def independent(a, b):
    """Independence of two random expressions

    Two expressions are independent if knowledge of one does not change
    computations on the other

    >>> from sympy.stats import Die, independent, Given
    >>> from sympy import Tuple

    >>> X, Y = Die(6), Die(6)
    >>> independent(X, Y)
    True
    >>> independent(2*X + Y, -Y)
    False
    >>> X, Y = Given(Tuple(X, Y), X>Y)
    >>> independent(X, Y)
    False

    See Also:
        dependent
    """
    return not dependent(a, b)

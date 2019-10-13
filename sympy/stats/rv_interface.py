from __future__ import print_function, division

from sympy import sqrt, log, exp, FallingFactorial
from .rv import (probability, expectation, density, where, given, pspace, cdf,
                 characteristic_function, sample, sample_iter, random_symbols, independent, dependent,
                 sampling_density, moment_generating_function, quantile)

__all__ = ['P', 'E', 'H', 'density', 'where', 'given', 'sample', 'cdf', 'characteristic_function', 'pspace',
        'sample_iter', 'variance', 'std', 'skewness', 'kurtosis', 'covariance',
        'dependent', 'independent', 'random_symbols', 'correlation', 'factorial_moment',
        'moment', 'cmoment', 'sampling_density', 'moment_generating_function', 'quantile']



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

def entropy(expr, condition=None, **kwargs):
    """
    Calculuates entropy of a probability distribution

    Parameters
    ==========

    expression : the random expression whose entropy is to be calculated
    condition : optional, to specify conditions on random expression
    b: base of the logarithm, optional
       By default, it is taken as Euler's number

    Returns
    =======

    result : Entropy of the expression, a constant

    Examples
    ========

    >>> from sympy.stats import Normal, Die, entropy
    >>> X = Normal('X', 0, 1)
    >>> entropy(X)
    log(2)/2 + 1/2 + log(pi)/2

    >>> D = Die('D', 4)
    >>> entropy(D)
    log(4)

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Entropy_(information_theory)
    .. [2] https://www.crmarsh.com/static/pdf/Charles_Marsh_Continuous_Entropy.pdf
    .. [3] http://www.math.uconn.edu/~kconrad/blurbs/analysis/entropypost.pdf
    """
    pdf = density(expr, condition, **kwargs)
    base = kwargs.get('b', exp(1))
    if hasattr(pdf, 'dict'):
            return sum([-prob*log(prob, base) for prob in pdf.dict.values()])
    return expectation(-log(pdf(expr), base))

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
    E(((X - mu)/sigma(X))**n)

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
    Measure of the asymmetry of the probability distribution.

    Positive skew indicates that most of the values lie to the right of
    the mean.

    skewness(X) = E(((X - E(X))/sigma)**3)

    Parameters
    ==========

    condition : Expr containing RandomSymbols
            A conditional expression. skewness(X, X>0) is skewness of X given X > 0

    Examples
    ========

    >>> from sympy.stats import skewness, Exponential, Normal
    >>> from sympy import Symbol
    >>> X = Normal('X', 0, 1)
    >>> skewness(X)
    0
    >>> skewness(X, X > 0) # find skewness given X > 0
    (-sqrt(2)/sqrt(pi) + 4*sqrt(2)/pi**(3/2))/(1 - 2/pi)**(3/2)

    >>> rate = Symbol('lambda', positive=True, real=True, finite=True)
    >>> Y = Exponential('Y', rate)
    >>> skewness(Y)
    2
    """
    return smoment(X, 3, condition=condition, **kwargs)

def kurtosis(X, condition=None, **kwargs):
    """
    Characterizes the tails/outliers of a probability distribution.

    Kurtosis of any univariate normal distribution is 3. Kurtosis less than
    3 means that the distribution produces fewer and less extreme outliers
    than the normal distribution.

    kurtosis(X) = E(((X - E(X))/sigma)**4)

    Parameters
    ==========

    condition : Expr containing RandomSymbols
            A conditional expression. kurtosis(X, X>0) is kurtosis of X given X > 0

    Examples
    ========

    >>> from sympy.stats import kurtosis, Exponential, Normal
    >>> from sympy import Symbol
    >>> X = Normal('X', 0, 1)
    >>> kurtosis(X)
    3
    >>> kurtosis(X, X > 0) # find kurtosis given X > 0
    (-4/pi - 12/pi**2 + 3)/(1 - 2/pi)**2

    >>> rate = Symbol('lamda', positive=True, real=True, finite=True)
    >>> Y = Exponential('Y', rate)
    >>> kurtosis(Y)
    9

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Kurtosis
    .. [2] http://mathworld.wolfram.com/Kurtosis.html
    """
    return smoment(X, 4, condition=condition, **kwargs)


def factorial_moment(X, n, condition=None, **kwargs):
    """
    The factorial moment is a mathematical quantity defined as the expectation
    or average of the falling factorial of a random variable.

    factorial_moment(X, n) = E(X*(X - 1)*(X - 2)*...*(X - n + 1))

    Parameters
    ==========

    n: A natural number, n-th factorial moment.

    condition : Expr containing RandomSymbols
            A conditional expression.

    Examples
    ========

    >>> from sympy.stats import factorial_moment, Poisson, Binomial
    >>> from sympy import Symbol, S
    >>> lamda = Symbol('lamda')
    >>> X = Poisson('X', lamda)
    >>> factorial_moment(X, 2)
    lamda**2
    >>> Y = Binomial('Y', 2, S.Half)
    >>> factorial_moment(Y, 2)
    1/2
    >>> factorial_moment(Y, 2, Y > 1) # find factorial moment for Y > 1
    2

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Factorial_moment
    .. [2] http://mathworld.wolfram.com/FactorialMoment.html
    """
    return expectation(FallingFactorial(X, n), condition=condition, **kwargs)


P = probability
E = expectation
H = entropy

"""
Continuous Random Variables - Prebuilt variables

Contains
========
Normal
LogNormal
Exponential
Uniform
Pareto
Weibull
Beta
Gamma
"""

from sympy import (exp, log, sqrt, pi, S, Dummy, Interval, S, sympify, gamma,
    Piecewise)
from sympy import beta as beta_fn
from crv import SingleContinuousPSpace
from sympy.core.decorators import _sympifyit
import random

oo = S.Infinity

__all__ = ['ContinuousRV', 'Normal', 'LogNormal', 'Exponential', 'Pareto',
    'Weibull', 'Beta', 'Gamma', 'Uniform']

def _value_check(condition, message):
    """
    Check a condition on input value.

    Raises ValueError with message if condition is not True
    """
    if condition is not True:
        raise ValueError(message)


def ContinuousRV(symbol, density, set=Interval(-oo,oo)):
    """
    Create a Continuous Random Variable given the following:

    -- a symbol
    -- a probability density function
    -- set on which the pdf is valid (defaults to entire real line)

    Returns a RandomSymbol.

    Many common continuous random variable types are already implemented.
    This function should be necessary only very rarely.

    Examples
    ========

    >>> from sympy import Symbol, sqrt, exp, pi
    >>> from sympy.stats import ContinuousRV, P, E

    >>> x = Symbol('x')
    >>> pdf = sqrt(2)*exp(-x**2/2)/(2*sqrt(pi)) # Normal distribution
    >>> X = ContinuousRV(x, pdf)

    >>> E(X)
    0
    >>> P(X>0)
    1/2
    """
    return SingleContinuousPSpace(symbol, density, set).value

class NormalPSpace(SingleContinuousPSpace):
    def __new__(cls, mean, std, symbol=None):
        _value_check(std > 0, "Standard deviation must be positive")

        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = exp(-(x-mean)**2 / (2*std**2)) / (sqrt(2*pi)*std)
        obj = SingleContinuousPSpace.__new__(cls, x, pdf)
        obj.mean = mean
        obj.std = std
        obj.variance = std**2
        return obj

    def sample(self):
        return {self.value: random.normalvariate(self.mean, self.std)}

def Normal(mean, std, symbol=None):
    """
    Create a Continuous Random Variable with a Normal distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Normal, Density, E, Std
    >>> from sympy import Symbol, simplify

    >>> X = Normal(0, 1, symbol=Symbol('x')) # Mean 0, standard deviation 1
    >>> Density(X)
    (x, sqrt(2)*exp(-x**2/2)/(2*sqrt(pi)))

    >>> E(2*X + 1)
    1

    >>> simplify(Std(2*X + 1))
    2
    """

    return NormalPSpace(mean, std, symbol).value

class LogNormalPSpace(SingleContinuousPSpace):
    def __new__(cls, mean, std, symbol=None):
        mean, std = sympify(mean), sympify(std)

        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = exp(-(log(x)-mean)**2 / (2*std**2)) / (x*sqrt(2*pi)*std)
        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set=Interval(0, oo))
        obj.mean = mean
        obj.std = std
        return obj

    def sample(self):
        return {self.value: random.lognormvariate(self.mean, self.std)}

def LogNormal(mean, std, symbol=None):
    """
    Create a Continuous Random Variable with a LogNormal distribution.

    Note: Only density and sampling work.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import LogNormal, Density, E, Std
    >>> from sympy import Symbol, simplify

    >>> X = LogNormal(0, 1, symbol=Symbol('x')) # Mean 0, standard deviation 1
    >>> Density(X)
    (x, sqrt(2)*exp(-log(x)**2/2)/(2*sqrt(pi)*x))
    """

    return LogNormalPSpace(mean, std, symbol).value


class ExponentialPSpace(SingleContinuousPSpace):
    def __new__(cls, rate, symbol=None):
        _value_check(rate > 0, "Rate must be positive.")

        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = rate * exp(-rate*x)
        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set=Interval(0, oo))
        obj.rate = rate
        return obj

    def sample(self):
        return {self.value: random.expovariate(self.rate)}

def Exponential(rate, symbol=None):
    """
    Create a Continuous Random Variable with an Exponential distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Exponential, Density, E, Std
    >>> from sympy import Symbol

    >>> X = Exponential(rate=10, symbol=Symbol('x')) # Decay rate equals 10
    >>> Density(X)
    (x, 10*exp(-10*x))

    >>> E(X)
    1/10

    >>> Std(X)
    1/10
    """

    return ExponentialPSpace(rate, symbol).value

class ParetoPSpace(SingleContinuousPSpace):
    def __new__(cls, xm, alpha, symbol=None):
        _value_check(xm > 0, "Xm must be positive")
        _value_check(alpha > 0, "Alpha must be positive")

        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = alpha * xm**alpha / x**(alpha+1)
        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set=Interval(xm, oo))
        obj.xm = xm
        obj.alpha = alpha
        return obj

    def sample(self):
        return {self.value: random.paretovariate(self.alpha)}

def Pareto(xm, alpha, symbol=None):
    """
    Create a Continuous Random Variable with the Pareto distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Pareto, Density, E, Std
    >>> from sympy import symbols

    >>> x, xm, beta = symbols('x xm beta', positive=True)
    >>> X = Pareto(xm, beta, symbol=x)
    >>> Density(X)
    (x, beta*x**(-beta - 1)*xm**beta)
    """

    return ParetoPSpace(xm, alpha, symbol).value

class WeibullPSpace(SingleContinuousPSpace):
    def __new__(cls, alpha, beta, symbol=None):
        _value_check(alpha > 0, "Alpha must be positive")
        _value_check(beta > 0, "Beta must be positive")

        alpha, beta = sympify(alpha), sympify(beta)

        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = beta * (x/alpha)**(beta-1) * exp(-(x/alpha)**beta) / alpha

        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set=Interval(0, oo))
        obj.alpha = alpha
        obj.beta = beta
        return obj

    def sample(self):
        return {self.value: random.weibullvariate(self.alpha, self.beta)}

def Weibull(alpha, beta, symbol=None):
    """
    Create a Continuous Random Variable with a Weibull distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Weibull, Density, E, Var
    >>> from sympy import symbols, simplify
    >>> x, a, b = symbols('x a b', positive=True)

    >>> X = Weibull(a, b, symbol=x)
    >>> Density(X)
    (x, b*(x/a)**(b - 1)*exp(-(x/a)**b)/a)
    >>> simplify(E(X))
    a*gamma(1 + 1/b)
    >>> simplify(Var(X))
    -a**2*(gamma(1 + 1/b)**2 - gamma(1 + 2/b))
    """
    return WeibullPSpace(alpha, beta, symbol).value

class BetaPSpace(SingleContinuousPSpace):
    def __new__(cls, alpha, beta, symbol=None):
        _value_check(alpha > 0, "Alpha must be positive")
        _value_check(beta > 0, "Beta must be positive")

        alpha, beta = sympify(alpha), sympify(beta)

        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = x**(alpha-1) * (1-x)**(beta-1) / beta_fn(alpha, beta)

        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set=Interval(0, 1))
        obj.alpha = alpha
        obj.beta = beta
        return obj

    def sample(self):
        return {self.value: random.betavariate(self.alpha, self.beta)}

def Beta(alpha, beta, symbol=None):
    """
    Create a Continuous Random Variable with a Beta distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Beta, Density, E, Std
    >>> from sympy import symbols
    >>> x, a, b = symbols('x a b', positive=True)

    >>> X = Beta(a, b, symbol=x)
    >>> Density(X)
    (x, x**(a - 1)*(-x + 1)**(b - 1)*gamma(a + b)/(gamma(a)*gamma(b)))
    """

    return BetaPSpace(alpha, beta, symbol).value

class GammaPSpace(SingleContinuousPSpace):
    def __new__(cls, k, theta, symbol=None):
        _value_check(k > 0, "k must be positive")
        _value_check(theta > 0, "Theta must be positive")

        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = x**(k-1) * exp(-x/theta) / (gamma(k)*theta**k)

        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set=Interval(0, oo))
        obj.k = k
        obj.theta = theta
        return obj

    def sample(self):
        return {self.value: random.gammavariate(self.k, self.theta)}

def Gamma(k, theta, symbol=None):
    """
    Create a Continuous Random Variable with a Gamma distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Gamma, Density, E, Std
    >>> from sympy import symbols
    >>> x, k, theta = symbols('x k theta', positive=True)

    >>> X = Gamma(k, theta, symbol=x)
    >>> Density(X)
    (x, theta**(-k)*x**(k - 1)*exp(-x/theta)/gamma(k))

    >>> E(X)
    theta*gamma(k + 1)/gamma(k)
    """

    return GammaPSpace(k, theta, symbol).value

class UniformPSpace(SingleContinuousPSpace):
    def __new__(cls, left, right, symbol=None):
        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = Piecewise(
                (S.Zero, x<left),
                (S.Zero, x>right),
                (S.One/(right-left), True))

        obj = SingleContinuousPSpace.__new__(cls, x, pdf)
        obj.left = left
        obj.right = right
        return obj

    def sample(self):
        return {self.value: random.uniform(self.left, self.right)}

def Uniform(left, right, symbol=None):
    """
    Create a Continuous Random Variable with a Uniform distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from sympy.stats import Uniform, Density, E, Var
    >>> from sympy import symbols, simplify
    >>> x, l, r = symbols('x l r')

    >>> X = Uniform(l, r, symbol=x)

    >>> Density(X)
    (x, Piecewise((0, x < l), (0, x > r), (1/(-l + r), True)))

    >>> simplify(E(X))
    l/2 + r/2

    >>> simplify(Var(X))
    l**2/12 - l*r/6 + r**2/12
    """

    return UniformPSpace(left, right, symbol).value

"""
Continuous Random Variables - Prebuilt variables

Contains
========
Normal
Exponential
Uniform
Pareto
Beta
Gamma
"""

from sympy import exp, sqrt, pi, S, Dummy, Interval, S, sympify, gamma, Piecewise
from sympy import beta as beta_fn
from crv import SingleContinuousPSpace
from sympy.core.decorators import _sympifyit
import random

oo = S.Infinity

__all__ = ['ContinuousRV', 'Normal', 'Exponential', 'Pareto', 'Beta', 'Gamma',
'Uniform']

def ContinuousRV(symbol, density, set=Interval(-oo,oo)):
    """
    Create a Continuous Random Variable given
    -- a symbol
    -- a probability density function
    -- set on which the pdf is valid (defaults to entire real line)

    Returns a RandomSymbol

    Many common continuous random variable types are already implemented.
    This function should be necessary only very rarely.

    >>> from sympy import Symbol
    >>> from sympy.stats import ContinuousRV, P, E

    >>> x = Symbol('x')
    >>> pdf = sqrt(2)*exp(-x**2/2)/(2*sqrt(pi)) # Normal distribution
    >>> X = ContinuousRV(x, density)

    >>> E(X)
    0
    >>> P(X>0)
    1/2
    """
    return SingleContinuousPSpace(symbol, density, set).value

class NormalPSpace(SingleContinuousPSpace):
    def __new__(cls, mean, std, symbol = None):

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
    Create a Continuous Random Varible with a Normal distribution
    Returns a RandomSymbol

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

class ExponentialPSpace(SingleContinuousPSpace):
    def __new__(cls, rate, symbol=None):
        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = rate * exp(-rate*x)
        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set = Interval(0, oo))
        obj.rate = rate
        return obj

    def sample(self):
        return {self.value: random.expovariate(self.rate)}

def Exponential(rate, symbol=None):
    """
    Create a Continuous Random Varible with an Exponential distribution
    Returns a RandomSymbol

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
        assert xm>0, "Xm must be positive"
        assert alpha>0, "Alpha must be positive"

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
    Create a Continuous Random Varible with the Pareto distribution
    Returns a RandomSymbol

    >>> from sympy.stats import Pareto, Density, E, Std
    >>> from sympy import symbols

    >>> x, xm, beta = symbols('x xm beta', positive=True)
    >>> X = Pareto(xm, beta, symbol=x)
    >>> Density(X)
    (x, beta*x**(-beta - 1)*xm**beta)
    """

    return ParetoPSpace(xm, alpha, symbol).value

class BetaPSpace(SingleContinuousPSpace):
    def __new__(cls, alpha, beta, symbol=None):
        assert alpha>0, "Alpha must be positive"
        assert beta>0, "Beta must be positive"

        alpha, beta = sympify(alpha), sympify(beta)

        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = x**(alpha-1) * (1-x)**(beta-1) / beta_fn(alpha, beta)

        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set = Interval(0, 1))
        obj.alpha = alpha
        obj.beta = beta
        return obj

    def sample(self):
        return {self.value: random.betavariate(self.alpha, self.beta)}

def Beta(alpha, beta, symbol=None):
    """
    Create a Continuous Random Varible with a Beta distribution
    Returns a RandomSymbol

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
        assert k>0, "k must be positive"
        assert theta>0, "theta must be positive"

        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = x**(k-1) * exp(-x/theta) / (gamma(k)*theta**k)

        obj = SingleContinuousPSpace.__new__(cls, x, pdf, set = Interval(0, oo))
        obj.k = k
        obj.theta = theta
        return obj

def Gamma(k, theta, symbol=None):
    """
    Create a Continuous Random Varible with a Gamma distribution
    Returns a RandomSymbol

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
    Create a Continuous Random Varible with a Uniform distribution
    Returns a RandomSymbol

    >>> from sympy.stats import Uniform, Density, E, Var
    >>> from sympy import symbols, simplify
    >>> x, l, r = symbols('x l r')

    >>> X = Uniform(l, r, symbol=x)

    >>> Density(X)
    (x, Piecewise((0, x < l), (0, r < x), (1/(-l + r), True)))

    >>> simplify(E(X))
    l/2 + r/2

    >>> simplify(Var(X))
    l**2/12 - l*r/6 + r**2/12
    """

    return UniformPSpace(left, right, symbol).value

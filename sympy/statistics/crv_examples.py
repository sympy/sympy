from sympy import exp, sqrt, pi, S, Dummy, Interval, S, sympify, gamma, Piecewise
from sympy import beta as beta_fn
from crv import SingleContinuousPSpace, integrate
from sympy.core.decorators import _sympifyit

oo = S.Infinity

class NormalPSpace(SingleContinuousPSpace):
    def __new__(cls, mean, std, symbol = None):

        x = symbol or SingleContinuousPSpace.create_symbol()
        pdf = exp(-(x-mean)**2 / (2*std**2)) / (sqrt(2*pi)*std)
        obj = SingleContinuousPSpace.__new__(cls, x, pdf)
        obj.mean = mean
        obj.std = std
        obj.variance = std**2
        return obj

def Normal(mean, std, symbol=None):
    """
    Create a Continuous Random Varible with a Normal Distribution
    Returns a RandomSymbol

    >>> from sympy.statistics import Normal, Density, E, std
    >>> from sympy import Symbol

    >>> X = Normal(0, 1, symbol=Symbol('x')) # Mean 0, standard deviation 1
    >>> Density(X)
    (x, 2**(1/2)*exp(-x**2/2)/(2*pi**(1/2)))

    >>> E(2*X + 1)
    1

    >>> std(2*X + 1)
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

def Exponential(rate, symbol=None):
    """
    Create a Continuous Random Varible with an Exponential Distribution
    Returns a RandomSymbol

    >>> from sympy.statistics import Exponential, Density, E, std
    >>> from sympy import Symbol

    >>> X = Exponential(rate=10, symbol=Symbol('x')) # Decay rate equals 10
    >>> Density(X)
    (x, 10*exp(-10*x))

    >>> E(X)
    1/10

    >>> std(X)
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

def Pareto(xm, alpha, symbol=None):
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

def Beta(alpha, beta, symbol=None):
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

def Uniform(left, right, symbol=None):
    return UniformPSpace(left, right, symbol).value

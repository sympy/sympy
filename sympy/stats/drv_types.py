from sympy.stats.drv import SingleDiscreteDistribution, SingleDiscretePSpace
from sympy import factorial, exp, Basic, Range, S, oo, sympify
from sympy.stats.rv import _value_check

def rv(symbol, cls, *args):
    args = map(sympify, args)
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
    return rv(name, PoissonDistribution, lamda)

class GeometricDistribution(SingleDiscreteDistribution):
    _argnames = ('p',)
    set = S.Naturals

    @staticmethod
    def check(p):
        _value_check(0 < p <= 1, "p must be between 0 and 1")

    def pdf(self, k):
        return (1 - self.p)**(k - 1) * self.p

def Geometric(name, p):
    return rv(name, GeometricDistribution, p)

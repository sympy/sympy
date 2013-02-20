from sympy.stats.drv import SingleDiscreteDistribution
from sympy import factorial, exp, Basic, Range, S, oo

class PoissonDistribution(SingleDiscreteDistribution):
    _argnames = ('lamda',)

    set = S.Naturals0

    @staticmethod
    def check(mean, std):
        _value_check(lamda > 0, "Standard deviation must be positive")

    def pdf(self, k):
        return self.lamda**k / factorial(k) * exp(-self.lamda)

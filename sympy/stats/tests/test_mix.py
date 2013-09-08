from sympy.stats import Poisson, Beta
from sympy.stats.rv import pspace, ProductPSpace, density
from sympy.stats.drv_types import PoissonDistribution
from sympy import Symbol, Eq

def test_density():
    x = Symbol('x')
    l = Symbol('l', positive=True)
    rate = Beta(l, 2, 3)
    X = Poisson(x, rate)
    assert isinstance(pspace(X), ProductPSpace)
    assert density(X, Eq(rate, rate.symbol)) == PoissonDistribution(l)

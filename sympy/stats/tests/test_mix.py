from sympy import Eq, Symbol
from sympy.stats import Beta, Poisson
from sympy.stats.drv_types import PoissonDistribution
from sympy.stats.rv import ProductPSpace, density, pspace


def test_density():
    x = Symbol('x')
    l = Symbol('l', positive=True)
    rate = Beta(l, 2, 3)
    X = Poisson(x, rate)
    assert isinstance(pspace(X), ProductPSpace)
    assert density(X, Eq(rate, rate.symbol)) == PoissonDistribution(l)

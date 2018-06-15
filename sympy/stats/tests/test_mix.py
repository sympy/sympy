from sympy.stats import Poisson, Beta, Exponential, P
from sympy.stats.rv import pspace, ProductPSpace, density
from sympy.stats.drv_types import PoissonDistribution
from sympy import Symbol, Eq, Ne

def test_density():
    x = Symbol('x')
    l = Symbol('l', positive=True)
    rate = Beta(l, 2, 3)
    X = Poisson(x, rate)
    assert isinstance(pspace(X), ProductPSpace)
    assert density(X, Eq(rate, rate.symbol)) == PoissonDistribution(l)

def test_mix_expression():
    Y, E = Poisson('Y', 1), Exponential('E', 1)
    assert P(Eq(Y + E, 1)) == 0
    assert P(Ne(Y + E, 2)) == 1
    assert str(P(E + Y < 2, evaluate=False)) == """Integral(Sum(exp(-1)*Integral"""\
+"""(exp(-E)*DiracDelta(-_z + E + Y - 2), (E, 0, oo))/factorial(Y), (Y, 0, oo)), (_z, -oo, 0))"""
    assert str(P(E + Y > 2, evaluate=False)) == """Integral(Sum(exp(-1)*Integral"""\
+"""(exp(-E)*DiracDelta(-_z + E + Y - 2), (E, 0, oo))/factorial(Y), (Y, 0, oo)), (_z, 0, oo))"""

from sympy.statistics import NormalPSpace, P, E, Where, Density, var
from sympy import Symbol, exp, S, pi, simplify

def test_single_normal():
    mu = Symbol('mu', real=True, bounded=True, finite=True)
    sigma = Symbol('sigma', real=True, positive=True, bounded=True)
    ps = NormalPSpace(0,1)
    X = ps.value

    assert E(mu + sigma*X) == mu
    assert simplify(var(mu + sigma*X)) == sigma**2
    x, pdf = Density(X*sigma + mu)
    assert pdf == 2**S.Half*exp(-(x - mu)**2/(2*sigma**2))/(2*pi**S.Half*sigma)

def test_multiple_normal():
    X, Y = NormalPSpace(0,1).value, NormalPSpace(0,1).value

    assert E(X+Y) == 0
    assert var(X+Y) == 2
    assert var(X+X) == 4



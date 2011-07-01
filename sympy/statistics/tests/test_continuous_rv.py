from sympy.statistics import NormalPSpace, P, E, Where, Density, var
from sympy import Symbol, exp, S, pi, simplify, Interval, erf

def test_single_normal():
    mu = Symbol('mu', real=True, bounded=True, finite=True)
    sigma = Symbol('sigma', real=True, positive=True, bounded=True)
    ps = NormalPSpace(0,1)
    X = ps.value
    Y = X*sigma + mu

    assert E(Y) == mu
    assert simplify(var(Y)) == sigma**2
    x, pdf = Density(Y)
    assert pdf == 2**S.Half*exp(-(x - mu)**2/(2*sigma**2))/(2*pi**S.Half*sigma)

    assert Where(X**2<=1).set == Interval(-1,1)
    assert Where(X**2<=1).symbol == X.symbol

    assert P(X**2<1) == erf(2**S.Half/2)




def test_multiple_normal():
    X, Y = NormalPSpace(0,1).value, NormalPSpace(0,1).value

    assert E(X+Y) == 0
    assert var(X+Y) == 2
    assert var(X+X) == 4



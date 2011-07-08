from sympy.statistics import (Normal, Exponential, P, E, Where,
        Density, var, covar, skewness)
from sympy import Symbol, exp, S, pi, simplify, Interval, erf, Eq

def test_single_normal():
    mu = Symbol('mu', real=True, bounded=True, finite=True)
    sigma = Symbol('sigma', real=True, positive=True, bounded=True)
    X = Normal(0,1)
    Y = X*sigma + mu

    assert E(Y) == mu
    assert simplify(var(Y)) == sigma**2
    x, pdf = Density(Y)
    assert pdf == 2**S.Half*exp(-(x - mu)**2/(2*sigma**2))/(2*pi**S.Half*sigma)

    assert Where(X**2<=1).set == Interval(-1,1)
    assert Where(X**2<=1).symbol == X.symbol

    assert P(X**2<1) == erf(2**S.Half/2)

    assert E(X, Eq(X, mu)) == mu

def test_multiple_normal():
    X, Y = Normal(0,1), Normal(0,1)

    assert E(X+Y) == 0
    assert var(X+Y) == 2
    assert var(X+X) == 4
    assert covar(X, Y) == 0
    assert covar(2*X + Y, -X) == -2*var(X)

    assert E(X, Eq(X+Y, 0)) == 0
    assert var(X, Eq(X+Y, 0)) == S.Half

def test_exponential():

    rate = Symbol('lambda', positive=True, real=True, finite=True)
    X = Exponential(rate)

    assert E(X) == 1/rate
    assert var(X) == 1/rate**2
    assert skewness(X) == 2
    assert P(X>0) == S(1)
    assert P(X>1) == exp(-rate)
    assert P(X>10) == exp(-10*rate)

    assert Where(X<=1).set == Interval(0,1)



from sympy.stats import (Normal, Exponential, P, E, Where,
        Density, var, covar, skewness, Gamma, Pareto, Beta, Given, pspace, CDF)
from sympy import (Symbol, exp, S, pi, simplify, Interval, erf, Eq, symbols,
        sqrt, And, gamma, beta)
from sympy.utilities.pytest import raises

oo = S.Infinity

def test_single_normal():
    mu = Symbol('mu', real=True, bounded=True, finite=True)
    sigma = Symbol('sigma', real=True, positive=True, bounded=True)
    X = Normal(0,1)
    Y = X*sigma + mu

    assert simplify(E(Y)) == mu
    assert simplify(var(Y)) == sigma**2
    x, pdf = Density(Y)
    assert pdf == 2**S.Half*exp(-(x - mu)**2/(2*sigma**2))/(2*pi**S.Half*sigma)

    assert P(X**2<1) == erf(2**S.Half/2)

    assert E(X, Eq(X, mu)) == mu

def test_conditional_1d():
    X = Normal(0,1)
    Y = Given(X, X>=0)

    assert Density(Y)[1] == 2 * Density(X)[1]

    assert Y.pspace.domain.set == Interval(0, oo)
    assert E(Y) == sqrt(2) / sqrt(pi)

    assert E(X**2) == E(Y**2)

def test_ContinuousDomain():
    X = Normal(0,1)
    assert Where(X**2<=1).set == Interval(-1,1)
    assert Where(X**2<=1).symbol == X.symbol
    Where(And(X**2<=1, X>=0)).set == Interval(0,1)

    Y = Given(X, X>=0)

    assert Y.pspace.domain.set == Interval(0, oo)





def test_multiple_normal():
    X, Y = Normal(0,1), Normal(0,1)

    assert E(X+Y) == 0
    assert var(X+Y) == 2
    assert var(X+X) == 4
    assert covar(X, Y) == 0
    assert covar(2*X + Y, -X) == -2*var(X)

    assert E(X, Eq(X+Y, 0)) == 0
    assert var(X, Eq(X+Y, 0)) == S.Half

def test_symbolic():
    mu1, mu2 = symbols('mu1 mu2', real=True, finite=True, bounded=True)
    s1, s2 = symbols('sigma1 sigma2', real=True, finite=True, positive=True)
    rate = Symbol('lambda', real=True, positive=True, bounded=True)
    X = Normal(mu1, s1)
    Y = Normal(mu2, s2)
    Z = Exponential(rate)
    a, b, c = symbols('a b c', real=True, finite=True)

    assert E(X) == mu1
    assert E(X+Y) == mu1+mu2
    assert E(a*X+b) == a*E(X)+b
    assert var(X) == s1**2
    assert simplify(var(X+a*Y+b)) == var(X) + a**2*var(Y)

    assert E(Z) == 1/rate
    assert E(a*Z+b) == a*E(Z)+b
    assert E(X+a*Z+b) == mu1 + a/rate + b

def test_CDF():
    X = Normal(0,1)
    x,d = CDF(X)

    assert P(X<1) == d.subs(x,1)
    assert d.subs(x,0) == S.Half

    x,d = CDF(X, X>0) # given X>0

    assert d.subs(x,0) == 0

    Y = Exponential(10)
    y,d = CDF(Y)

    assert d.subs(y, -5) == 0
    assert P(Y>3) == 1 - d.subs(y, 3)

    raises(ValueError, "CDF(X+Y)")

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

def test_pareto():

    xm, beta = symbols('xm beta', real=True, positive=True)
    alpha = beta + 5
    X = Pareto(xm, alpha)

    assert simplify(E(X)) == alpha*xm/(alpha-1)
    assert simplify(var(X)) == xm**2*alpha / ((alpha-1)**2*(alpha-2))

def test_gamma():
    k, theta = symbols('k theta', real=True, finite=True, positive=True)
    X = Gamma(k, theta)

    assert simplify(E(X)) == k*theta
    # can't get things to simplify on this one so we use subs
    assert var(X).subs(k,5) == (k*theta**2).subs(k, 5)
    assert simplify(skewness(X)).subs(k, 5) == (2/sqrt(k)).subs(k, 5)

def test_beta():
    a, b = symbols('alpha beta', positive=True)

    B = Beta(a, b)

    assert pspace(B).domain.set == Interval(0, 1)

    x, dens = Density(B)
    assert dens == x**(a-1)*(1-x)**(b-1) / beta(a,b)

    assert E(B) == a / (a + b)
    assert var(B) == (a*b) / ((a+b)**2 * (a+b+1))

def test_uniform():
    l = Symbol('l', real=True, bounded=True)
    w = Symbol('w', positive=True, bounded=True)
    X = Uniform(l, l+w)

    assert simplify(E(X)) == l + w/2
    assert simplify(E(X)) == w**2/12

    assert P(X<l) == 0 and P(X>l+w) == 0

from sympy.stats import (Normal, LogNormal, Exponential, P, E, Where, Density,
        Var, Covar, Skewness, Gamma, Pareto, Weibull, Beta, Uniform, Given, pspace, CDF, ContinuousRV, Sample)
from sympy import (Symbol, exp, S, N, pi, simplify, Interval, erf, Eq, symbols,
        sqrt, And, gamma, beta, Piecewise)
from sympy.utilities.pytest import raises

oo = S.Infinity

def test_single_normal():
    mu = Symbol('mu', real=True, bounded=True)
    sigma = Symbol('sigma', real=True, positive=True, bounded=True)
    X = Normal(0,1)
    Y = X*sigma + mu

    assert simplify(E(Y)) == mu
    assert simplify(Var(Y)) == sigma**2
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
    assert Var(X+Y) == 2
    assert Var(X+X) == 4
    assert Covar(X, Y) == 0
    assert Covar(2*X + Y, -X) == -2*Var(X)

    assert E(X, Eq(X+Y, 0)) == 0
    assert Var(X, Eq(X+Y, 0)) == S.Half

def test_symbolic():
    mu1, mu2 = symbols('mu1 mu2', real=True, bounded=True)
    s1, s2 = symbols('sigma1 sigma2', real=True, bounded=True, positive=True)
    rate = Symbol('lambda', real=True, positive=True, bounded=True)
    X = Normal(mu1, s1)
    Y = Normal(mu2, s2)
    Z = Exponential(rate)
    a, b, c = symbols('a b c', real=True, bounded=True)

    assert E(X) == mu1
    assert E(X+Y) == mu1+mu2
    assert E(a*X+b) == a*E(X)+b
    assert Var(X) == s1**2
    assert simplify(Var(X+a*Y+b)) == Var(X) + a**2*Var(Y)

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

    Z = Exponential(1)
    z, cdf = CDF(Z)
    assert cdf == Piecewise((0, z < 0), (1 - exp(-z), True))

def test_sample():
    z = Symbol('z')
    Z = ContinuousRV(z, exp(-z), set=Interval(0,oo))
    assert Sample(Z) in Z.pspace.domain.set
    sym, val = Z.pspace.sample().items()[0]
    assert sym == Z and val in Interval(0, oo)

def test_ContinuousRV():
    x = Symbol('x')
    pdf = sqrt(2)*exp(-x**2/2)/(2*sqrt(pi)) # Normal distribution
    # X and Y should be equivalent
    X = ContinuousRV(x, pdf)
    Y = Normal(0, 1)

    assert Var(X) == Var(Y)
    assert P(X>0) == P(Y>0)

def test_lognormal():
    mean = Symbol('mu', real=True, bounded=True)
    std = Symbol('sigma', positive=True, real=True, bounded=True)
    X = LogNormal(mean, std)
    # The sympy integrator can't do this too well
    #assert E(X) == exp(mean+std**2/2)
    #assert Var(X) == (exp(std**2)-1) * exp(2*mean + std**2)

    # Right now, only density function and sampling works
    # Test sampling: Only e^mean in sample std of 0
    for i in range(3):
        X = LogNormal(i, 0)
        assert S(Sample(X)) == N(exp(i))
    # The sympy integrator can't do this too well
    #assert E(X) ==

def test_exponential():
    rate = Symbol('lambda', positive=True, real=True, bounded=True)
    X = Exponential(rate)

    assert E(X) == 1/rate
    assert Var(X) == 1/rate**2
    assert Skewness(X) == 2
    assert P(X>0) == S(1)
    assert P(X>1) == exp(-rate)
    assert P(X>10) == exp(-10*rate)

    assert Where(X<=1).set == Interval(0,1)

def test_pareto_numeric():
    xm, beta = 3, 2
    alpha = beta + 5
    X = Pareto(xm, alpha)

    assert E(X) == alpha*xm/S(alpha-1)
    assert Var(X) == xm**2*alpha / S(((alpha-1)**2*(alpha-2)))

def test_pareto():

    xm, beta = symbols('xm beta', positive=True, bounded=True)
    alpha = beta + 5
    X = Pareto(xm, alpha)

    x, density = Density(X)
    assert density == x**(-(alpha+1))*xm**(alpha)*(alpha)

    # These fail because SymPy can not deduce that 1/xm != 0
    # assert simplify(E(X)) == alpha*xm/(alpha-1)
    # assert simplify(Var(X)) == xm**2*alpha / ((alpha-1)**2*(alpha-2))

def test_weibull_numeric():
    # Test for integers and rationals
    a = 1
    bvals = [S.Half, 1, S(3)/2, 5]
    for b in bvals:
        X = Weibull(a, b)
        assert simplify(E(X)) == simplify(a * gamma(1 + 1/S(b)))
        assert simplify(Var(X)) == simplify(a**2 * gamma(1 + 2/S(b)) - E(X)**2)
        # Not testing Skew... it's slow with int/frac values > 3/2

def test_weibull():
    a, b = symbols('a b', positive=True)
    X = Weibull(a, b)

    assert simplify(E(X)) == simplify(a * gamma(1 + 1/b))
    assert simplify(Var(X)) == simplify(a**2 * gamma(1 + 2/b) - E(X)**2)
    # Skewness tests too slow. Try shortcutting function?

def test_gamma():
    k, theta = symbols('k theta', real=True, bounded=True, positive=True)
    X = Gamma(k, theta)

    assert simplify(E(X)) == k*theta
    # can't get things to simplify on this one so we use subs
    assert Var(X).subs(k,5) == (k*theta**2).subs(k, 5)
    # The following is too slow
    # assert simplify(Skewness(X)).subs(k, 5) == (2/sqrt(k)).subs(k, 5)

def test_beta():
    a, b = symbols('alpha beta', positive=True)

    B = Beta(a, b)

    assert pspace(B).domain.set == Interval(0, 1)

    x, dens = Density(B)
    assert dens == x**(a-1)*(1-x)**(b-1) / beta(a,b)

    # This is too slow
    # assert E(B) == a / (a + b)
    # assert Var(B) == (a*b) / ((a+b)**2 * (a+b+1))

    # Full symbolic solution is too much, test with numeric version
    a, b = 1, 2
    B = Beta(a, b)
    assert E(B) == a / S(a + b)
    assert Var(B) == (a*b) / S((a+b)**2 * (a+b+1))

def test_uniform():
    l = Symbol('l', real=True, bounded=True)
    w = Symbol('w', positive=True, bounded=True)
    X = Uniform(l, l+w)

    assert simplify(E(X)) == l + w/2
    assert simplify(Var(X)) == w**2/12

    assert P(X<l) == 0 and P(X>l+w) == 0

    # With numbers all is well
    X = Uniform(3, 5)
    assert P(X<3) == 0 and P(X>5) == 0
    assert P(X<4) == P(X>4) == S.Half

def test_prefab_sampling():
    N = Normal(0, 1)
    L = LogNormal(0, 1)
    E = Exponential(1)
    P = Pareto(1, 3)
    W = Weibull(1, 1)
    U = Uniform(0, 1)
    B = Beta(2,5)
    G = Gamma(1,3)

    variables = [N,L,E,P,W,U,B,G]
    niter = 10
    for var in variables:
        for i in xrange(niter):
            assert Sample(var) in var.pspace.domain.set

def test_input_value_assertions():
    a, b = symbols('a b')
    p, q = symbols('p q', positive=True)

    raises(ValueError, "Normal(3, 0)")
    raises(ValueError, "Normal(a, b)")
    Normal(a, p) # No error raised
    raises(ValueError, "Exponential(a)")
    Exponential(p) # No error raised
    for fn_name in ['Pareto', 'Weibull', 'Beta', 'Gamma']:
        raises(ValueError, "%s(a, p)" % fn_name)
        raises(ValueError, "%s(p, a)" % fn_name)
        eval("%s(p, q)" % fn_name) # No error raised

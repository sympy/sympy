from sympy import symbols, S, erf, sqrt, pi, exp, gamma, Interval, oo
from sympy.stats import Normal, P, E, density, Gamma, Poisson, Rayleigh, variance
from sympy.stats.compound_rv import CompoundDistribution, compound_pspace
from sympy.stats.crv_types import NormalDistribution
from sympy.testing.pytest import raises
from sympy.stats.joint_rv_types import MultivariateNormalDistribution

x = symbols('x')
def test_normal_CompoundDist():
    X = Normal('X', 1, 2)
    Y = Normal('X', X, 4)
    assert density(Y)(x).simplify() == sqrt(10)*exp(-x**2/40 + x/20 - S(1)/40)/(20*sqrt(pi))
    assert E(Y) == 1 # it is always equal to mean of X
    assert P(Y > 1) == S(1)/2 # as 1 is the mean
    assert P(Y > 5).simplify() ==  S(1)/2 - erf(sqrt(10)/5)/2
    assert variance(Y) == variance(X) + 4**2 # 2**2 + 4**2
    # https://math.stackexchange.com/questions/1484451/
    # (Contains proof of E and variance computation)

def test_poisson_CompoundDist():
    k, t, y = symbols('k t y', positive=True, real=True)
    G = Gamma('G', k, t)
    D = Poisson('P', G)
    assert density(D)(y).simplify() == t**y*(t + 1)**(-k - y)*gamma(k + y)/(gamma(k)*gamma(y + 1))
    # https://en.wikipedia.org/wiki/Negative_binomial_distribution#Gamma%E2%80%93Poisson_mixture
    assert E(D).simplify() == k*t # mean of NegativeBinomialDistribution

def test_unevaluated_CompoundDist():
    # these tests need to be removed once they work with evaluation as they are currently not
    # evaluated completely in sympy.
    R = Rayleigh('R', 4)
    X = Normal('X', 3, R)
    assert str(density(X)(x).simplify()) == ("Piecewise((exp(3/4 - x/4)/8, 2*Abs(arg(x - 3))"
    " <= pi/2), (sqrt(2)*Integral(exp(-(16*(x - 3)**2 + R**4)/(32*R**2)), "
    "(R, 0, oo))/(32*sqrt(pi)), True))")
    assert str(E(X, evaluate=False)) == \
    ("Integral(Piecewise((X*(-sqrt(pi)*sinh(X/4 - 3/4) + sqrt(pi)*cosh(X/4 - 3/4))/(8*sqrt(pi)), "
    "2*Abs(arg(X - 3)) <= pi/2), (X*Integral(sqrt(2)*exp(-(X - 3)**2/(2*R**2))*exp(-R**2/32)"
    "/(32*sqrt(pi)), (R, 0, oo)), True)), (X, -oo, oo))")

def test_Compound_Distribution():
    X = Normal('X', 2, 4)
    N = NormalDistribution(X, 4)
    C = CompoundDistribution(N)
    assert C.is_Continuous
    assert C.set == Interval(-oo, oo)
    assert C.pdf(x).simplify() == exp(-x**2/64 + x/16 - S(1)/16)/(8*sqrt(pi))

    raises(ValueError, lambda: CompoundDistribution(NormalDistribution(2, 3)))
    M = MultivariateNormalDistribution([1, 2], [[2, 1], [1, 2]])
    raises(NotImplementedError, lambda: CompoundDistribution(M))

def test_compound_pspace():
    X = Normal('X', 2, 4)
    Y = Normal('Y', 3, 6)
    N = NormalDistribution(X, Y)
    raises(ValueError, lambda: compound_pspace('C', N))
    raises(TypeError, lambda: compound_pspace(['C'], N))
    C = CompoundDistribution(N)
    raises(NotImplementedError, lambda: compound_pspace('C', C))

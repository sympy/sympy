from sympy import (symbols, S, erf, sqrt, pi, exp, gamma, Interval, oo, beta,
                    Eq, Piecewise, Integral, Abs, arg, Dummy, Sum, factorial)
from sympy.stats import (Normal, P, E, density, Gamma, Poisson, Rayleigh,
                        variance, Bernoulli, Beta, Uniform, cdf)
from sympy.stats.compound_rv import CompoundDistribution, CompoundPSpace
from sympy.stats.crv_types import NormalDistribution
from sympy.stats.drv_types import PoissonDistribution
from sympy.stats.frv_types import BernoulliDistribution
from sympy.testing.pytest import raises, ignore_warnings
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


def test_bernoulli_CompoundDist():
    X = Beta('X', 1, 2)
    Y = Bernoulli('Y', X)
    assert density(Y).dict == {0: S(2)/3, 1: S(1)/3}
    assert E(Y) == P(Eq(Y, 1)) == S(1)/3
    assert variance(Y) == S(2)/9
    assert cdf(Y) == {0: S(2)/3, 1: 1}

    # test issue 8128
    a = Bernoulli('a', S(1)/2)
    b = Bernoulli('b', a)
    assert density(b).dict == {0: S(1)/2, 1: S(1)/2}
    assert P(b > 0.5) == S(1)/2

    X = Uniform('X', 0, 1)
    Y = Bernoulli('Y', X)
    assert E(Y) == S(1)/2
    assert P(Eq(Y, 1)) == E(Y)


def test_unevaluated_CompoundDist():
    # these tests need to be removed once they work with evaluation as they are currently not
    # evaluated completely in sympy.
    R = Rayleigh('R', 4)
    X = Normal('X', 3, R)
    _k = Dummy('k')
    exprd = Piecewise((exp(S(3)/4 - x/4)/8, 2*Abs(arg(x - 3)) <= pi/2),
    (sqrt(2)*Integral(exp(-(_k**4 + 16*(x - 3)**2)/(32*_k**2)),
    (_k, 0, oo))/(32*sqrt(pi)), True))
    assert (density(X)(x).simplify()).dummy_eq(exprd.simplify())

    expre = Integral(_k*Integral(sqrt(2)*exp(-_k**2/32)*exp(-(_k - 3)**2/(2*_k**2)
    )/(32*sqrt(pi)), (_k, 0, oo)), (_k, -oo, oo))
    with ignore_warnings(UserWarning): ### TODO: Restore tests once warnings are removed
        assert E(X, evaluate=False).rewrite(Integral).dummy_eq(expre)

    X = Poisson('X', 1)
    Y = Poisson('Y', X)
    Z = Poisson('Z', Y)
    exprd = exp(-1)*Sum(exp(-Y)*Y**x*Sum(exp(-X)*X**Y/(factorial(X)*factorial(Y)
                ), (X, 0, oo)), (Y, 0, oo))/factorial(x)
    assert density(Z)(x).simplify() == exprd

    N = Normal('N', 1, 2)
    M = Normal('M', 3, 4)
    D = Normal('D', M, N)
    exprd = Integral(sqrt(2)*exp(-(_k - 1)**2/8)*Integral(exp(-(-_k + x
    )**2/(2*_k**2))*exp(-(_k - 3)**2/32)/(8*pi*_k)
    , (_k, -oo, oo))/(4*sqrt(pi)), (_k, -oo, oo))
    assert density(D, evaluate=False)(x).dummy_eq(exprd)

def test_Compound_Distribution():
    X = Normal('X', 2, 4)
    N = NormalDistribution(X, 4)
    C = CompoundDistribution(N)
    assert C.is_Continuous
    assert C.set == Interval(-oo, oo)
    assert C.pdf(x, evaluate=True).simplify() == exp(-x**2/64 + x/16 - S(1)/16)/(8*sqrt(pi))

    assert not isinstance(CompoundDistribution(NormalDistribution(2, 3)),
                            CompoundDistribution)
    M = MultivariateNormalDistribution([1, 2], [[2, 1], [1, 2]])
    raises(NotImplementedError, lambda: CompoundDistribution(M))

    X = Beta('X', 2, 4)
    B = BernoulliDistribution(X, 1, 0)
    C = CompoundDistribution(B)
    assert C.is_Finite
    assert C.set == {0, 1}
    y = symbols('y', negative=False, integer=True)
    assert C.pdf(y, evaluate=True) == Piecewise((S(1)/(30*beta(2, 4)), Eq(y, 0)),
                (S(1)/(60*beta(2, 4)), Eq(y, 1)), (0, True))

    k, t, z = symbols('k t z', positive=True, real=True)
    G = Gamma('G', k, t)
    X = PoissonDistribution(G)
    C = CompoundDistribution(X)
    assert C.is_Discrete
    assert C.set == S.Naturals0
    assert C.pdf(z, evaluate=True).simplify() == t**z*(t + 1)**(-k - z)*gamma(k \
                    + z)/(gamma(k)*gamma(z + 1))


def test_compound_pspace():
    X = Normal('X', 2, 4)
    Y = Normal('Y', 3, 6)
    assert not isinstance(Y.pspace, CompoundPSpace)
    N = NormalDistribution(1, 2)
    D = PoissonDistribution(3)
    B = BernoulliDistribution(0.2, 1, 0)
    pspace1 = CompoundPSpace('N', N)
    pspace2 = CompoundPSpace('D', D)
    pspace3 = CompoundPSpace('B', B)
    assert not isinstance(pspace1, CompoundPSpace)
    assert not isinstance(pspace2, CompoundPSpace)
    assert not isinstance(pspace3, CompoundPSpace)
    M = MultivariateNormalDistribution([1, 2], [[2, 1], [1, 2]])
    raises(ValueError, lambda: CompoundPSpace('M', M))
    Y = Normal('Y', X, 6)
    assert isinstance(Y.pspace, CompoundPSpace)
    assert Y.pspace.distribution == CompoundDistribution(NormalDistribution(X, 6))
    assert Y.pspace.domain.set == Interval(-oo, oo)

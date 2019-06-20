from sympy import (Symbol, Eq, Ne, simplify, sqrt, exp, pi, symbols,
                Piecewise, factorial, gamma, IndexedBase)
from sympy.stats import (Poisson, Beta, Exponential, P,
                        Multinomial, MultivariateBeta)
from sympy.stats.crv_types import Normal
from sympy.stats.drv_types import PoissonDistribution
from sympy.stats.joint_rv import JointPSpace, CompoundDistribution, MarginalDistribution
from sympy.stats.rv import pspace, density

def test_density():
    x = Symbol('x')
    l = Symbol('l', positive=True)
    rate = Beta(l, 2, 3)
    X = Poisson(x, rate)
    assert isinstance(pspace(X), JointPSpace)
    assert density(X, Eq(rate, rate.symbol)) == PoissonDistribution(l)
    N1 = Normal('N1', 0, 1)
    N2 = Normal('N2', N1, 2)
    assert density(N2)(0).doit() == sqrt(10)/(10*sqrt(pi))
    assert simplify(density(N2, Eq(N1, 1))(x)) == \
        sqrt(2)*exp(-(x - 1)**2/8)/(4*sqrt(pi))

def test_MarginalDistribution():
    a1, p1, p2, p3 = symbols('a1 p1 p2 p3', positive=True)
    C = Multinomial('C', 3, p1, p2, p3)
    B = MultivariateBeta('B', a1, C[0])
    MGR = MarginalDistribution(B, C[0])
    assert str(MGR.pdf(C)) == \
    str(B*Piecewise((6*p1**C[0]*p2**C[1]*p3**C[2]/(factorial(C[0])*factorial(C[1])*
    factorial(C[2])), Eq(C[0] + C[1] + C[2], 3)), (0, True))*gamma(a1 + C[0])*
    B[0]**(a1 - 1)*B[1]**(C[0] - 1)/(gamma(a1)*gamma(C[0])))

def test_compound_distribution():
    Y = Poisson('Y', 1)
    Z = Poisson('Z', Y)
    assert isinstance(pspace(Z), JointPSpace)
    assert isinstance(pspace(Z).distribution, CompoundDistribution)
    assert Z.pspace.distribution.pdf(1).doit() == exp(-2)*exp(exp(-1))

def test_mix_expression():
    Y, E = Poisson('Y', 1), Exponential('E', 1)
    assert P(Eq(Y + E, 1)) == 0
    assert P(Ne(Y + E, 2)) == 1
    assert str(P(E + Y < 2, evaluate=False)) == """Integral(Sum(exp(-1)*Integral"""\
+"""(exp(-E)*DiracDelta(-_z + E + Y - 2), (E, 0, oo))/factorial(Y), (Y, 0, oo)), (_z, -oo, 0))"""
    assert str(P(E + Y > 2, evaluate=False)) == """Integral(Sum(exp(-1)*Integral"""\
+"""(exp(-E)*DiracDelta(-_z + E + Y - 2), (E, 0, oo))/factorial(Y), (Y, 0, oo)), (_z, 0, oo))"""

from sympy.stats.drv_types import (PoissonDistribution, GeometricDistribution,
        Poisson, Geometric)
from sympy.abc import x
from sympy import S, Sum
from sympy.stats import (P, E, variance, density, characteristic_function,
        where)
from sympy.stats.rv import sample
from sympy.core.relational import Eq, Ne
from sympy.functions.elementary.exponential import exp
from sympy.sets.fancysets import Range
from sympy.logic.boolalg import Or

def test_PoissonDistribution():
    l = 3
    p = PoissonDistribution(l)
    assert abs(p.cdf(10).evalf() - 1) < .001
    assert p.expectation(x, x) == l
    assert p.expectation(x**2, x) - p.expectation(x, x)**2 == l

def test_Poisson():
    l = 3
    x = Poisson('x', l)
    assert E(x) == l
    assert variance(x) == l
    assert density(x) == PoissonDistribution(l)
    assert isinstance(E(x, evaluate=False), Sum)
    assert isinstance(E(2*x, evaluate=False), Sum)
    assert characteristic_function(x)(0).doit() == 1

def test_GeometricDistribution():
    p = S.One / 5
    d = GeometricDistribution(p)
    assert d.expectation(x, x) == 1/p
    assert d.expectation(x**2, x) - d.expectation(x, x)**2 == (1-p)/p**2
    assert abs(d.cdf(20000).evalf() - 1) < .001
    assert d.characteristic_function(0).doit() == 1

def test_sample():
    X, Y, Z = Geometric('X', S(1)/2), Poisson('Y', 4), Poisson('Z', 1000)
    W = Poisson('W', S(1)/100)
    assert sample(X) in X.pspace.domain.set
    assert sample(Y) in Y.pspace.domain.set
    assert sample(Z) in Z.pspace.domain.set
    assert sample(W) in W.pspace.domain.set

def test_discrete_probability():
    X = Geometric('X', S(1)/5)
    Y = Poisson('Y', 4)
    G = Geometric('e', x)
    assert P(Eq(X, 3)) == S(16)/125
    assert P(X < 3) == S(9)/25
    assert P(X > 3) == S(64)/125
    assert P(X >= 3) == S(16)/25
    assert P(X <= 3) == S(61)/125
    assert P(Ne(X, 3)) == S(109)/125
    assert P(Eq(Y, 3)) == 32*exp(-4)/3
    assert P(Y < 3) == 13*exp(-4)
    assert P(Y > 3).equals(32*(-S(71)/32 + 3*exp(4)/32)*exp(-4)/3)
    assert P(Y >= 3).equals(32*(-39/32 + 3*exp(4)/32)*exp(-4)/3)
    assert P(Y <= 3) == 71*exp(-4)/3
    assert P(Ne(Y, 3)).equals(
        13*exp(-4) + 32*(-71/32 + 3*exp(4)/32)*exp(-4)/3)
    assert P(X < S.Infinity) is S.One
    assert P(X > S.Infinity) is S.Zero
    assert P(G < 3) == x*(-x + 1) + x
    assert P(Eq(G, 3)) == x*(-x + 1)**2

def test_Or():
    X = Geometric('X', S(1)/2)
    P(Or(X < 3, X > 4)) == S(13)/16
    P(Or(X > 2, X > 1)) == P(X > 1)
    P(Or(X >= 3, X < 3)) == 1

def test_where():
    X = Geometric('X', S(1)/5)
    Y = Poisson('Y', 4)
    assert where(X**2 > 4).set == Range(3, S.Infinity, 1)
    assert where(X**2 >= 4).set == Range(2, S.Infinity, 1)
    assert where(Y**2 < 9).set == Range(0, 3, 1)
    assert where(Y**2 <= 9).set == Range(0, 4, 1)

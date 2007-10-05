from sympy import *
from sympy.numerics import Float
from sympy.functions import erf
from sympy.statistics import *
import operator # XXX weird abs/sympy.abs conflict

def test_normal():
    Float.store()
    Float.setdps(20)
    N = Normal(0, 1)
    assert N.mean == 0
    assert N.variance == 1
    assert N.probability(-1, 1) == erf(1/sqrt(2))
    assert N.probability(-1, 0) == erf(1/sqrt(2))/2
    N = Normal(2, 4)
    assert N.mean == 2
    assert N.variance == 16
    assert N.confidence(1) == (-oo, oo)
    assert N.probability(1, 3) == erf(1/sqrt(32))
    for p in [0.1, 0.3, 0.7, 0.9, 0.995]:
        a, b = N.confidence(p)
        assert operator.abs(float(N.probability(a, b).evalf()) - p) < 1e-10
    Float.revert()

def test_uniform():
    U = Uniform(-3, -1)
    assert U.mean == -2
    assert U.confidence(1) == (-3, -1)
    assert U.confidence(Rational(1,2)) == (Rational(-5,2), Rational(-3,2))
    assert U.pdf(-4) == 0
    assert U.pdf(-Rational(3,2)) == Rational(1,2)
    assert U.pdf(0) == 0
    assert U.cdf(-4) == 0
    assert U.cdf(-Rational(3,2)) == Rational(3,4)
    assert U.cdf(0) == 1

def test_fit():
    import random
    random.seed(1234)
    n = Normal.fit(Uniform.fit(Normal(2, 1.5).random(1000)))
    #print n.mean
    #print n.stddev
    assert abs(n.mean - 2) < 0.3
    assert abs(n.stddev - 1.5) < 0.3

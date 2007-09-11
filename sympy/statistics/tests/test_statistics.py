from sympy import *
from sympy.statistics import *

def test_normal():
    N = Normal(2, 4)
    assert N.mean == 2
    assert N.variance == 16
    assert N.confidence(1) == (-oo, oo)
    for p in [0.1, 0.3, 0.7, 0.9, 0.995]:
        a, b = N.confidence(p)
        assert abs(float(N.probability(a, b).evalf()) - p) < 1e-10

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

from sympy import sqrt, Rational, oo, Symbol, exp, pi
from sympy.functions import erf
from sympy.statistics import Normal, Uniform
from sympy.statistics.distributions import PDF
import operator # XXX weird abs/sympy.abs conflict

from sympy.mpmath import mp

def test_normal():
    dps, mp.dps = mp.dps, 20

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
    assert N.pdf(1).evalf() == (exp(Rational(-1,32)) / (4*sqrt(2*pi))).evalf()
    for p in [0.1, 0.3, 0.7, 0.9, 0.995]:
        a, b = N.confidence(p)
        assert operator.abs(float(N.probability(a, b).evalf()) - p) < 1e-10

    N = Normal(0, 2/sqrt(2*pi))
    assert N.pdf(0) == Rational(1,2)
    mp.dps = dps

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

def test_sample():
    from sympy.statistics.distributions import Sample
    s = Sample([0,1])
    assert s.mean == Rational(1,2)
    assert s.median == Rational(1,2)
    s = Sample([4,2,3])
    assert s == Sample([2, 3, 4])
    assert s.median == 3
    s = Sample([4,2,3,1])
    assert s.median == Rational(5,2)

def test_PDF():
    a = Symbol('a', positive=True)
    x = Symbol('x', real=True)
    exponential = PDF(exp(-x/a), (x,0,oo))
    exponential = exponential.normalize()
    assert exponential.pdf(x) == 1/a*exp(-x/a)
    assert exponential.cdf(x) == 1 - exp(-x/a)
    assert exponential.mean == a
    assert exponential.variance == a**2
    assert exponential.stddev == a

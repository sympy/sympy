from sympy import sqrt, Rational, oo, Symbol, exp, pi, symbols
from sympy.functions import erf

from operator import abs

from sympy.mpmath import mp

from sympy.utilities.tests.test_pickling import check

# Disable sympy.statistics deprecation warning for the tests
# The warning is in __init__.py, so we only need to disable it for imports

import warnings
from sympy.utilities.exceptions import SymPyDeprecationWarning
warnings.filterwarnings("ignore", message="sympy.statistics has been deprecated since SymPy 0.7.2",
    category=SymPyDeprecationWarning)

from sympy.statistics.distributions import (Normal, Uniform, Sample, PDF,
    ContinuousProbability)

warnings.simplefilter("error", category=SymPyDeprecationWarning)

x, y, z = symbols('x y z')

def test_normal():
    dps, mp.dps = mp.dps, 20

    N = Normal(0, 1)
    assert N.random()
    assert N.mean == 0
    assert N.variance == 1
    assert N.probability(-1, 1) == erf(1/sqrt(2))
    assert N.probability(-1, 0) == erf(1/sqrt(2))/2
    N = Normal(2, 4)
    assert N.mean == 2
    assert N.variance == 16
    assert N.confidence(1) == (-oo, oo)
    assert N.probability(1, 3) == erf(1/sqrt(32))
    assert N.pdf(1).evalf() == (exp(Rational(-1, 32)) / (4*sqrt(2*pi))).evalf()
    for p in [0.1, 0.3, 0.7, 0.9, 0.995]:
        a, b = N.confidence(p)
        assert abs(float(N.probability(a, b).evalf()) - p) < 1e-10

    N = Normal(0, 2/sqrt(2*pi))
    assert N.pdf(0) == Rational(1, 2)
    mp.dps = dps


def test_uniform():
    U = Uniform(-3, -1)
    assert str(U) == "Uniform(-3, -1)"
    assert repr(U) == "Uniform(-3, -1)"
    x = U.random()
    assert x < -1 and x > -3
    assert U.mean == -2
    assert U.confidence(1) == (-3, -1)
    assert U.confidence(Rational(1, 2)) == (Rational(-5, 2), Rational(-3, 2))
    assert U.pdf(-4) == 0
    assert U.pdf(-Rational(3, 2)) == Rational(1, 2)
    assert U.pdf(0) == 0
    assert U.cdf(-4) == 0
    assert U.cdf(-Rational(3, 2)) == Rational(3, 4)
    assert U.cdf(0) == 1


def test_fit():
    import random
    random.seed(1234)
    n = Normal.fit(Uniform.fit(Normal(2, 1.5).random(1000)))
    #print n.mean
    #print n.stddev
    assert abs(n.mean - 2) < 0.3
    assert abs(n.stddev - 1.5) < 0.3
    n = Normal.fit([1, 2, 3, 4, 5])
    assert n.mean == 3
    assert n.stddev == sqrt(2)
    n = Uniform.fit([1, 2, 3, 4, 5])
    assert n.mean == 3
    assert n.stddev == sqrt(2)


def test_sample():
    s = Sample([0, 1])
    assert str(s) == "Sample([0, 1])"
    assert repr(s) == "Sample([0, 1])"
    assert s.mean == Rational(1, 2)
    assert s.median == Rational(1, 2)
    s = Sample([4, 2, 3])
    assert s == Sample([2, 3, 4])
    assert s.median == 3
    s = Sample([4, 2, 3, 1])
    assert s.median == Rational(5, 2)


def test_PDF():
    a = Symbol('a', positive=True)
    x = Symbol('x', real=True)
    exponential = PDF(exp(-x/a), (x, 0, oo))
    exponential = exponential.normalize()
    assert exponential.pdf(x) == 1/a*exp(-x/a)
    assert exponential.cdf(x) == 1 - exp(-x/a)
    assert exponential.mean == a
    assert exponential.variance == a**2
    assert exponential.stddev == a
    exponential = PDF(exp(-x/a), x)
    assert exponential.pdf(x) == exp(-x/a)
    assert exponential.cdf(x) == -a*exp(-x/a) + oo
    assert exponential.mean == -oo
    exponential = PDF(1, (x, 1, 2))
    assert exponential.normalize() == exponential
    assert exponential._get_stddev() == sqrt(3)/6
    assert exponential._get_stddev() == sqrt(3)/6
    #This test is intentionally repeated to test PDF._get_stddev() properly.
    exponential = exponential.transform(x, x)
    assert exponential.pdf(x) == 1
    assert exponential.cdf(x) == x - 1

# These last two tests are here instead of test_str.py and test_pickling.py
# because this module is deprecated.

def test_printing():
    assert str(Normal(x + y, z)) == "Normal(x + y, z)"
    assert str(Sample([x, y, 1])) in [
        "Sample([x, y, 1])",
        "Sample([y, 1, x])",
        "Sample([1, x, y])",
        "Sample([y, x, 1])",
        "Sample([x, 1, y])",
        "Sample([1, y, x])",
    ]
    assert str(Uniform(x, y)) == "Uniform(x, y)"
    assert str(Uniform(x + y, y)) == "Uniform(x + y, y)"

def test_pickling():
    for c in (ContinuousProbability, ContinuousProbability(), Normal,
              Normal(x, y), Sample, Sample([1, 3, 4]), Uniform, Uniform(x, y)):
        check(c)

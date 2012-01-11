from sympy import symbols, erf, nan, oo, Float, sqrt, pi, O, I

from sympy.core.function import ArgumentIndexError

from sympy.utilities.pytest import raises

x, y = symbols('x,y')

def test_erf():
    assert erf(nan) == nan

    assert erf(oo) == 1
    assert erf(-oo) == -1

    assert erf(0) == 0

    assert erf(-2) == -erf(2)
    assert erf(-x*y) == -erf(x*y)
    assert erf(-x - y) == -erf(x + y)

    assert erf(I).is_real == False
    assert erf(0).is_real == True

    assert erf(x).as_leading_term(x) == x
    assert erf(1/x).as_leading_term(x) == erf(1/x)

    raises(ArgumentIndexError, 'erf(x).fdiff(2)')

def test_erf_series():
    assert erf(x).series(x, 0, 7) == 2*x/sqrt(pi) - \
        2*x**3/3/sqrt(pi) + x**5/5/sqrt(pi) + O(x**7)

def test_erf_evalf():
    assert abs( erf(Float(2.0)) - 0.995322265 )  <  1E-8  # XXX

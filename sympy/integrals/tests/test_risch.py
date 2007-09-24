
from sympy import *
from sympy.integrals.risch import components

x, y = symbols('xy')

def test_risch_norman_polynomials():
    assert risch_norman(1, x) == x
    assert risch_norman(x, x) == x**2/2
    assert risch_norman(x**17, x) == x**18/18

def test_risch_norman_fractions():
    assert risch_norman(1/x, x) == log(x)
    assert risch_norman(1/(2 + x), x) == log(x + 2)

    assert risch_norman(5*x**5/(2*x**6 + 5), x) == 5*log(5 + 2*x**6) / 12

    assert risch_norman(1/x**2, x) == -1/x
    assert risch_norman(-1/x**5, x) == 1/(4*x**4)

def test_risch_norman_log():
    assert risch_norman(log(x), x) == x*log(x) - x
    assert risch_norman(log(3*x), x) == x*log(3*x) - x
    assert risch_norman(log(x**2), x) == x*log(x**2) - 2*x

def test_risch_norman_exp():
    assert risch_norman(exp(x), x) == exp(x)
    assert risch_norman(exp(-x), x) == -exp(-x)
    assert risch_norman(exp(17*x), x) == exp(17*x) / 17
    assert risch_norman(x*exp(x), x) == x*exp(x) - exp(x)
    assert risch_norman(x*exp(x**2), x) == exp(x**2) / 2

    assert risch_norman(exp(-x**2), x) is None

def test_risch_norman_trigonometric():
    assert risch_norman(sin(x), x) == -cos(x)
    assert risch_norman(cos(x), x) == sin(x)

    assert risch_norman(sin(x)*cos(y), x) == -cos(x)*cos(y)

    assert risch_norman(sin(x)*cos(x), x) == sin(x)**2 / 2
    assert risch_norman(cos(x)/sin(x), x) == log(sin(x))

    assert risch_norman(x*sin(7*x), x) == sin(7*x) / 49 - x*cos(7*x) / 7
    assert risch_norman(x**2*cos(x), x) == x**2*sin(x) - 2*sin(x) + 2*x*cos(x)

def test_risch_norman_hyperbolic():
    assert risch_norman(sinh(x), x) == cosh(x)
    assert risch_norman(cosh(x), x) == sinh(x)

    assert risch_norman(x*sinh(x), x) == x*cosh(x) - sinh(x)
    assert risch_norman(x*cosh(x), x) == x*sinh(x) - cosh(x)

def test_risch_norman_mixed():
    assert risch_norman(sin(x)*exp(x), x) == exp(x)*sin(x)/2 - exp(x)*cos(x)/2

def test_risch_norman_special():
    assert risch_norman(erf(x), x) == x*erf(x) + exp(-x**2)/sqrt(pi)
    assert risch_norman(exp(-x**2)*erf(x), x) == sqrt(pi)*erf(x)**2 / 4

def test_components():
    assert components(x*y) == set([y, x])
    assert components(sin(x)) == set([sin(x), x])
    assert components(sin(x)*cos(x)**2) == set([sin(x), cos(x), x])
    assert components(sin(x)*sqrt(log(x))) == set([sin(x), sqrt(log(x)), log(x), x])
    assert components(x*sin(exp(x)*y)) == set([y, sin(y*exp(x)), x, exp(x)])

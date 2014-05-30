from sympy import (Symbol, Rational, ln, exp, log, sqrt, E, O, pi, I, oo, sinh,
    sin, cosh, cos, tanh, coth, asinh, acosh, atanh, acoth, tan, cot, Integer,
    PoleError, floor, ceiling, asin, symbols, limit, Piecewise, Eq, sign,
    Derivative)
from sympy.abc import x, y, z

from sympy.utilities.pytest import raises, XFAIL

def test_simple():
    e = sin(1/x + exp(-x)) - sin(1/x)
    assert e.aseries(x) == (1/(24*x**4) - 1/(2*x**2) + 1 + O(x**(-6), (x, oo)))*exp(-x)
    assert e.series(x, oo) == (1/(24*x**4) - 1/(2*x**2) + 1 + O(x**(-6), (x, oo)))*exp(-x)

    e = exp(exp(x)) * (exp(sin(1/x + 1/exp(exp(x)))) - exp(sin(1/x)))
    assert e.aseries(x, n=4) == -1/(2*x**3) + 1/x + 1 + O(x**(-4), (x, oo))
    e = exp(x) * (exp(1/x + exp(-x)) - exp(1/x))
    assert e.aseries(x, n=4) == 1/(6*x**3) + 1/(2*x**2) + 1/x + 1 + O(x**(-4), (x, oo))
    e = exp(sin(1/x + exp(-exp(x)))) - exp(sin(1/x))
    assert e.aseries(x, n=4) == (-1/(2*x**3) + 1/x + 1 + O(x**(-4), (x, oo)))*exp(-exp(x))

    e = exp(exp(x)/(1 - 1/x))
    assert e.aseries(x) == exp(exp(x)/(1 - 1/x))
    assert e.aseries(x, bound=3) == exp(exp(x)/x**2)*exp(exp(x)/x)*exp(-exp(x) + exp(x)/(1 - 1/x) - \
            exp(x)/x - exp(x)/x**2)*exp(exp(x))

    n = Symbol('n', integer=True)
    e = (sqrt(n)*log(n)**2*exp(sqrt(log(n))*log(log(n))**2*exp(sqrt(log(log(n)))*log(log(log(n)))**3)))/n
    assert e.aseries(n) == exp(exp(sqrt(log(log(n)))*log(log(log(n)))**3)*sqrt(log(n))*log(log(n))**2)*log(n)**2/sqrt(n)


# TODO: add test cases involving special functions
def test_special():
    pass


def test_hierarchical():
    e = sin(1/x + exp(-x))
    assert e.aseries(x, n=3, hir=True) == -exp(-2*x)*sin(1/x)/2 + \
            exp(-x)*cos(1/x) + sin(1/x) + O(exp(-3*x), (x, oo))
    e = sin(x) * cos(exp(-x))
    assert e.aseries(x, hir=True) == exp(-4*x)*sin(x)/24 - \
            exp(-2*x)*sin(x)/2 + sin(x) + O(exp(-6*x), (x, oo))
    raises(PoleError, lambda: e.aseries(x))

    a, b = symbols('a b', integer=True)
    e = exp(1/x + exp(-x**2) * (exp(a*x) - exp(b*x))) - exp(1/x)
    assert e.aseries(x, n=3, hir=True) == (exp(2*a*x + 1/x)/2 + exp(2*b*x + 1/x)/2 - \
            exp(a*x + b*x + 1/x))*exp(-2*x**2) + (exp(a*x + 1/x) - exp(b*x + 1/x))*exp(-x**2) + O(exp(-3*x**2), (x, oo))

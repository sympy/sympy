from sympy import (
    Symbol, gamma, I, oo, nan, zoo, factorial, sqrt, Rational, log,
    polygamma, EulerGamma, pi, uppergamma, S, expand_func, loggamma, sin,
    cos, O, cancel, lowergamma, exp, erf, beta, exp_polar, harmonic, zeta,
    factorial)
from sympy.core.function import ArgumentIndexError
from sympy.utilities.randtest import (test_derivative_numerically as td,
                                      random_complex_number as randcplx,
                                      test_numerically as tn)
from sympy.utilities.pytest import raises

x = Symbol('x')
y = Symbol('y')
n = Symbol('n', integer=True)


def test_beta():
    x, y = Symbol('x'), Symbol('y')
    assert beta(x, y) == gamma(x)*gamma(y)/gamma(x + y)
    assert beta(x, y) == beta(y, x)  # Symmetric
    assert beta(x, y) == beta(x,y+1) + beta(x+1,y)

    assert diff(beta(x, y), x) == beta(x, y)*(digamma(x) - digamma(x+y))
    assert diff(beta(x, y), y) == beta(x, y)*(digamma(y) - digamma(x+y))

    
    
from sympy import (
    Symbol, gamma, I, oo, nan, zoo, factorial, sqrt, Rational, log, beta,
    polygamma, EulerGamma, pi, uppergamma, S, expand_func, loggamma, sin,
    cos, O, cancel, lowergamma, exp, erf, beta, exp_polar, harmonic, zeta,
    factorial, digamma)
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
    assert isinstance(beta(x, y), beta)
    assert expand_func(beta(x, y)) == gamma(x)*gamma(y)/gamma(x + y)
    assert expand_func(beta(x, y) - beta(y, x)) == 0 # Symmetric
    assert expand_func(beta(x, y)) == expand_func(beta(x,y + 1) + beta(x + 1,y)).simplify()
    assert diff(beta(x, y), x) == beta(x, y)*(digamma(x) - digamma(x + y))
    assert diff(beta(x, y), y) == beta(x, y)*(digamma(y) - digamma(x + y))
    assert beta(x, y) == (beta(x,y + 1) + beta(x + 1,y)).simplify()

    
    

from sympy import (Symbol, S, gamma, expand_func, beta, digamma, diff, conjugate, Integral,hyper, betainc, log)
from sympy.functions.special.gamma_functions import polygamma
from sympy.core.function import ArgumentIndexError
from sympy.testing.pytest import raises



def test_beta():
    x, y = Symbol('x'), Symbol('y')
    assert isinstance(beta(x, y), beta)
    assert expand_func(beta(x, y)) == gamma(x)*gamma(y)/gamma(x + y)
    assert expand_func(beta(x, y) - beta(y, x)) == 0  # Symmetric
    assert expand_func(beta(x, y)) == expand_func(beta(x, y + 1) + beta(x + 1, y)).simplify()
    assert diff(beta(x, y), x) == beta(x, y)*(polygamma(0, x) - polygamma(0, x + y))
    assert diff(beta(x, y), y) == beta(x, y)*(polygamma(0, y) - polygamma(0, x + y))

    assert conjugate(beta(x, y)) == beta(conjugate(x), conjugate(y))

    raises(ArgumentIndexError, lambda: beta(x, y).fdiff(3))

    assert beta(x, y).rewrite(gamma) == gamma(x)*gamma(y)/gamma(x + y)


def test_betainc():
    x, y, z, t = Symbol('x'), Symbol('y'), Symbol('z'), Symbol('t')
    assert isinstance(betainc(z, x, y), betainc)
    assert diff(betainc(z,x,y),z) == (1-z)**(y-1)*z**(x-1)
    assert betainc(1,x,y) == beta(x,y)
    assert betainc(0,x,y) == S.Zero
    assert betainc(z, x, y).rewrite(Integral) == Integral(t**(x - 1)*(1 - t)**(y - 1),(t, 0, z))
    assert betainc(z, x, y).rewrite(hyper) == gamma(x)*z**x*hyper([x, 1-y], [x+1], z)
    assert diff(betainc(z,x,y),y) == gamma(y)**2*(1 - z)**y*hyper([y, y, 1-x], [y+1, y+1], 1-z) - (log(1 - z)*betainc(z,x,y).func(1 - z, x, y)) + (digamma(y) - digamma(x + y))*beta(x, y)
    assert diff(betainc(z,x,y),z) == (1 - z)**(y - 1)*z**(x - 1)
    assert diff(betainc(z,x,y),x) == betainc(z, x, y)*log(z)-z**x*gamma(x)**2*hyper([x, x, 1-y], [x+1, x+1], z)
    raises(ArgumentIndexError, lambda: betainc(z,x, y).fdiff(4))
    raises(ArgumentIndexError, lambda: betainc(z,x, y).fdiff(0))
    assert conjugate(betainc(z, x, y)) == betainc(conjugate(z), conjugate(x), conjugate(y))

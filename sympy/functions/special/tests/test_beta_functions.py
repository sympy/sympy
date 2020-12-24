from sympy import (symbols, gamma, expand_func, beta, betainc,
        hyper, diff, conjugate, Integral, Dummy, I,
        regularized_betainc)
from sympy.functions.special.gamma_functions import polygamma
from sympy.core.function import ArgumentIndexError
from sympy.testing.pytest import raises


def test_beta():
    x, y = symbols('x y')
    t = Dummy('t')

    assert isinstance(beta(x, y), beta)

    assert expand_func(beta(x, y)) == gamma(x)*gamma(y)/gamma(x + y)
    assert expand_func(beta(x, y) - beta(y, x)) == 0  # Symmetric
    assert expand_func(beta(x, y)) == expand_func(beta(x, y + 1) + beta(x + 1, y)).simplify()

    assert diff(beta(x, y), x) == beta(x, y)*(polygamma(0, x) - polygamma(0, x + y))
    assert diff(beta(x, y), y) == beta(x, y)*(polygamma(0, y) - polygamma(0, x + y))

    assert conjugate(beta(x, y)) == beta(conjugate(x), conjugate(y))

    raises(ArgumentIndexError, lambda: beta(x, y).fdiff(3))

    assert beta(x, y).rewrite(gamma) == gamma(x)*gamma(y)/gamma(x + y)
    assert beta(x).rewrite(gamma) == gamma(x)**2/gamma(2*x)
    assert beta(x, y).rewrite(Integral).dummy_eq(Integral(t**(x - 1) * (1 - t)**(y - 1), (t, 0, 1)))

def test_betainc():
    a, b, x1, x2 = symbols('a b x1 x2')

    assert isinstance(betainc(a, b, x1, x2), betainc)
    assert isinstance(betainc(a, b, 0, x1), betainc)

    assert betainc(1, 2, 0, 2 + I).is_real == False
    assert conjugate(betainc(I, 2, 3 - I, 1 + 4*I)) == betainc(-I, 2, 3 + I, 1 - 4*I)

    assert betainc(a, b, 0, 1).rewrite(Integral).dummy_eq(beta(a, b).rewrite(Integral))
    assert betainc(1, 2, 0, x2).rewrite(hyper) == x2*hyper((1, -1), (2,), x2)

def test_regularized_betainc():
    a, b, x1, x2 = symbols('a b x1 x2')

    assert isinstance(regularized_betainc(a, b, x1, x2), regularized_betainc)
    assert isinstance(regularized_betainc(a, b, 0, x1), regularized_betainc)

    assert regularized_betainc(3, 5, 0, -1).is_real == True
    assert conjugate(regularized_betainc(3*I, 1, 2 + I, 1 + 2*I)) == regularized_betainc(-3*I, 1, 2 - I, 1 - 2*I)

    assert regularized_betainc(a, b, 0, 1).rewrite(Integral) == 1
    assert regularized_betainc(1, 2, x1, x2).rewrite(hyper) == 2*x2*hyper((1, -1), (2,), x2) - 2*x1*hyper((1, -1), (2,), x1)

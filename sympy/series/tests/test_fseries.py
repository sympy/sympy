from sympy import sin, exp, asin, Derivative, Symbol
from sympy.core import Function
from sympy.functions import airyai
from sympy.series.formal import simpleDE, DEtoRE
from sympy.abc import x


def test_simpleDE():
    f = Function('f')

    assert simpleDE(x, sin(x), f) == f(x) + Derivative(f(x), x, x)
    assert simpleDE(x, exp(x), f) == -f(x) + Derivative(f(x), x)
    assert simpleDE(x, sin(x)*exp(x), f) == 2*f(x) - 2*Derivative(f(x), x) + Derivative(f(x), x, x)
    assert simpleDE(x, airyai(x), f) == -x*f(x) + Derivative(f(x), x, x)
    assert simpleDE(x, asin(x), f) == x*Derivative(f(x), x) + (x**2 - 1)*Derivative(f(x), x, x)
    assert simpleDE(x, asin(x**5), f) == x**10*Derivative(f(x), x) + x*(x**10 - 1)*Derivative(f(x), x, x) + 4*Derivative(f(x), x)
    assert simpleDE(x, asin(x)**3, f) == x**4*Derivative(f(x), x, x, x, x) + 6*x**3*Derivative(f(x), x, x, x) + \
            7*x**2*Derivative(f(x), x, x) - 2*x**2*Derivative(f(x), x, x, x, x) + x*Derivative(f(x), x) - \
            6*x*Derivative(f(x), x, x, x) - 4*Derivative(f(x), x, x) + Derivative(f(x), x, x, x, x)


def test_DEtoRE():
    f = Function('f')
    r = Function('r')
    n = Symbol('n', integer=True)

    DE = f(x) + Derivative(f(x), x, x)
    assert DEtoRE(x, DE, r, n) == (n + 1)*(n + 2)*r(n + 2) + r(n)
    DE = 2*f(x) - 2*Derivative(f(x), x) + Derivative(f(x), x, x)
    assert DEtoRE(x, DE, r, n) == (-2*n - 2)*r(n + 1) + (n + 1)*(n + 2)*r(n + 2) + 2*r(n)

from sympy import (sin, cos, tan, sec, csc, cot, log, exp, atan, asin, acos,
                   Symbol, Integral, integrate, pi, Dummy, Derivative,
                   diff, I, sqrt, erf, Piecewise, Eq, Ne, symbols, Rational,
                   And, Heaviside, S, asinh, acosh, atanh, acoth, expand,
                   Function, jacobi, gegenbauer, chebyshevt, chebyshevu,
                   legendre, hermite, laguerre, assoc_laguerre, uppergamma, li,
                   Ei, Ci, Si, Chi, Shi, fresnels, fresnelc, polylog, erf, erfi,
                   sinh, cosh, elliptic_f, elliptic_e)
from sympy.integrals.manualintegrate import (manualintegrate, find_substitutions,
    _parts_rule, integral_steps, contains_dont_know, manual_subs)
from sympy.utilities.pytest import slow

x, y, z, u, n, a, b, c = symbols('x y z u n a b c')
f = Function('f')

def test_find_substitutions():
    assert find_substitutions((cot(x)**2 + 1)**2*csc(x)**2*cot(x)**2, x, u) == \
        [(cot(x), 1, -u**6 - 2*u**4 - u**2)]
    assert find_substitutions((sec(x)**2 + tan(x) * sec(x)) / (sec(x) + tan(x)),
                              x, u) == [(sec(x) + tan(x), 1, 1/u)]
    assert find_substitutions(x * exp(-x**2), x, u) == [(-x**2, -S.Half, exp(u))]


def test_manualintegrate_parts():
    # Make sure _parts_rule doesn't pick u = constant but can pick dv =
    # constant if necessary, e.g. for integrate(atan(x))
    assert _parts_rule(cos(x), x) == None
    assert _parts_rule(exp(x), x) == None
    assert _parts_rule(x**2, x) == None
    result = _parts_rule(atan(x), x)
    assert result[0] == atan(x) and result[1] == 1


def test_issue_6799():
    r, x, phi = map(Symbol, 'r x phi'.split())
    n = Symbol('n', integer=True, positive=True)

    integrand = (cos(n*(x-phi))*cos(n*x))
    limits = (x, -pi, pi)
    assert r * integrate(integrand, limits).trigsimp() / pi == r * cos(n * phi)
    assert not integrate(integrand, limits).has(Dummy)


def test_issue_3796():
    assert integrate(x * exp(x**4), x, risch=False) == -I*sqrt(pi)*erf(I*x**2)/4


def test_manual_true():
    assert integrate(exp(x) * sin(x), x, manual=True) == \
        (exp(x) * sin(x)) / 2 - (exp(x) * cos(x)) / 2
    assert integrate(sin(x) * cos(x), x, manual=True) in \
        [sin(x) ** 2 / 2, -cos(x)**2 / 2]


def test_issue_9462():
    assert not contains_dont_know(integral_steps(sin(2*x)*exp(x), x))


def test_manual_subs():
    x, y = symbols('x y')
    expr = log(x) + exp(x)
    # if log(x) is y, then exp(y) is x
    assert manual_subs(expr, log(x), y) == y + exp(exp(y))
    # if exp(x) is y, then log(y) need not be x
    assert manual_subs(expr, exp(x), y) == log(x) + y

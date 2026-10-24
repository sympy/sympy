"""Exact and numerical high-order derivatives using Taylor coefficients."""

from sympy import (Derivative, Rational, Subs, Symbol, cos, factorial,
                   sin, sqrt, symbols)


def test_point_derivative_exact():
    x = Symbol('x')
    f = sqrt(1 + x**2)
    for n in (1, 2, 4, 8, 12):
        d = Subs(Derivative(f, (x, n)), x, 2)
        assert (d.doit() - f.diff(x, n).subs(x, 2)).equals(0)

    d = Subs(Derivative(f, (x, 30)), x, 2)
    expected = Rational(118850792325494031296768942403977241,
                        11920928955078125)*sqrt(5)
    assert d.doit() == expected
    assert d.evalf(50) == expected.evalf(50)
    assert Derivative(f, (x, 30)).doit_numerically(2) == expected.evalf()


def test_point_derivative_parameters_and_polynomials():
    x, a, b = symbols('x a b')
    f = sqrt(1 + a*x**2) + b*x
    for n in (1, 2, 6, 12):
        d = Subs(Derivative(f, (x, n)), (x, a, b), (2, 3, 4))
        exact = f.diff(x, n).subs({x: 2, a: 3, b: 4})
        assert (d.doit() - exact).equals(0)

    assert Subs(Derivative(x**9, (x, 9)), x, 0).doit() == factorial(9)
    assert Subs(Derivative(x**9, (x, 10)), x, 0).doit() == 0
    assert Subs(Derivative(1/(x + 1), (x, 20)), x, 0).doit() == factorial(20)


def test_point_derivative_optical():
    x, c = symbols('x c')
    f = x/c*sqrt(7 + 9/((c/x)**2 - 5) + 5/(6 - (c/x)**2) + 8/(c/x))
    d = Subs(Derivative(f, (x, 12)), (x, c), (2, 10))
    assert abs(d.doit().evalf(20) + 2600.6802537683420886) < 1e-10
    d30 = Subs(Derivative(f, (x, 30)), (x, c), (2, 10))
    assert abs(d30.doit().evalf(15) + 4.339172027035365e21) < 1e8


def test_point_derivative_unsupported_fallbacks():
    x, y = symbols('x y')
    assert Subs(Derivative(sin(x), (x, 3)), x, 2).doit() == -cos(2)
    assert Subs(Derivative(x*y, x, y), (x, y), (2, 3)).doit() == 1
    singular = Subs(Derivative(sqrt(x), (x, 2)), x, 0)
    assert singular.doit() == sqrt(x).diff(x, 2).subs(x, 0)
    f = sqrt(1 + y*x**2)
    d = Subs(Derivative(f, (x, 2)), x, 2)
    assert (d.doit() - f.diff(x, 2).subs(x, 2)).equals(0)
    d = Subs(Derivative(x**3, (x, 2)), (x, y), (2, x + 1))
    assert d.doit() == 12


def test_point_derivative_symbolic_order():
    from sympy import S, rf
    from sympy.core._point_derivative import _subs_derivative_at_point

    x, y, n = symbols('x y n')
    f = 1/sqrt(x + y) - sqrt(1/x)

    d = Subs(Derivative(f, (x, n)), (n, y, x), (33, 2, 10))

    expected = rf(S.Half, 33)*(
        S(10)**(-S(67)/2) - S(12)**(-S(67)/2))

    assert d.doit() == expected
    assert _subs_derivative_at_point(
        d.expr, d.variables, d.point) == expected

    d = Subs(Derivative(f, (x, n + 1)), (n, y, x), (32, 2, 10))
    assert _subs_derivative_at_point(
        d.expr, d.variables, d.point) == expected

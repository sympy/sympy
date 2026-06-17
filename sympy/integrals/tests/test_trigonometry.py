from __future__ import annotations
from sympy.core import Ne, Rational, Symbol
from sympy.functions import sin, cos, tan, csc, sec, cot, log, Piecewise
from sympy.integrals.trigonometry import trigintegrate

x = Symbol('x')


def test_trigintegrate_odd():
    assert trigintegrate(Rational(1), x) == x
    assert trigintegrate(x, x) is None
    assert trigintegrate(x**2, x) is None

    assert trigintegrate(sin(x), x) == -cos(x)
    assert trigintegrate(cos(x), x) == sin(x)

    assert trigintegrate(sin(3*x), x) == -cos(3*x)/3
    assert trigintegrate(cos(3*x), x) == sin(3*x)/3

    y = Symbol('y')
    assert trigintegrate(sin(y*x), x) == Piecewise(
        (-cos(y*x)/y, Ne(y, 0)), (0, True))
    assert trigintegrate(cos(y*x), x) == Piecewise(
        (sin(y*x)/y, Ne(y, 0)), (x, True))
    assert trigintegrate(sin(y*x)**2, x) == Piecewise(
        ((x*y/2 - sin(x*y)*cos(x*y)/2)/y, Ne(y, 0)), (0, True))
    assert trigintegrate(sin(y*x)*cos(y*x), x) == Piecewise(
        (sin(x*y)**2/(2*y), Ne(y, 0)), (0, True))
    assert trigintegrate(cos(y*x)**2, x) == Piecewise(
        ((x*y/2 + sin(x*y)*cos(x*y)/2)/y, Ne(y, 0)), (x, True))

    y = Symbol('y', positive=True)
    # TODO: remove conds='none' below. For this to work we would have to rule
    #       out (e.g. by trying solve) the condition y = 0, incompatible with
    #       y.is_positive being True.
    assert trigintegrate(sin(y*x), x, conds='none') == -cos(y*x)/y
    assert trigintegrate(cos(y*x), x, conds='none') == sin(y*x)/y

    assert trigintegrate(sin(x)*cos(x), x) == sin(x)**2/2
    assert trigintegrate(sin(x)*cos(x)**2, x) == -cos(x)**3/3
    assert trigintegrate(sin(x)**2*cos(x), x) == sin(x)**3/3

    # check if it selects right function to substitute,
    # so the result is kept simple
    assert trigintegrate(sin(x)**7 * cos(x), x) == sin(x)**8/8
    assert trigintegrate(sin(x) * cos(x)**7, x) == -cos(x)**8/8

    assert trigintegrate(sin(x)**7 * cos(x)**3, x) == \
        -sin(x)**10/10 + sin(x)**8/8
    assert trigintegrate(sin(x)**3 * cos(x)**7, x) == \
        cos(x)**10/10 - cos(x)**8/8

    # both n, m are odd and -ve, and not necessarily equal
    assert trigintegrate(sin(x)**-1*cos(x)**-1, x) == \
        -log(sin(x)**2 - 1)/2 + log(sin(x))


def test_trigintegrate_even():
    assert trigintegrate(sin(x)**2, x) == x/2 - cos(x)*sin(x)/2
    assert trigintegrate(cos(x)**2, x) == x/2 + cos(x)*sin(x)/2

    assert trigintegrate(sin(3*x)**2, x) == x/2 - cos(3*x)*sin(3*x)/6
    assert trigintegrate(cos(3*x)**2, x) == x/2 + cos(3*x)*sin(3*x)/6
    assert trigintegrate(sin(x)**2 * cos(x)**2, x) == \
        x/8 - sin(2*x)*cos(2*x)/16

    assert trigintegrate(sin(x)**4 * cos(x)**2, x) == \
        x/16 - sin(x) *cos(x)/16 - sin(x)**3*cos(x)/24 + \
        sin(x)**5*cos(x)/6

    assert trigintegrate(sin(x)**2 * cos(x)**4, x) == \
        x/16 + cos(x) *sin(x)/16 + cos(x)**3*sin(x)/24 - \
        cos(x)**5*sin(x)/6

    assert trigintegrate(sin(x)**(-4), x) == -2*cos(x)/(3*sin(x)) \
        - cos(x)/(3*sin(x)**3)

    assert trigintegrate(cos(x)**(-6), x) == sin(x)/(5*cos(x)**5) \
        + 4*sin(x)/(15*cos(x)**3) + 8*sin(x)/(15*cos(x))


def test_trigintegrate_mixed():
    assert trigintegrate(sin(x)*sec(x), x) == -log(cos(x))
    assert trigintegrate(sin(x)*csc(x), x) == x
    assert trigintegrate(sin(x)*cot(x), x) == sin(x)

    assert trigintegrate(cos(x)*sec(x), x) == x
    assert trigintegrate(cos(x)*csc(x), x) == log(sin(x))
    assert trigintegrate(cos(x)*tan(x), x) == -cos(x)
    assert trigintegrate(cos(x)*cot(x), x) == log(cos(x) - 1)/2 \
        - log(cos(x) + 1)/2 + cos(x)
    assert trigintegrate(cot(x)*cos(x)**2, x) == log(sin(x)) - sin(x)**2/2


def test_trigintegrate_symbolic():
    n = Symbol('n', integer=True)
    assert trigintegrate(cos(x)**n, x) is None
    assert trigintegrate(sin(x)**n, x) is None
    assert trigintegrate(cot(x)**n, x) is None


def test_trigintegrate_rational_exponents():
    """Regression test for issue #29882:
    trigintegrate should handle Rational exponents via Incomplete Beta."""
    from sympy.functions.special.beta_functions import betainc

    # Basic case from the issue report: sin(x)**(1/2) * cos(x)**(3/2)
    result = trigintegrate(sin(x)**Rational(1, 2) * cos(x)**Rational(3, 2), x)
    assert result is not None
    expected = betainc(Rational(3, 4), Rational(5, 4), 0, sin(x)**2) / 2
    assert result == expected

    # Both exponents rational: sin(x)**(1/3) * cos(x)**(2/3)
    result2 = trigintegrate(sin(x)**Rational(1, 3) * cos(x)**Rational(2, 3), x)
    assert result2 is not None
    expected2 = betainc(Rational(2, 3), Rational(5, 6), 0, sin(x)**2) / 2
    assert result2 == expected2

    # One integer, one rational: sin(x)**2 * cos(x)**(1/2)
    result3 = trigintegrate(sin(x)**2 * cos(x)**Rational(1, 2), x)
    assert result3 is not None
    expected3 = betainc(Rational(3, 2), Rational(3, 4), 0, sin(x)**2) / 2
    assert result3 == expected3

    # Negative rational exponent (still > -1): sin(x)**(-1/2) * cos(x)**(3/2)
    result4 = trigintegrate(sin(x)**Rational(-1, 2) * cos(x)**Rational(3, 2), x)
    assert result4 is not None
    expected4 = betainc(Rational(1, 4), Rational(5, 4), 0, sin(x)**2) / 2
    assert result4 == expected4

    # Integer exponents should NOT be rerouted (still handled by original path)
    assert trigintegrate(sin(x)**2 * cos(x)**3, x) == \
        sin(x)**3/3 - sin(x)**5/5


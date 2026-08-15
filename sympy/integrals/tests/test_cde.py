"""Most of these tests come from the examples in Bronstein's book."""
from sympy.core.numbers import oo
from sympy.core.singleton import S
from sympy.polys.polytools import Poly, cancel
from sympy.simplify.simplify import simplify
from sympy.integrals.risch import (DifferentialExtension, derivation,
    NonElementaryIntegralException)
from sympy.integrals.cde import (coupled_DE_system, coupled_DE_cancel_prim,
    coupled_DE_cancel_exp, coupled_DE_cancel_tan)

from sympy.testing.pytest import raises
from sympy.abc import x, t


def _check(q1, q2, b1e, b2e, c1e, c2e, DE):
    """(q1, q2) solves Dq1 + b1*q1 - b2*q2 == c1, Dq2 + b2*q1 + b1*q2 == c2."""
    q1e, q2e = q1.as_expr(), q2.as_expr()
    assert simplify(derivation(q1, DE).as_expr() +
        b1e*q1e - b2e*q2e - c1e) == 0
    assert simplify(derivation(q2, DE).as_expr() +
        b2e*q1e + b1e*q2e - c2e) == 0


def test_coupled_DE_system():
    DE = DifferentialExtension(extension={'D': [Poly(1, x)]})
    one = Poly(1, x)
    # From the integration of tan(x**2) (Section 5.10):
    # CoupledDESystem(0, 1, 2/x, 0) has no solution, so tan(x**2) has
    # no elementary integral
    raises(NonElementaryIntegralException, lambda: coupled_DE_system(
        (Poly(0, x), one), (one, one), (Poly(2, x), Poly(x, x)),
        (Poly(0, x), one), DE))
    # From Example 6.6.2: CoupledDESystem(1, -1, -2*x, -2*(x + 1))
    # == (0, -2*x)
    (y1a, y1d), (y2a, y2d) = coupled_DE_system((one, one), (Poly(-1, x), one),
        (Poly(-2*x, x), one), (Poly(-2*(x + 1), x), one), DE)
    assert cancel(y1a.as_expr()/y1d.as_expr()) == 0
    assert cancel(y2a.as_expr()/y2d.as_expr()) == -2*x
    # The inner system of Example 8.4.1:
    # CoupledDESystem(0, 4*x - 2, 2 - 8*x**2, 4 - 4*x) == (-1, 2*x + 1)
    (y1a, y1d), (y2a, y2d) = coupled_DE_system((Poly(0, x), one),
        (Poly(4*x - 2, x), one), (Poly(2 - 8*x**2, x), one),
        (Poly(4 - 4*x, x), one), DE)
    assert cancel(y1a.as_expr()/y1d.as_expr()) == -1
    assert cancel(y2a.as_expr()/y2d.as_expr()) == 2*x + 1


def test_coupled_DE_cancel_prim():
    # t = log(x)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)]})
    # b == b1 + b2*sqrt(-1) == Dz/z for z == x + sqrt(-1); the system
    # was built from the solution (t, 1), and the returned solution
    # differs from it by the homogeneous solution -x/(x**2 + 1) +
    # sqrt(-1)/(x**2 + 1) == -1/(x + sqrt(-1))
    b1e = x/(x**2 + 1)
    b2e = -1/(x**2 + 1)
    c1e = cancel(1/x + b1e*t - b2e)
    c2e = cancel(b2e*t + b1e)
    q1, q2 = coupled_DE_cancel_prim(Poly(b1e, t), Poly(b2e, t),
        Poly(c1e, t), Poly(c2e, t), DE, 1)
    _check(q1, q2, b1e, b2e, c1e, c2e, DE)
    assert q1.degree(t) == 1
    # b == sqrt(-1) is not a logarithmic derivative of k(sqrt(-1))*,
    # so this exercises the coefficient-descent loop; the solution
    # (t, 0) is unique
    q1, q2 = coupled_DE_cancel_prim(Poly(0, t), Poly(1, t),
        Poly(1/x, t), Poly(t, t), DE, 1)
    assert (q1.as_expr(), q2.as_expr()) == (t, S.Zero)
    # no solution: Dq1 == 1/x + t needs q1 == t*(t + ...)/2 of degree 2
    raises(NonElementaryIntegralException, lambda: coupled_DE_cancel_prim(
        Poly(0, t), Poly(1, t), Poly(1/x + t, t), Poly(t, t), DE, 1))


def test_coupled_DE_cancel_exp():
    # t = exp(x)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    # b == Dz/z + m*Dt/t for z == x + sqrt(-1) and m == 1; the system
    # was built from the solution (t, 0), which is unique here
    b1e = cancel(x/(x**2 + 1) + 1)
    b2e = -S(1)/(x**2 + 1)
    c1e = cancel(t + b1e*t)
    c2e = cancel(b2e*t)
    q1, q2 = coupled_DE_cancel_exp(Poly(b1e, t), Poly(b2e, t),
        Poly(c1e, t), Poly(c2e, t), DE, 1)
    assert (q1.as_expr(), q2.as_expr()) == (t, S.Zero)
    _check(q1, q2, b1e, b2e, c1e, c2e, DE)
    # b == sqrt(-1) is not of the form Dz/z + m*Dt/t, so this
    # exercises the coefficient-descent loop; the solution (t, t) is
    # unique
    q1, q2 = coupled_DE_cancel_exp(Poly(0, t), Poly(1, t),
        Poly(0, t), Poly(2*t, t), DE, 1)
    assert (q1.as_expr(), q2.as_expr()) == (t, t)
    # no solution: the t-free parts require Ds1 - s2 == 2/x and
    # Ds2 + s1 == 0 in k, the system whose unsolvability proves that
    # tan(x**2) has no elementary integral (Section 5.10)
    raises(NonElementaryIntegralException, lambda: coupled_DE_cancel_exp(
        Poly(0, t), Poly(1, t), Poly(2/x, t), Poly(2*t, t), DE, 1))


def test_coupled_DE_cancel_tan():
    # Example 8.4.1: t = tan(x), from the integration of (8.12); the
    # system (8.15) with b0 == 0, b2 == 4*x, n == 2 has the solution
    # (t - 1, 2*x)
    DE = DifferentialExtension(extension={'D': [Poly(1, x),
        Poly(t**2 + 1, t)]})
    b2e = 4*x
    c1e = -t**2 + 2*t - 8*x**2 + 1
    c2e = 2*(1 - 2*x)
    q1, q2 = coupled_DE_cancel_tan(Poly(0, t), Poly(b2e, t),
        Poly(c1e, t), Poly(c2e, t), DE, 2)
    assert (q1.as_expr(), q2.as_expr()) == (t - 1, 2*x)
    # the diagonal entries of the system are b0 - n*eta*t (the matrix
    # in the book's CoupledDECancelTan specification misprints the
    # bottom-right entry as b0 + n*eta*t)
    _check(q1, q2, -2*t, b2e, c1e, c2e, DE)
    # n == 0: the system for constant-degree solutions is a coupled
    # system over k; (x, 1) is the unique solution
    q1, q2 = coupled_DE_cancel_tan(Poly(0, t), Poly(1, t),
        Poly(0, t), Poly(x, t), DE, 0)
    assert (q1.as_expr(), q2.as_expr()) == (x, S.One)
    _check(q1, q2, S.Zero, S.One, S.Zero, x, DE)
    # n == 0 with deg(c) > 0 has no solution
    raises(NonElementaryIntegralException, lambda: coupled_DE_cancel_tan(
        Poly(0, t), Poly(1, t), Poly(t, t), Poly(x, t), DE, 0))
    # the recursion descends from n to 0, so n == oo is not allowed
    raises(ValueError, lambda: coupled_DE_cancel_tan(
        Poly(0, t), Poly(1, t), Poly(0, t), Poly(x, t), DE, oo))

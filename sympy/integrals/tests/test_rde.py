"""Most of these tests come from the examples in Bronstein's book."""
from __future__ import annotations
from sympy.core.numbers import (I, Rational, oo)
from sympy.core.symbol import symbols
from sympy.polys.polytools import Poly, cancel
from sympy.integrals.risch import (DifferentialExtension, derivation,
    NonElementaryIntegralException)
from sympy.integrals.rde import (order_at, order_at_oo, weak_normalizer,
    normal_denom, special_denom, _special_denom_cancel_bound, bound_degree,
    spde, solve_poly_rde, no_cancel_equal, cancel_primitive, cancel_exp,
    rischDE)

from sympy.testing.pytest import raises
from sympy.abc import x, t, z, n

t0, t1, t2, k = symbols('t:3 k')


def test_order_at():
    a = Poly(t**4, t)
    b = Poly((t**2 + 1)**3*t, t)
    c = Poly((t**2 + 1)**6*t, t)
    d = Poly((t**2 + 1)**10*t**10, t)
    e = Poly((t**2 + 1)**100*t**37, t)
    p1 = Poly(t, t)
    p2 = Poly(1 + t**2, t)
    assert order_at(a, p1, t) == 4
    assert order_at(b, p1, t) == 1
    assert order_at(c, p1, t) == 1
    assert order_at(d, p1, t) == 10
    assert order_at(e, p1, t) == 37
    assert order_at(a, p2, t) == 0
    assert order_at(b, p2, t) == 3
    assert order_at(c, p2, t) == 6
    assert order_at(d, p1, t) == 10
    assert order_at(e, p2, t) == 100
    assert order_at(Poly(0, t), Poly(t, t), t) is oo
    assert order_at_oo(Poly(t**2 - 1, t), Poly(t + 1), t) == \
        order_at_oo(Poly(t - 1, t), Poly(1, t), t) == -1
    assert order_at_oo(Poly(0, t), Poly(1, t), t) is oo

def test_weak_normalizer():
    a = Poly((1 + x)*t**5 + 4*t**4 + (-1 - 3*x)*t**3 - 4*t**2 + (-2 + 2*x)*t, t)
    d = Poly(t**4 - 3*t**2 + 2, t)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    r = weak_normalizer(a, d, DE, z)
    assert r == (Poly(t**5 - t**4 - 4*t**3 + 4*t**2 + 4*t - 4, t, domain='ZZ[x]'),
        (Poly((1 + x)*t**2 + x*t, t, domain='ZZ[x]'),
         Poly(t + 1, t, domain='ZZ[x]')))
    assert weak_normalizer(r[1][0], r[1][1], DE) == (Poly(1, t), r[1])
    r = weak_normalizer(Poly(1 + t**2), Poly(t**2 - 1, t), DE, z)
    assert r == (Poly(t**4 - 2*t**2 + 1, t), (Poly(-3*t**2 + 1, t), Poly(t**2 - 1, t)))
    assert weak_normalizer(r[1][0], r[1][1], DE, z) == (Poly(1, t), r[1])
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1 + t**2)]})
    r = weak_normalizer(Poly(1 + t**2), Poly(t, t), DE, z)
    assert r == (Poly(t, t), (Poly(0, t), Poly(1, t)))
    assert weak_normalizer(r[1][0], r[1][1], DE, z) == (Poly(1, t), r[1])


def test_normal_denom():
    DE = DifferentialExtension(extension={'D': [Poly(1, x)]})
    raises(NonElementaryIntegralException, lambda: normal_denom(Poly(1, x), Poly(1, x),
    Poly(1, x), Poly(x, x), DE))
    fa, fd = Poly(t**2 + 1, t), Poly(1, t)
    ga, gd = Poly(1, t), Poly(t**2, t)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t**2 + 1, t)]})
    assert normal_denom(fa, fd, ga, gd, DE) == \
        (Poly(t, t), (Poly(t**3 - t**2 + t - 1, t), Poly(1, t)), (Poly(1, t),
        Poly(1, t)), Poly(t, t))


def test_special_denom():
    # TODO: add more tests here
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    assert special_denom(Poly(1, t), Poly(t**2, t), Poly(1, t), Poly(t**2 - 1, t),
    Poly(t, t), DE) == \
        (Poly(1, t), Poly(t**2 - 1, t), Poly(t**2 - 1, t), Poly(t, t))
#    assert special_denom(Poly(1, t), Poly(2*x, t), Poly((1 + 2*x)*t, t), DE) == 1

    # issue 3940
    # Note, this isn't a very good test, because the denominator is just 1,
    # but at least it tests the exp cancellation case
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(-2*x*t0, t0),
        Poly(I*k*t1, t1)]})
    DE.decrement_level()
    assert special_denom(Poly(1, t0), Poly(I*k, t0), Poly(1, t0), Poly(t0, t0),
    Poly(1, t0), DE) == \
        (Poly(1, t0, domain='ZZ'), Poly(I*k, t0, domain='ZZ_I[k,x]'),
                Poly(t0, t0, domain='ZZ'), Poly(1, t0, domain='ZZ'))


    assert special_denom(Poly(1, t), Poly(t**2, t), Poly(1, t), Poly(t**2 - 1, t),
    Poly(t, t), DE, case='tan') == \
           (Poly(1, t, t0, domain='ZZ'), Poly(t**2, t0, t, domain='ZZ[x]'),
            Poly(t, t, t0, domain='ZZ'), Poly(1, t0, domain='ZZ'))

    raises(ValueError, lambda: special_denom(Poly(1, t), Poly(t**2, t), Poly(1, t), Poly(t**2 - 1, t),
    Poly(t, t), DE, case='unrecognized_case'))

    # Example 6.2.1
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    a = Poly(t**2 + 2*x*t + x**2, t)
    b = Poly((1 + 1/x**2)*t**2 + (2*x - 1 + 2/x)*t + x**2, t)
    c = Poly(t/x**2 - 1 + 2/x, t)
    assert special_denom(a, b, Poly(1, t), c, Poly(1, t), DE) == \
        (a, Poly(t**2/x**2 + (2/x - 1)*t, t), Poly(t**2/x**2 + (2/x - 1)*t, t), Poly(t, t))

    # Example 6.2.2, adding a hypertangent case where the special denominator
    # is not required in the reduction of the system, hence returning h = 1
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1 + t**2, t)]})
    a = Poly(t, t, domain='ZZ')
    ba = Poly((t - 1)*(t**2 + 1), t, domain='ZZ')
    bd = Poly(1, t, domain='ZZ')
    ca = Poly(1, t, domain='ZZ')
    cd = Poly(1, t, domain='ZZ')
    expected_a = Poly(t, t, domain='ZZ')
    expected_b = Poly(t**3 - t**2 + t - 1, t, domain='ZZ')
    expected_c = Poly(1, t, domain='ZZ')
    expected_h = Poly(1, t, domain='ZZ')
    assert special_denom(a, ba, bd, ca, cd, DE) == \
        (expected_a, expected_b, expected_c, expected_h)

    # A hypertangent case with a nontrivial special denominator:
    # t*Dq + t*q == (t - 2*t**2)/(t**2 + 1) has the solution
    # q == 1/(t**2 + 1), so with h == t**2 + 1, r == q*h == 1 must satisfy
    # A*Dr + B*r == C.
    A, B, C, h = special_denom(Poly(t, t), Poly(t, t), Poly(1, t),
        Poly(-2*t**2 + t, t), Poly(t**2 + 1, t), DE)
    assert (A, B, C, h) == (Poly(t, t), Poly(-2*t**2 + t, t),
        Poly(-2*t**2 + t, t), Poly(t**2 + 1, t))
    # Verify the contract explicitly (r == q*h == 1, so Dr == 0)
    assert C == B  # A*D(1) + B*1 == C

    # The possible-cancellation sharpening itself: a == 1, b == 2*t over
    # Dt == t**2 + 1 admits the solution q == 1/(t**2 + 1) of
    # Dq + 2*t*q == 0, so the bound must drop from 0 to -1 (the earlier
    # hypertangent cases enter the helper with n already at -1, which
    # cannot detect a broken adjustment).
    assert _special_denom_cancel_bound(Poly(1, t), Poly(2*t, t),
        Poly(1, t), 0, DE, 'tan') == -1

def test_bound_degree_fail():
    # Primitive
    DE = DifferentialExtension(extension={'D': [Poly(1, x),
        Poly(t0/x**2, t0), Poly(1/x, t)]})
    assert bound_degree(Poly(t**2, t), Poly(-(1/x**2*t**2 + 1/x), t),
        Poly((2*x - 1)*t**4 + (t0 + x)/x*t**3 - (t0 + 4*x**2)/2*x*t**2 + x*t,
        t), DE) == 3


def test_bound_degree():
    # Base
    DE = DifferentialExtension(extension={'D': [Poly(1, x)]})
    assert bound_degree(Poly(1, x), Poly(-2*x, x), Poly(1, x), DE) == 0

    # Primitive (see above test_bound_degree_fail)
    # TODO: Add test for when the degree bound becomes larger after limited_integrate
    # TODO: Add test for db == da - 1 case

    # Exp
    # TODO: Add tests
    # TODO: Add test for when the degree becomes larger after parametric_log_deriv()

    # Nonlinear
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t**2 + 1, t)]})
    assert bound_degree(Poly(t, t), Poly((t - 1)*(t**2 + 1), t), Poly(1, t), DE) == 0


def test_spde():
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t**2 + 1, t)]})
    raises(NonElementaryIntegralException, lambda: spde(Poly(t, t), Poly((t - 1)*(t**2 + 1), t), Poly(1, t), 0, DE))
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    assert spde(Poly(t**2 + x*t*2 + x**2, t), Poly(t**2/x**2 + (2/x - 1)*t, t),
        Poly(t**2/x**2 + (2/x - 1)*t, t), 0, DE) == \
        (Poly(0, t), Poly(0, t), 0, Poly(0, t), Poly(1, t, domain='ZZ(x)'))
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t0/x**2, t0), Poly(1/x, t)]})
    assert spde(Poly(t**2, t), Poly(-t**2/x**2 - 1/x, t),
    Poly((2*x - 1)*t**4 + (t0 + x)/x*t**3 - (t0 + 4*x**2)/(2*x)*t**2 + x*t, t), 3, DE) == \
        (Poly(0, t), Poly(0, t), 0, Poly(0, t),
        Poly(t0*t**2/2 + x**2*t**2 - x**2*t, t, domain='ZZ(x,t0)'))
    DE = DifferentialExtension(extension={'D': [Poly(1, x)]})
    assert spde(Poly(x**2 + x + 1, x), Poly(-2*x - 1, x), Poly(x**5/2 +
    3*x**4/4 + x**3 - x**2 + 1, x), 4, DE) == \
        (Poly(0, x, domain='QQ'), Poly(x/2 - Rational(1, 4), x), 2, Poly(x**2 + x + 1, x), Poly(x*Rational(5, 4), x))
    assert spde(Poly(x**2 + x + 1, x), Poly(-2*x - 1, x), Poly(x**5/2 +
    3*x**4/4 + x**3 - x**2 + 1, x), n, DE) == \
        (Poly(0, x, domain='QQ'), Poly(x/2 - Rational(1, 4), x), -2 + n, Poly(x**2 + x + 1, x), Poly(x*Rational(5, 4), x))
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1, t)]})
    raises(NonElementaryIntegralException, lambda: spde(Poly((t - 1)*(t**2 + 1)**2, t), Poly((t - 1)*(t**2 + 1), t), Poly(1, t), 0, DE))
    DE = DifferentialExtension(extension={'D': [Poly(1, x)]})
    assert spde(Poly(x**2 - x, x), Poly(1, x), Poly(9*x**4 - 10*x**3 + 2*x**2, x), 4, DE) == \
        (Poly(0, x, domain='ZZ'), Poly(0, x), 0, Poly(0, x), Poly(3*x**3 - 2*x**2, x, domain='QQ'))
    assert spde(Poly(x**2 - x, x), Poly(x**2 - 5*x + 3, x), Poly(x**7 - x**6 - 2*x**4 + 3*x**3 - x**2, x), 5, DE) == \
        (Poly(1, x, domain='QQ'), Poly(x + 1, x, domain='QQ'), 1, Poly(x**4 - x**3, x), Poly(x**3 - x**2, x, domain='QQ'))

def test_solve_poly_rde_no_cancel():
    # deg(b) large
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1 + t**2, t)]})
    assert solve_poly_rde(Poly(t**2 + 1, t), Poly(t**3 + (x + 1)*t**2 + t + x + 2, t),
    oo, DE) == Poly(t + x, t)
    # deg(b) small
    DE = DifferentialExtension(extension={'D': [Poly(1, x)]})
    assert solve_poly_rde(Poly(0, x), Poly(x/2 - Rational(1, 4), x), oo, DE) == \
        Poly(x**2/4 - x/4, x)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t**2 + 1, t)]})
    assert solve_poly_rde(Poly(2, t), Poly(t**2 + 2*t + 3, t), 1, DE) == \
        Poly(t + 1, t, x)
    # deg(b) == deg(D) - 1
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t**2 + 1, t)]})
    assert no_cancel_equal(Poly(1 - t, t),
    Poly(t**3 + t**2 - 2*x*t - 2*x, t), oo, DE) == \
        (Poly(t**2, t), 1, Poly((-2 - 2*x)*t - 2*x, t))
    # The m == 0 branch: Dq + t*q == t has the solution q == 1 (this used
    # to return an Expr instead of a Poly)
    assert no_cancel_equal(Poly(t, t), Poly(t, t), 5, DE) == Poly(1, t)


def test_solve_poly_rde_cancel():
    # exp
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    assert cancel_exp(Poly(2*x, t), Poly(2*x, t), 0, DE) == \
        Poly(1, t)
    assert cancel_exp(Poly(2*x, t), Poly((1 + 2*x)*t, t), 1, DE) == \
        Poly(t, t)
    # TODO: Add more exp tests, including tests that require is_deriv_in_field()

    # primitive
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)]})

    # If the DecrementLevel context manager is working correctly, this shouldn't
    # cause any problems with the further tests.
    raises(NonElementaryIntegralException, lambda: cancel_primitive(Poly(1, t), Poly(t, t), oo, DE))

    # b == 1/(2*x) has radical index 2 (2*b == Dx/x, so
    # is_log_deriv_k_t_radical_in_field() returns (2, x)); the polynomial
    # loop must keep the caller's degree bound n == 3 rather than
    # clobbering it with the radical index (which would wrongly raise
    # NonElementaryIntegralException for the solution q == t**3).
    b = Poly(1/(2*x), t, field=True)
    q0 = Poly(t**3, t)
    c = (derivation(q0, DE) + b*q0).to_field()
    assert cancel_primitive(b, c, 3, DE).as_expr() == t**3

    assert cancel_primitive(Poly(1, t), Poly(t + 1/x, t), 2, DE) == \
        Poly(t, t)
    assert cancel_primitive(Poly(4*x, t), Poly(4*x*t**2 + 2*t/x, t), 3, DE) == \
        Poly(t**2, t)

    # The degree bound n must not be clobbered by the index returned from
    # is_log_deriv_k_t_radical_in_field(): b == 1/(2*x) has index 2
    # (2*b == D(x)/x), which used to replace n == 3 and wrongly reject
    # the degree 3 solution q == t**3 of Dq + b*q == c.
    b = Poly(Rational(1, 2)/x, t)
    q = Poly(t**3, t)
    c = derivation(q, DE) + b*q
    assert cancel_primitive(b, c, 3, DE) == q

    # TODO: Add more primitive tests, including tests that require is_deriv_in_field()


def test_rischDE():
    # TODO: Add more tests for rischDE, including ones from the text
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    DE.decrement_level()
    assert rischDE(Poly(-2*x, x), Poly(1, x), Poly(1 - 2*x - 2*x**2, x),
    Poly(1, x), DE) == \
        (Poly(x + 1, x), Poly(1, x))

    # See issue 28407
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    fa = Poly(-t + 1, t)
    fd = Poly(1, t)
    ga = Poly(1, t)
    gd = Poly(1, t)
    assert rischDE(fa, fd, ga, gd, DE) == (Poly(-1, t), Poly(t, t))

    # A derivation with an 'other_linear' monomial raises a graceful
    # NotImplementedError (it used to be an uncaught ValueError)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t + 1, t)]})
    raises(NotImplementedError, lambda: rischDE(Poly(1, t), Poly(1, t),
        Poly(t, t), Poly(1, t), DE))

    # Dy + y/x == 1.  f == 1/x is not weakly normalized; the weak normalizer
    # q == x must be applied to g and divided back out of the solution.
    # rischDE() used to discard it and return y == x, which solves
    # Dy + y/x == 2, not 1.  The correct answer here is y == x/2 (plus any
    # multiple of the homogeneous solution 1/x).
    DE = DifferentialExtension(extension={'D': [Poly(1, x)]})
    ya, yd = rischDE(Poly(1, x), Poly(x, x), Poly(1, x), Poly(1, x), DE)
    assert (ya, yd) == (Poly(x**2/2, x), Poly(x, x))
    # Explicitly verify Dy + f*y == g for the returned solution
    y = ya.as_expr()/yd.as_expr()
    assert cancel(y.diff(x) + y/x - 1) == 0

"""Most of these tests come from the examples in Bronstein's book."""
from __future__ import annotations
from sympy.core.numbers import (I, Rational, oo)
from sympy.core.singleton import S
from sympy.core.symbol import symbols
from sympy.polys.polytools import Poly, cancel
from sympy.integrals.risch import (DifferentialExtension, derivation,
    NonElementaryIntegralException)
from sympy.integrals.rde import (order_at, order_at_oo, weak_normalizer,
    normal_denom, special_denom, _special_denom_cancel_bound, bound_degree,
    spde, solve_poly_rde, no_cancel_equal, cancel_primitive, cancel_exp,
    cancel_tan, rischDE)

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


def test_bound_degree_rational_z():
    # Primitive case with deg(a) == deg(b) and alpha == 1/(x*(x + 1)) ==
    # Dz/z for z == x/(x + 1) -- a proper ratio in k*, so derivation(z, DE)
    # inside bound_degree() needs basic=True (z is not polynomial in the
    # tower generators; this used to crash with SympifyError).
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)],
        'exts': ['log'], 'extargs': [x]})
    assert bound_degree(Poly(t, t, field=True),
        Poly(-t/(x*(x + 1)), t, field=True), Poly(1, t, field=True), DE) == 0
    # The same alpha end-to-end: Dy - y/(x*(x + 1)) == (x + 1 - t)/(x*(x + 1))
    # has the solution y == t (plus c*x/(x + 1) for any constant c, so check
    # the defining equation rather than the exact result).
    ya, yd = rischDE(Poly(-1, t, field=True), Poly(x*(x + 1), t, field=True),
        Poly(x + 1 - t, t, field=True), Poly(x*(x + 1), t, field=True), DE)
    y = ya.as_expr()/yd.as_expr()
    assert cancel(derivation(Poly(y, t, field=True), DE).as_expr() -
        y/(x*(x + 1)) - (x + 1 - t)/(x*(x + 1))) == 0


def test_bound_degree_undecidable():
    # Exp case with deg(a) == deg(b) and alpha == -lc(b)/lc(a) == 1/(x + 1):
    # deciding whether alpha == m*Dt/t + Dz/z requires log(x + 1), which is
    # not in the tower, so parametric_log_deriv() cannot decide and
    # bound_degree() must raise rather than return an unsound bound.
    # param_rischDE() propagates this (it used to guess n = 5 instead,
    # which could truncate the parametric solution basis and turn into a
    # false proof of nonelementarity downstream in is_deriv_in_field()).
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t0),
        Poly(t1, t1)], 'exts': ['log', 'exp'], 'extargs': [x, x]})
    raises(NotImplementedError, lambda: bound_degree(Poly(t1, t1),
        Poly(-t1/(x + 1), t1, field=True), Poly(1, t1), DE))
    raises(NotImplementedError, lambda: bound_degree(Poly(t1, t1),
        Poly(-t1/(x + 1), t1, field=True), [Poly(1, t1)], DE, parametric=True))


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
    # n == oo (from rischDE() when bound_degree() cannot decide) must not
    # loop forever when deg(a) > 0: a == t, b == 1, c == x has no
    # solution and gcd(a, b) == 1 stays trivial, so no other exit is
    # ever taken.  This used to spin indefinitely.
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    raises(NotImplementedError, lambda: spde(Poly(t, t), Poly(1, t),
        Poly(x, t, domain='ZZ[x]'), oo, DE))
    # ...while a solvable positive-degree case with n == oo still reduces
    # and returns (q == 1 here; c becomes zero after one pass)...
    assert spde(Poly(t, t), Poly(1, t), Poly(1, t), oo, DE) == \
        (Poly(0, t), Poly(0, t), 0, Poly(0, t), Poly(1, t, domain='QQ'))
    # A SymPy Integer bound must behave like the equal plain int
    raises(NonElementaryIntegralException, lambda: spde(Poly(t, t),
        Poly(1, t), Poly(x, t, domain='ZZ[x]'), S(1), DE))
    # ...and with deg(a) == 0 it returns immediately, so n == oo is fine.
    assert spde(Poly(2, t), Poly(t, t), Poly(t, t), oo, DE) == \
        (Poly(t/2, t, domain='QQ'), Poly(t/2, t, domain='QQ'), oo,
         Poly(1, t, domain='ZZ'), Poly(0, t, domain='ZZ'))

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
    # Immediate u == 0 return: -lc(b)/lc(Dt) == 3 is hit on the first
    # iteration, reducing to a smaller problem of degree at most 3
    assert no_cancel_equal(Poly(-3*t, t), Poly(t**2 + 1, t), oo, DE) == \
        (Poly(0, t), 3, Poly(t**2 + 1, t))


def test_solve_poly_rde_cancel():
    # exp
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    assert cancel_exp(Poly(2*x, t), Poly(2*x, t), 0, DE) == \
        Poly(1, t)
    assert cancel_exp(Poly(2*x, t), Poly((1 + 2*x)*t, t), 1, DE) == \
        Poly(t, t)
    # The is_deriv_in_field() branch: b == (2*x + 1)/x == D(x)/x + 2*Dt/t
    # (parametric_log_deriv() gives (1, 2, x)), and Dq + b*q == c has the
    # solution q == 1 for c == (2*x + 1)/x.
    assert cancel_exp(Poly((2*x + 1)/x, t), Poly((2*x + 1)/x, t), 0, DE) == \
        Poly(1, t)
    # ... but 1/(2*x) also works for c == 1/x (the code finds solutions
    # with denominators via p/(z*t**m)):
    assert cancel_exp(Poly((2*x + 1)/x, t), Poly(1/x, t), 0, DE) == \
        Poly(1/(2*x), t)
    # c == 1/x**2 requires a' + 2*a == 1/x over QQ(x), which has no
    # rational solution (pole order mismatch at 0), so no solution exists
    raises(NonElementaryIntegralException, lambda: cancel_exp(
        Poly((2*x + 1)/x, t), Poly(1/x**2, t), 0, DE))

    # The t**(-m) stripping in the m < 0 branch is only valid when the
    # stripped coefficient is C/z for a constant C.  Here b == 1/x - 2 ==
    # D(x)/x - 2*Dt/t (z == x, m == -2), and the candidate q == 1 + t**2
    # has t**2-coefficient 1, with z*1 == x nonconstant, so no choice of
    # the integration constant removes the term: there is no solution
    # with n == 0.
    raises(NonElementaryIntegralException, lambda: cancel_exp(
        Poly(1/x - 2, t, field=True),
        Poly(1/x - 2 + t**2/x, t, field=True), 0, DE))
    # ...while the same b with c == b has the genuine solution q == 1.
    assert cancel_exp(Poly(1/x - 2, t, field=True),
        Poly(1/x - 2, t, field=True), 0, DE) == Poly(1, t)

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

    # The is_deriv_in_field() branch: b == 1/x == D(x)/x
    # (is_log_deriv_k_t_radical_in_field() gives (1, x)), and
    # Dq + b*q == c has the solution q == t for c == (t + 1)/x.
    assert cancel_primitive(Poly(1/x, t), Poly((t + 1)/x, t), 1, DE) == \
        Poly(t, t)
    # z*c == t/(x + 1) has antiderivative t*log(x + 1) - log(x + 1) + x,
    # which is not in QQ(x, log(x)), so no solution exists
    raises(NonElementaryIntegralException, lambda: cancel_primitive(
        Poly(1/x, t), Poly(t/(x*(x + 1)), t), 1, DE))

    # The degree bound n must not be clobbered by the index returned from
    # is_log_deriv_k_t_radical_in_field(): b == 1/(2*x) has index 2
    # (2*b == D(x)/x), which used to replace n == 3 and wrongly reject
    # the degree 3 solution q == t**3 of Dq + b*q == c.
    b = Poly(Rational(1, 2)/x, t)
    q = Poly(t**3, t)
    c = derivation(q, DE) + b*q
    assert cancel_primitive(b, c, 3, DE) == q

    # The b == 0 cancellation case (Dq == c by in-field integration):
    # primitive (t == log(x)): D(t**2) == 2*t/x
    assert solve_poly_rde(Poly(0, t), Poly(2*t/x, t), 2, DE) == Poly(t**2, t)
    # (note Dq == 1/x has the in-field solution q == t == log(x); use an
    # integrand whose antiderivative needs log(x + 1), which is not in
    # the field)
    raises(NonElementaryIntegralException, lambda: solve_poly_rde(
        Poly(0, t), Poly(t/(x + 1), t), 2, DE))
    # exp (t == exp(x)): D(x*t) == (x + 1)*t
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    assert solve_poly_rde(Poly(0, t), Poly((x + 1)*t, t), 1, DE) == \
        Poly(x*t, t)
    raises(NonElementaryIntegralException, lambda: solve_poly_rde(
        Poly(0, t), Poly(t/x, t), 1, DE))


def test_cancel_tan():
    # t = tan(x)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t**2 + 1, t)]})
    # Example 6.6.1: Dq + (1 - t)*q == -2*(x + 1)*t - 2*x with n == 1
    # (b0 == 1) has the solution q == -2*x*t
    assert cancel_tan(Poly(1, t), Poly(-2*(x + 1)*t - 2*x, t), 1, DE) == \
        Poly(-2*x*t, t)
    # ... reached from solve_poly_rde() via no_cancel_equal()'s
    # hand-off (Example 6.5.3), and from rischDE(): equation (6.21),
    # Dy + (1 - t)*y == t**3 + t**2 - 2*x*t - 2*x, has the solution
    # y == t**2 - 2*x*t
    assert solve_poly_rde(Poly(1 - t, t), Poly(t**3 + t**2 - 2*x*t - 2*x, t),
        2, DE) == Poly(t**2 - 2*x*t, t)
    ya, yd = rischDE(Poly(1 - t, t), Poly(1, t),
        Poly(t**3 + t**2 - 2*x*t - 2*x, t), Poly(1, t), DE)
    assert cancel(ya.as_expr()/yd.as_expr()) == t**2 - 2*x*t
    # n == 0 with b0 != 0 is a Risch differential equation over k
    assert cancel_tan(Poly(1, t), Poly(x + 1, t), 0, DE) == Poly(x, t)
    # n == 0 with b0 == 0 is in-field integration in k
    assert cancel_tan(Poly(0, t), Poly(2*x, t), 0, DE) == Poly(x**2, t)
    raises(NonElementaryIntegralException, lambda: cancel_tan(
        Poly(0, t), Poly(1/x, t), 0, DE))
    # The book's pseudocode returns u + v*t unverified when n == 1,
    # which is a wrong answer here: the t**2 coefficient of
    # Dq + (1 - t)*q vanishes identically for deg(q) <= 1, so
    # c == t**2 + 1 has no solution, while the projected coupled
    # system is solvable (by (u, v) == (0, 0), so the book returns
    # q == 0 with Dq + (1 - t)*q == 0 != c)
    raises(NonElementaryIntegralException, lambda: cancel_tan(
        Poly(1, t), Poly(t**2 + 1, t), 1, DE))
    # the degree bound appears in the equation itself, so n == oo is
    # not allowed
    raises(ValueError, lambda: cancel_tan(Poly(1, t), Poly(1, t), oo, DE))
    # A nonconstant eta (Dt == x*(t**2 + 1)) with a constant
    # -lc(b)/lc(Dt) must pass the symbolic-ratio guard in
    # solve_poly_rde() (which used to test lc(b) itself):
    # Dq + (1 - x*t)*q == x*t**3 + t**2 + 2*x*t has the solution t**2
    DE = DifferentialExtension(extension={'D': [Poly(1, x),
        Poly(x*(t**2 + 1), t)]})
    assert solve_poly_rde(Poly(1 - x*t, t), Poly(x*t**3 + t**2 + 2*x*t, t),
        2, DE) == Poly(t**2, t)
    # ... including when the ratio only becomes a number after
    # cancellation ((2*x + 2)/(x + 1) does not auto-simplify, so the
    # uncancelled dispatch comparison used to raise TypeError)
    DE = DifferentialExtension(extension={'D': [Poly(1, x),
        Poly((x + 1)*(t**2 + 1), t)]})
    b = Poly(1 - (2*x + 2)*t, t)
    q0 = Poly(t**2, t)
    c = derivation(q0, DE) + b*q0
    assert solve_poly_rde(b, c, 2, DE) == q0
    # A ratio that is a nonconstant element of k ((x**2 - 1)/(x - 1)
    # == x + 1) never equals an integer degree, so no cancellation
    # occurs and no_cancel_equal() decides the equation (this also
    # used to raise TypeError from the dispatch comparison)
    DE = DifferentialExtension(extension={'D': [Poly(1, x),
        Poly((x - 1)*(t**2 + 1), t)]})
    b = Poly(1 - (x**2 - 1)*t, t)
    q0 = Poly(t, t)
    c = derivation(q0, DE) + b*q0
    assert solve_poly_rde(b, c, 1, DE) == q0
    # A numeric ratio that is not a positive integer (here 3/2) can
    # never be a cancellation degree, so no_cancel_equal() decides
    # the equation even with n below the ratio (this used to fall
    # into the cancellation dispatch and raise NotImplementedError)
    DE = DifferentialExtension(extension={'D': [Poly(1, x),
        Poly(t**2 + 1, t)]})
    b = Poly(1 - Rational(3, 2)*t, t, field=True)
    q0 = Poly(t, t, field=True)
    c = derivation(q0, DE) + b*q0
    assert solve_poly_rde(b, c, 1, DE).as_expr() == t
    # A positive-integer ratio above the degree bound is likewise
    # unattainable: Dq + (1 - 2*t)*q == -t**2 + t + 1 with n == 1 has
    # the solution t (this used to raise NotImplementedError from the
    # cancellation dispatch, and no_cancel_equal() used to turn it
    # into a false no-solution proof)
    b = Poly(1 - 2*t, t)
    assert solve_poly_rde(b, Poly(-t**2 + t + 1, t), 1, DE).as_expr() == t
    # A non-real numeric ratio (here sqrt(-1)) admits no ordered
    # comparison with n but can never be a cancellation degree
    b = Poly(1 - I*t, t)
    q0 = Poly(t, t)
    c = derivation(q0, DE) + b*q0
    assert solve_poly_rde(b, c, 1, DE).as_expr().expand() == t
    # A parameter-mixed ratio whose tower-dependent part survives
    # every specialization (x + k: the coefficient of x is 1) is
    # accepted...
    b = Poly(1 - (x + k)*t, t)
    q0 = Poly(t, t, field=True)
    c = derivation(q0, DE) + b*q0
    assert solve_poly_rde(b, c, 1, DE).as_expr() == t
    # ... while a ratio depending only on a parameter of the constant
    # field is still refused (the cancellation degree would be
    # parameter-dependent, making no_cancel_equal()'s divisions only
    # generically valid)
    raises(TypeError, lambda: solve_poly_rde(Poly(1 - k*t, t),
        Poly(1, t), 2, DE))


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

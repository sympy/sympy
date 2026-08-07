"""Most of these tests come from the examples in Bronstein's book."""
from __future__ import annotations
from sympy.integrals.risch import DifferentialExtension, derivation, frac_in
from sympy.integrals.prde import (prde_normal_denom, prde_special_denom,
    prde_linear_constraints, constant_system, prde_spde, prde_no_cancel_b_large,
    prde_no_cancel_b_small, prde_no_cancel_b_equal, limited_integrate_reduce,
    limited_integrate,
    is_deriv_k, is_log_deriv_k_t_radical, parametric_log_deriv_heu,
    parametric_log_deriv, parametric_log_deriv_structure,
    is_log_deriv_k_t_radical_in_field, param_poly_rischDE, param_rischDE,
    is_deriv_in_field,
    prde_cancel_liouvillian)

from sympy.polys.polymatrix import PolyMatrix as Matrix

from sympy.testing.pytest import raises

from sympy.core import Add
from sympy.core.numbers import Rational
from sympy.functions.elementary.exponential import exp
from sympy.core.singleton import S
from sympy.core.symbol import symbols
from sympy.polys.domains.rationalfield import QQ
from sympy.polys.polytools import Poly, cancel
from sympy.abc import x, t, n

t0, t1, t2, t3, k = symbols('t:4 k')


def test_prde_normal_denom():
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1 + t**2, t)]})
    fa = Poly(1, t)
    fd = Poly(x, t)
    G = [(Poly(t, t), Poly(1 + t**2, t)), (Poly(1, t), Poly(x + x*t**2, t))]
    assert prde_normal_denom(fa, fd, G, DE) == \
        (Poly(x, t, domain='ZZ(x)'), (Poly(1, t, domain='ZZ(x)'), Poly(1, t,
            domain='ZZ(x)')), [(Poly(x*t, t, domain='ZZ(x)'),
         Poly(t**2 + 1, t, domain='ZZ(x)')), (Poly(1, t, domain='ZZ(x)'),
             Poly(t**2 + 1, t, domain='ZZ(x)'))], Poly(1, t, domain='ZZ(x)'))
    G = [(Poly(t, t), Poly(t**2 + 2*t + 1, t)), (Poly(x*t, t),
        Poly(t**2 + 2*t + 1, t)), (Poly(x*t**2, t), Poly(t**2 + 2*t + 1, t))]
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    assert prde_normal_denom(Poly(x, t), Poly(1, t), G, DE) == \
        (Poly(t + 1, t), (Poly((-1 + x)*t + x, t), Poly(1, t, domain='ZZ[x]')), [(Poly(t, t),
        Poly(1, t)), (Poly(x*t, t), Poly(1, t, domain='ZZ[x]')), (Poly(x*t**2, t),
        Poly(1, t, domain='ZZ[x]'))], Poly(t + 1, t))


def test_prde_special_denom():
    a = Poly(t + 1, t)
    ba = Poly(t**2, t)
    bd = Poly(1, t)
    G = [(Poly(t, t), Poly(1, t)), (Poly(t**2, t), Poly(1, t)), (Poly(t**3, t), Poly(1, t))]
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    assert prde_special_denom(a, ba, bd, G, DE) == \
        (Poly(t + 1, t), Poly(t**2, t), [(Poly(t, t), Poly(1, t)),
        (Poly(t**2, t), Poly(1, t)), (Poly(t**3, t), Poly(1, t))], Poly(1, t))
    G = [(Poly(t, t), Poly(1, t)), (Poly(1, t), Poly(t, t))]
    assert prde_special_denom(Poly(1, t), Poly(t**2, t), Poly(1, t), G, DE) == \
        (Poly(1, t), Poly(t**2 - 1, t), [(Poly(t**2, t), Poly(1, t)),
        (Poly(1, t), Poly(1, t))], Poly(t, t))
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(-2*x*t0, t0)]})
    DE.decrement_level()
    G = [(Poly(t, t), Poly(t**2, t)), (Poly(2*t, t), Poly(t, t))]
    assert prde_special_denom(Poly(5*x*t + 1, t), Poly(t**2 + 2*x**3*t, t), Poly(t**3 + 2, t), G, DE) == \
        (Poly(5*x*t + 1, t), Poly(0, t, domain='ZZ[x]'), [(Poly(t, t), Poly(t**2, t)),
        (Poly(2*t, t), Poly(t, t))], Poly(1, x))
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly((t**2 + 1)*2*x, t)]})
    G = [(Poly(t + x, t), Poly(t*x, t)), (Poly(2*t, t), Poly(x**2, x))]
    assert prde_special_denom(Poly(5*x*t + 1, t), Poly(t**2 + 2*x**3*t, t), Poly(t**3, t), G, DE) == \
        (Poly(5*x*t + 1, t), Poly(0, t, domain='ZZ[x]'), [(Poly(t + x, t), Poly(x*t, t)),
        (Poly(2*t, t, x), Poly(x**2, t, x))], Poly(1, t))
    assert prde_special_denom(Poly(t + 1, t), Poly(t**2, t), Poly(t**3, t), G, DE) == \
        (Poly(t + 1, t), Poly(0, t, domain='ZZ[x]'), [(Poly(t + x, t), Poly(x*t, t)), (Poly(2*t, t, x),
        Poly(x**2, t, x))], Poly(1, t))

    # Hypertangent case with a nontrivial special denominator (parametric
    # version of the corresponding test in test_rde.py)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t**2 + 1, t)]})
    assert prde_special_denom(Poly(t, t), Poly(t, t), Poly(1, t),
        [(Poly(-2*t**2 + t, t), Poly(t**2 + 1, t))], DE) == \
        (Poly(t, t), Poly(-2*t**2 + t, t),
        [(Poly(-2*t**2 + t, t), Poly(1, t))], Poly(t**2 + 1, t))


def test_prde_linear_constraints():
    DE = DifferentialExtension(extension={'D': [Poly(1, x)]})
    G = [(Poly(2*x**3 + 3*x + 1, x), Poly(x**2 - 1, x)), (Poly(1, x), Poly(x - 1, x)),
        (Poly(1, x), Poly(x + 1, x))]
    assert prde_linear_constraints(Poly(1, x), Poly(0, x), G, DE) == \
        ((Poly(2*x, x, domain='QQ'), Poly(0, x, domain='QQ'), Poly(0, x, domain='QQ')),
            Matrix([[1, 1, -1], [5, 1, 1]], x))
    G = [(Poly(t, t), Poly(1, t)), (Poly(t**2, t), Poly(1, t)), (Poly(t**3, t), Poly(1, t))]
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    assert prde_linear_constraints(Poly(t + 1, t), Poly(t**2, t), G, DE) == \
        ((Poly(t, t, domain='QQ'), Poly(t**2, t, domain='QQ'), Poly(t**3, t, domain='QQ')),
            Matrix(0, 3, [], t))
    G = [(Poly(2*x, t), Poly(t, t)), (Poly(-x, t), Poly(t, t))]
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)]})
    assert prde_linear_constraints(Poly(1, t), Poly(0, t), G, DE) == \
        ((Poly(0, t, domain='QQ[x]'), Poly(0, t, domain='QQ[x]')), Matrix([[2*x, -x]], t))


def test_constant_system():
    A = Matrix([[-(x + 3)/(x - 1), (x + 1)/(x - 1), 1],
                [-x - 3, x + 1, x - 1],
                [2*(x + 3)/(x - 1), 0, 0]], t)
    u = Matrix([[(x + 1)/(x - 1)], [x + 1], [0]], t)
    DE = DifferentialExtension(extension={'D': [Poly(1, x)]})
    R = QQ.frac_field(x)[t]
    assert constant_system(A, u, DE) == \
        (Matrix([[1, 0, 0],
                 [0, 1, 0],
                 [0, 0, 0],
                 [0, 0, 1]], ring=R), Matrix([0, 1, 0, 0], ring=R))


def test_prde_spde():
    D = [Poly(x, t), Poly(-x*t, t)]
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)]})
    # TODO: when bound_degree() can handle this, test degree bound from that too
    assert prde_spde(Poly(t, t), Poly(-1/x, t), D, n, DE) == \
        (Poly(t, t), Poly(0, t, domain='ZZ(x)'),
        [Poly(2*x, t, domain='ZZ(x)'), Poly(-x, t, domain='ZZ(x)')],
        [Poly(-x**2, t, domain='ZZ(x)'), Poly(0, t, domain='ZZ(x)')], n - 1)


def test_prde_no_cancel():
    # b large
    DE = DifferentialExtension(extension={'D': [Poly(1, x)]})
    assert prde_no_cancel_b_large(Poly(1, x), [Poly(x**2, x), Poly(1, x)], 2, DE) == \
        ([Poly(x**2 - 2*x + 2, x), Poly(1, x)], Matrix([[1, 0, -1, 0],
                                                        [0, 1, 0, -1]], x))
    assert prde_no_cancel_b_large(Poly(1, x), [Poly(x**3, x), Poly(1, x)], 3, DE) == \
        ([Poly(x**3 - 3*x**2 + 6*x - 6, x), Poly(1, x)], Matrix([[1, 0, -1, 0],
                                                                 [0, 1, 0, -1]], x))
    assert prde_no_cancel_b_large(Poly(x, x), [Poly(x**2, x), Poly(1, x)], 1, DE) == \
        ([Poly(x, x, domain='ZZ'), Poly(0, x, domain='ZZ')], Matrix([[1, -1,  0,  0],
                                                                    [1,  0, -1,  0],
                                                                    [0,  1,  0, -1]], x))
    # b small
    # XXX: Is there a better example of a monomial with D.degree() > 2?
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t**3 + 1, t)]})

    # My original q was t**4 + t + 1, but this solution implies q == t**4
    # (c1 = 4), with some of the ci for the original q equal to 0.
    G = [Poly(t**6, t), Poly(x*t**5, t), Poly(t**3, t), Poly(x*t**2, t), Poly(1 + x, t)]
    R = QQ.frac_field(x)[t]
    assert prde_no_cancel_b_small(Poly(x*t, t), G, 4, DE) == \
        ([Poly(t**4/4 - x/12*t**3 + x**2/24*t**2 + (Rational(-11, 12) - x**3/24)*t + x/24, t),
        Poly(x/3*t**3 - x**2/6*t**2 + (Rational(-1, 3) + x**3/6)*t - x/6, t), Poly(t, t),
        Poly(0, t), Poly(0, t)], Matrix([[1, 0,              -1, 0, 0,  0,  0,  0,  0,  0],
                                         [0, 1, Rational(-1, 4), 0, 0,  0,  0,  0,  0,  0],
                                         [0, 0,               0, 0, 0,  0,  0,  0,  0,  0],
                                         [0, 0,               0, 1, 0,  0,  0,  0,  0,  0],
                                         [0, 0,               0, 0, 1,  0,  0,  0,  0,  0],
                                         [1, 0,               0, 0, 0, -1,  0,  0,  0,  0],
                                         [0, 1,               0, 0, 0,  0, -1,  0,  0,  0],
                                         [0, 0,               1, 0, 0,  0,  0, -1,  0,  0],
                                         [0, 0,               0, 1, 0,  0,  0,  0, -1,  0],
                                         [0, 0,               0, 0, 1,  0,  0,  0,  0, -1]], ring=R))

    # TODO: Add test for deg(b) <= 0 with b small
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1 + t**2, t)]})
    b = Poly(-1/x**2, t, field=True)  # deg(b) == 0
    q = [Poly(x**i*t**j, t, field=True) for i in range(2) for j in range(3)]
    h, A = prde_no_cancel_b_small(b, q, 3, DE)
    V = A.nullspace()
    R = QQ.frac_field(x)[t]
    assert len(V) == 1
    assert V[0] == Matrix([Rational(-1, 2), 0, 0, 1, 0, 0]*3, ring=R)
    assert (Matrix([h])*V[0][6:, :])[0] == Poly(x**2/2, t, domain='QQ(x)')
    assert (Matrix([q])*V[0][:6, :])[0] == Poly(x - S.Half, t, domain='QQ(x)')


def test_prde_no_cancel_b_equal():
    # deg(b) == delta(t) - 1, with t == tan-like (Dt == t**2 + 1)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t**2 + 1, t)]})
    # -lc(b)/lc(Dt) == -1 is not a positive integer, so the loop runs to
    # N == 0: Dq + t*q == 2*t**2 + 1 has the solution q == t
    assert prde_no_cancel_b_equal(Poly(t, t), [Poly(2*t**2 + 1, t)], 1, DE) == \
        ([Poly(t, t)], Matrix([[1, -1]], DE.t))
    # A solution with a nonconstant coefficient and two right hand sides:
    # q == t**2 solves Dq + t*q == c1 for c == D(t**2) + t*t**2
    b = Poly(t, t)
    q = Poly(t**2, t)
    c = derivation(q, DE) + b*q
    h, A = prde_no_cancel_b_equal(b, [c, Poly(0, t)], 2, DE)
    V = A.nullspace()
    sols = 0
    for v in V:
        if v[0] != 0:
            y = Add(*[(v[2 + j]*h[j]).as_expr() for j in range(len(h))])
            yp = Poly(y/v[0].as_expr(), t, field=True)
            if cancel(derivation(yp, DE).as_expr() + (b*yp).as_expr()
                    - c.as_expr()) == 0:
                sols += 1
    assert sols >= 1
    # When the possible-cancellation degree -lc(b)/lc(Dt) is reached, the
    # rest is delegated to the cancellation algorithms, which are not yet
    # implemented for delta(t) >= 2
    raises(NotImplementedError, lambda: prde_no_cancel_b_equal(
        Poly(-2*t, t), [Poly(t**4 + 3*t**2, t)], 3, DE))


def test_prde_cancel_liouvillian():
    ### 1. case == 'primitive'
    # used when integrating f = log(x) - log(x - 1)
    # Not taken from 'the' book
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)]})
    p0 = Poly(0, t, field=True)
    p1 = Poly((x - 1)*t, t, domain='ZZ(x)')
    p2 = Poly(x - 1, t, domain='ZZ(x)')
    p3 = Poly(-x**2 + x, t, domain='ZZ(x)')
    h, A = prde_cancel_liouvillian(Poly(-1/(x - 1), t), [Poly(-x + 1, t), Poly(1, t)], 1, DE)
    V = A.nullspace()
    assert h == [p0, p0, p1, p0, p0, p0, p0, p0, p0, p0, p2, p3, p0, p0, p0, p0]
    assert A.rank() == 16
    assert (Matrix([h])*V[0][:16, :]) == Matrix([[Poly(0, t, domain='QQ(x)')]])

    ### 2. case == 'exp'
    # used when integrating log(x/exp(x) + 1)
    # Not taken from book
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(-t, t)]})
    assert prde_cancel_liouvillian(Poly(0, t, domain='QQ[x]'), [Poly(1, t, domain='QQ(x)')], 0, DE) == \
            ([Poly(1, t, domain='QQ'), Poly(x, t, domain='ZZ(x)')], Matrix([[-1, 0, 1]], DE.t))

    ### 3. case == 'exp' with b != 0 and n > 0.  This exercises the level at
    ### which eta == Dt/t is computed and the sign of the residual update
    ### Fi == -(D(h) + b*h), both of which used to be wrong.
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    b = Poly(1/x, t, field=True)
    Q = [Poly((x + 1)*t/x, t, field=True)]
    h, A = prde_cancel_liouvillian(b, Q, 1, DE)
    V = A.nullspace()
    # Dy + y/x == c1*(x + 1)*t/x must have the solution y == t (c1 == 1)
    found = False
    for v in V:
        if v[0] == 0:
            continue
        y = Add(*[(v[1 + j]*h[j]).as_expr() for j in range(len(h))])/v[0].as_expr()
        yp = Poly(y, t, field=True)
        if cancel(derivation(yp, DE).as_expr() + b.as_expr()*y
                - Q[0].as_expr()) == 0:
            found = True
    assert found


def test_param_poly_rischDE():
    DE = DifferentialExtension(extension={'D': [Poly(1, x)]})
    a = Poly(x**2 - x, x, field=True)
    b = Poly(1, x, field=True)
    q = [Poly(x, x, field=True), Poly(x**2, x, field=True)]
    h, A = param_poly_rischDE(a, b, q, 3, DE)

    assert A.nullspace() == [Matrix([0, 1, 1, 1], DE.t)]  # c1, c2, d1, d2
    # Solution of a*Dp + b*p = c1*q1 + c2*q2 = q2 = x**2
    # is d1*h1 + d2*h2 = h1 + h2 = x.
    assert h[0] + h[1] == Poly(x, x, domain='QQ')
    # a*Dp + b*p = q1 = x has no solution.

    a = Poly(x**2 - x, x, field=True)
    b = Poly(x**2 - 5*x + 3, x, field=True)
    q = [Poly(1, x, field=True), Poly(x, x, field=True),
         Poly(x**2, x, field=True)]
    h, A = param_poly_rischDE(a, b, q, 3, DE)

    assert A.nullspace() == [Matrix([3, -5, 1, -5, 1, 1], DE.t)]
    p = -Poly(5, DE.t)*h[0] + h[1] + h[2]  # Poly(1, x)
    assert a*derivation(p, DE) + b*p == Poly(x**2 - 5*x + 3, x, domain='QQ')


def test_param_rischDE():
    DE = DifferentialExtension(extension={'D': [Poly(1, x)]})
    p1, px = Poly(1, x, field=True), Poly(x, x, field=True)
    G = [(p1, px), (p1, p1), (px, p1)]  # [1/x, 1, x]
    h, A = param_rischDE(-p1, Poly(x**2, x, field=True), G, DE)
    assert len(h) == 3
    p = [hi[0].as_expr()/hi[1].as_expr() for hi in h]
    V = A.nullspace()
    assert len(V) == 2
    assert V[0] == Matrix([-1, 1, 0, -1, 1, 0], DE.t)
    y = -p[0] + p[1] + 0*p[2]  # x
    assert y.diff(x) - y/x**2 == 1 - 1/x  # Dy + f*y == -G0 + G1 + 0*G2

    # the below test computation takes place while computing the integral
    # of 'f = log(log(x + exp(x)))'
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    G = [(Poly(t + x, t, domain='ZZ(x)'), Poly(1, t, domain='QQ')), (Poly(0, t, domain='QQ'), Poly(1, t, domain='QQ'))]
    h, A = param_rischDE(Poly(-t - 1, t, field=True), Poly(t + x, t, field=True), G, DE)
    assert len(h) == 5
    p = [hi[0].as_expr()/hi[1].as_expr() for hi in h]
    V = A.nullspace()
    assert len(V) == 3
    assert V[0] == Matrix([0, 0, 0, 0, 1, 0, 0], DE.t)
    y = 0*p[0] + 0*p[1] + 1*p[2] + 0*p[3] + 0*p[4]
    assert y.diff(t) - y/(t + x) == 0   # Dy + f*y = 0*G0 + 0*G1


def test_limited_integrate_reduce():
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)]})
    assert limited_integrate_reduce(Poly(x, t), Poly(t**2, t), [(Poly(x, t),
    Poly(t, t))], DE) == \
        (Poly(t, t), Poly(-1/x, t), Poly(t, t), 1, (Poly(x, t), Poly(1, t, domain='ZZ[x]')),
        [(Poly(-x*t, t), Poly(1, t, domain='ZZ[x]'))])

    # An exp case with a nontrivial special part (hs != 1).  This
    # distinguishes the first return component, which must be hn rather
    # than a == hn*hs (the book returns a, which does not satisfy the
    # stated contract): check the contract directly for the known
    # solution v == 1/(t*x), c1 == 3 of f == Dv + c1*w1.
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    fa, fd = frac_in((3*x**2*t - x - 1)/(x**2*t), t)
    G = [(Poly(1, t), Poly(1, t))]
    A, b, h, N, g, V = limited_integrate_reduce(fa, fd, G, DE)
    p = Poly(cancel(h.as_expr()/(t*x)), t, field=True)  # p == v*h
    assert p.degree(t) <= N
    lhs = (A*derivation(p, DE) + b*p).as_expr()
    rhs = cancel(g[0].as_expr()/g[1].as_expr() +
        3*V[0][0].as_expr()/V[0][1].as_expr())
    assert cancel(lhs - rhs) == 0


def test_limited_integrate():
    DE = DifferentialExtension(extension={'D': [Poly(1, x)]})
    G = [(Poly(x, x), Poly(x + 1, x))]
    assert limited_integrate(Poly(-(1 + x + 5*x**2 - 3*x**3), x),
    Poly(1 - x - x**2 + x**3, x), G, DE) == \
        ((Poly(x**2 - x + 2, x), Poly(x - 1, x, domain='QQ')), [2])
    G = [(Poly(1, x), Poly(x, x))]
    assert limited_integrate(Poly(5*x**2, x), Poly(3, x), G, DE) == \
        ((Poly(5*x**3/9, x), Poly(1, x, domain='QQ')), [0])
    # An empty list of special elements (the is_deriv_in_field() case)
    assert limited_integrate(Poly(2*x, x), Poly(1, x), [], DE) == \
        ((Poly(x**2, x), Poly(1, x, domain='QQ')), [])


def test_is_log_deriv_k_t_radical():
    DE = DifferentialExtension(extension={'D': [Poly(1, x)], 'exts': [],
        'extargs': []})
    assert is_log_deriv_k_t_radical(Poly(2*x, x), Poly(1, x), DE) is None

    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(2*t1, t1), Poly(1/x, t2)],
        'exts': ['exp', 'log'], 'extargs': [2*x, x]})
    assert is_log_deriv_k_t_radical(Poly(x + t2/2, t2), Poly(1, t2), DE) == \
        ([(t1, 1), (x, 1)], t1*x, 2, 0)
    # TODO: Add more tests

    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t0, t0), Poly(1/x, t)],
        'exts': ['exp', 'log'], 'extargs': [x, x]})
    assert is_log_deriv_k_t_radical(Poly(x + t/2 + 3, t), Poly(1, t), DE) == \
        ([(t0, 2), (x, 1)], x*t0**2, 2, 3)


def test_structure_theorem_guards():
    # A tower with an unlabeled top monomial (len(exts) < len(D)) is not
    # supported by the structure theorems.  When every primitive monomial
    # in it is a logarithm, it is reported as a nonelementary tower.
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t1),
        Poly(t2, t2)], 'exts': ['log'], 'extargs': [x]})
    for func in (is_deriv_k, is_log_deriv_k_t_radical):
        try:
            func(Poly(t2, t2), Poly(1, t2), DE)
        except NotImplementedError as e:
            assert "Nonelementary extensions" in str(e)
        else:
            raise AssertionError("NotImplementedError was not raised")

    # A primitive monomial that is not a logarithm (here an unlabeled
    # arctangent-like monomial) needs the real version of the structure
    # theorems instead.
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t1),
        Poly(1/(x**2 + 1), t2)], 'exts': ['log'], 'extargs': [x]})
    for func in (is_deriv_k, is_log_deriv_k_t_radical):
        try:
            func(Poly(t2, t2), Poly(1, t2), DE)
        except NotImplementedError as e:
            assert "hypertangent" in str(e)
        else:
            raise AssertionError("NotImplementedError was not raised")


def test_is_deriv_k():
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t1), Poly(1/(x + 1), t2)],
        'exts': ['log', 'log'], 'extargs': [x, x + 1]})
    assert is_deriv_k(Poly(2*x**2 + 2*x, t2), Poly(1, t2), DE) == \
        ([(t1, 1), (t2, 1)], t1 + t2, 2)

    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t1), Poly(t2, t2)],
        'exts': ['log', 'exp'], 'extargs': [x, x]})
    assert is_deriv_k(Poly(x**2*t2**3, t2), Poly(1, t2), DE) == \
        ([(x, 3), (t1, 2)], 2*t1 + 3*x, 1)
    # TODO: Add more tests, including ones with exponentials

    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(2/x, t1)],
        'exts': ['log'], 'extargs': [x**2]})
    assert is_deriv_k(Poly(x, t1), Poly(1, t1), DE) == \
        ([(t1, S.Half)], t1/2, 1)

    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(2/(1 + x), t0)],
        'exts': ['log'], 'extargs': [x**2 + 2*x + 1]})
    assert is_deriv_k(Poly(1 + x, t0), Poly(1, t0), DE) == \
        ([(t0, S.Half)], t0/2, 1)

    # Issue 10798
    # DE = DifferentialExtension(log(1/x), x)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(-1/x, t)],
        'exts': ['log'], 'extargs': [1/x]})
    assert is_deriv_k(Poly(1, t), Poly(x, t), DE) == ([(t, 1)], t, 1)


def test_is_log_deriv_k_t_radical_in_field():
    # NOTE: any potential constant factor in the second element of the result
    # doesn't matter, because it cancels in Da/a.
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)]})
    assert is_log_deriv_k_t_radical_in_field(Poly(5*t + 1, t), Poly(2*t*x, t), DE) == \
        (2, t*x**5)
    assert is_log_deriv_k_t_radical_in_field(Poly(2 + 3*t, t), Poly(5*x*t, t), DE) == \
        (5, x**3*t**2)

    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(-t/x**2, t)]})
    assert is_log_deriv_k_t_radical_in_field(Poly(-(1 + 2*t), t),
    Poly(2*x**2 + 2*x**2*t, t), DE) == \
        (2, t + t**2)
    assert is_log_deriv_k_t_radical_in_field(Poly(-1, t), Poly(x**2, t), DE) == \
        (1, t)
    assert is_log_deriv_k_t_radical_in_field(Poly(1, t), Poly(2*x**2, t), DE) == \
        (2, 1/t)


def test_parametric_log_deriv_structure():
    # The heuristic fails on all of these (z in k at a primitive level);
    # the structure-theorem method (equation (7.44)) decides them.
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)],
        'exts': ['log'], 'extargs': [x]})
    # f == 1/x == Dx/x, w == 1: 1*f == Dx/x + 0*w
    assert parametric_log_deriv_structure(Poly(1, t), Poly(x, t),
        Poly(1, t), Poly(1, t), DE) == (1, 0, x)
    # f == 1 + 1/(2*x), w == 1: 2*f == Dx/x + 2*w  (through the wrapper,
    # exercising the heuristic -> structure fallback)
    assert parametric_log_deriv(Poly(2*x + 1, t), Poly(2*x, t),
        Poly(1, t), Poly(1, t), DE) == (2, 2, x)
    # f == 1/(x + 1), w == 1 has the solution (1, 0, x + 1), but x + 1 is
    # not in the tower, so the structure method is inconclusive (None) and
    # the wrapper must raise rather than claim no solution exists.
    assert parametric_log_deriv_structure(Poly(1, t), Poly(x + 1, t),
        Poly(1, t), Poly(1, t), DE) is None
    raises(NotImplementedError, lambda: parametric_log_deriv(
        Poly(1, t), Poly(x + 1, t), Poly(1, t), Poly(1, t), DE))

    # Nontrivial n and m: f == 1/(2*x) + 3, w == 2:
    # 2*f == Dx/x + 3*w
    assert parametric_log_deriv_structure(Poly(6*x + 1, t), Poly(2*x, t),
        Poly(2, t), Poly(1, t), DE) == (2, 3, x)

    # Underdetermined system (w lies in the span of the generators) and a
    # w == 0 query.  Which solution is returned depends on the free
    # parameter choice, so check the defining equation n*f == Dv/v + m*w
    # instead of an exact tuple.
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t),
        Poly(t1, t1)], 'exts': ['log', 'exp'], 'extargs': [x, x]})
    for f, w in [(1/x + 2, S(3)), (1/x + 5, S.Zero)]:
        fa, fd = frac_in(f, t1)
        wa, wd = frac_in(w, t1)
        A = parametric_log_deriv_structure(fa, fd, wa, wd, DE)
        assert A is not None
        n, m, v = A
        assert n > 0 and n.is_Integer and m.is_Integer
        vv = v.subs(t1, exp(x))  # t1 == exp(x) in this extension
        assert cancel(n*f - m*w - vv.diff(x)/vv) == 0


def test_is_deriv_in_field():
    DE = DifferentialExtension(extension={'D': [Poly(1, x)]})
    assert is_deriv_in_field(Poly(2*x, x), Poly(1, x), DE) == \
        (Poly(x**2, x), Poly(1, x, domain='QQ'))
    assert is_deriv_in_field(Poly(1, x), Poly(x, x), DE) is None
    # t == log(x): D(x*t) == t + 1; 1/(x*t) == D(log(log(x))) is not in the field
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)]})
    A = is_deriv_in_field(Poly(t + 1, t), Poly(1, t), DE)
    assert A is not None
    va, vd = A
    vp = Poly(cancel(va.as_expr()/vd.as_expr()), t, field=True)
    assert cancel(derivation(vp, DE).as_expr() - (t + 1)) == 0
    assert is_deriv_in_field(Poly(1, t), Poly(x*t, t), DE) is None
    # t == exp(x): D(x*t**2) == (2*x + 1)*t**2; t is not a derivative
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    A = is_deriv_in_field(Poly((2*x + 1)*t**2, t), Poly(1, t), DE)
    assert A is not None
    va, vd = A
    vp = Poly(cancel(va.as_expr()/vd.as_expr()), t, field=True)
    assert cancel(derivation(vp, DE).as_expr() - (2*x + 1)*t**2) == 0


def test_parametric_log_deriv():
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)]})
    assert parametric_log_deriv_heu(Poly(5*t**2 + t - 6, t), Poly(2*x*t**2, t),
    Poly(-1, t), Poly(x*t**2, t), DE) == \
        (2, 6, t*x**5)

    DE = DifferentialExtension(extension={'D': [Poly(1, x)]})
    assert parametric_log_deriv_heu(Poly(-1, x), Poly(1, x), Poly(1, x),
    Poly(1, x), DE) ==\
        (1, -1, Poly(1, x))

    assert parametric_log_deriv_heu(Poly(-2, x), Poly(5, x), Poly(2, x),
    Poly(3, x), DE) ==\
        (5, -3, Poly(1, x))

    # f == w == 0: any n works, with m == 0 and v == 1 (the rational
    # shortcut used to raise ZeroDivisionError on this input, which
    # arises from the degenerate hypertangent case in test_special_denom())
    assert parametric_log_deriv_heu(Poly(0, x), Poly(1, x), Poly(0, x),
    Poly(1, x), DE) == \
        (1, 0, Poly(1, x))

    # Cases where z is in k, so the residue equations are vacuous.  These
    # used to return None ("proven no solution"), which was wrong: at the
    # base level they are completely decidable from the polynomial parts.
    # f == -1 + 1/(x + 1), w == 1
    assert parametric_log_deriv_heu(Poly(-x, x), Poly(x + 1, x), Poly(1, x),
    Poly(1, x), DE) == \
        (1, -1, x + 1)
    # f == 1/x, w == 1 (arises while integrating log(log(x + exp(x))))
    assert parametric_log_deriv_heu(Poly(1, x), Poly(x, x), Poly(1, x),
    Poly(1, x), DE) == \
        (1, 0, x)
    # f == 1/(x**2 + 1), w == 1: proven no solution (complex residues)
    assert parametric_log_deriv_heu(Poly(1, x), Poly(x**2 + 1, x), Poly(1, x),
    Poly(1, x), DE) is None

    # In a non-base field the z in k case cannot be decided by the
    # heuristic ("failed" in the book), so it must raise
    # NotImplementedError, not return None.
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)]})
    raises(NotImplementedError, lambda: parametric_log_deriv_heu(
        Poly(1, t), Poly(t + 1, t), Poly(1, t), Poly(1, t), DE))

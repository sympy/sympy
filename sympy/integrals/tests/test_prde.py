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
    is_deriv_in_field, is_deriv_k_atan, is_log_deriv_k_t_radical_tan,
    prde_cancel_liouvillian, prde_cancel_tan)

from sympy.polys.polymatrix import PolyMatrix as Matrix

from sympy.testing.pytest import raises

from sympy.core import Add, Dummy
from sympy.matrices import MutableDenseMatrix
from sympy.core.numbers import I, Rational, oo
from sympy.functions.elementary.exponential import exp
from sympy.core.singleton import S
from sympy.core.symbol import symbols
from sympy.polys.domains.rationalfield import QQ
from sympy.polys.polytools import Poly, cancel
from sympy.abc import x, t, n, y

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

    # Multiple rows with symbolic constants: pivoting on y is sound, since
    # y is a nonzero element of the constant field QQ(y)
    dum = Dummy()
    A = Matrix([[y, 1, x], [0, y, 1 + x]], dum)
    u = Matrix([[x + y], [x + 2]], dum)
    B, v = constant_system(A, u, DE)
    assert (B.to_Matrix(), v.to_Matrix()) == \
        (MutableDenseMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
         MutableDenseMatrix([[(y**2 - 1)/y**2], [1/y], [1]]))

    # An entry whose derivation-quotient rows still contain tower
    # variables, requiring more than one pass of the elimination loop
    # ("while A is not constant" in the book's ConstantSystem); the
    # output must be fully constant
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t0, t0),
        Poly((t0 + x*t0)*t1, t1)], 'exts': ['exp', 'exp'],
        'extargs': [x, x*t0]})
    A = Matrix([[t1, t0]], dum)
    u = Matrix([[0]], dum)
    B, v = constant_system(A, u, DE)
    assert not any(B.to_Matrix()[i, j].has(t0, t1)
        for i in range(B.rows) for j in range(B.cols))


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
    # rest is delegated to the cancellation algorithms (this used to
    # raise NotImplementedError before prde_cancel_tan() existed):
    # Dq - 2*t*q == c1*(t**4 + 3*t**2) has the solutions
    # q == c1*t**3 + d*(t**2 + 1) (t**2 + 1 solves the homogeneous
    # equation)
    b = Poly(-2*t, t)
    Q = [Poly(t**4 + 3*t**2, t)]
    h, A = prde_no_cancel_b_equal(b, Q, 3, DE)
    matched = set()
    for v in A.nullspace():
        c1 = v[0].as_expr()
        y = Poly(Add(*[(v[1 + j]*h[j]).as_expr() for j in range(len(h))]),
            t, field=True)
        assert cancel((derivation(y, DE) + b*y).as_expr()
            - c1*(t**4 + 3*t**2)) == 0
        matched.add((c1, y.as_expr()))
    assert (1, t**3) in matched
    assert (0, t**2 + 1) in matched
    # A cancellation degree above the bound n is unattainable, so the
    # descent runs through N == 0 instead of handing off with a bound
    # above the caller's: Dq + (1 - 2*t)*q == c1*(-t**2 + t + 1) with
    # n == 1 has the solutions q == c1*t
    b = Poly(1 - 2*t, t)
    Q = [Poly(-t**2 + t + 1, t)]
    h, A = prde_no_cancel_b_equal(b, Q, 1, DE)
    matched = set()
    for v in A.nullspace():
        c1 = v[0].as_expr()
        y = Poly(Add(*[(v[1 + j]*h[j]).as_expr() for j in range(len(h))]),
            t, field=True)
        assert cancel((derivation(y, DE) + b*y).as_expr()
            - c1*(-t**2 + t + 1)) == 0
        matched.add((c1, y.as_expr()))
    assert (1, t) in matched


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


def test_prde_cancel_tan():
    # t = tan(x)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t**2 + 1, t)]})
    b = Poly(1 - t, t)
    Q = [Poly(t**3 + t**2 - 2*x*t - 2*x, t)]

    def solutions(h, A, b, Q):
        # Every nullspace vector of A must produce a solution of
        # Dq + b*q == Sum(ci*qi); return the (c, q) pairs.
        sols = []
        for v in A.nullspace():
            cs = [v[i].as_expr() for i in range(len(Q))]
            q = Poly(Add(*[(v[len(Q) + j]*h[j]).as_expr()
                for j in range(len(h))]), t, field=True)
            lhs = derivation(q, DE) + b*q
            rhs = Add(*[ci*qi.as_expr() for ci, qi in zip(cs, Q)])
            assert cancel(lhs.as_expr() - rhs) == 0
            sols.append((cs, q))
        return sols

    # The parametric version of Example 6.6.1, reached through the
    # full param_poly_rischDE() hand-off (the parametric Example
    # 6.5.3): Dq + (1 - t)*q == c1*(t**3 + t**2 - 2*x*t - 2*x) with
    # n == 2 has exactly the solutions q == c1*(t**2 - 2*x*t)
    h, A = param_poly_rischDE(Poly(1, t), b, Q, 2, DE)
    sols = solutions(h, A, b, Q)
    assert any(cs == [1] and q.as_expr() == t**2 - 2*x*t for cs, q in sols)

    # Homogeneous solutions are found: D(t**2 + 1) == 2*t*(t**2 + 1),
    # so q == t**2 + 1 solves Dq - 2*t*q == 0 (b0 == 0, n == 2)
    b0 = Poly(0, t)
    h, A = prde_cancel_tan(b0, [Poly(0, t)], 2, DE)
    sols = solutions(h, A, Poly(-2*t, t), [Poly(0, t)])
    assert any(q.as_expr() == t**2 + 1 for cs, q in sols)

    # n == 0: q in k, so the t-coefficients of the right hand side
    # must vanish: Dq + q == c1*(x + 1) + c2*t forces c2 == 0, with
    # the solution q == c1*x
    h, A = prde_cancel_tan(Poly(1, t), [Poly(x + 1, t), Poly(t, t)], 0, DE)
    sols = solutions(h, A, Poly(1, t), [Poly(x + 1, t), Poly(t, t)])
    assert sols and all(cs[1] == 0 for cs, q in sols)
    assert any(cs == [1, 0] and q.as_expr() == x for cs, q in sols)

    # The degree bound appears in the equation itself, so n == oo is
    # not allowed
    raises(ValueError, lambda: prde_cancel_tan(Poly(1, t), [Poly(1, t)],
        oo, DE))


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
    # ... and with a symbolic constant coefficient
    assert limited_integrate(Poly(2*y*x, x), Poly(1, x), [], DE) == \
        ((Poly(y*x**2, x), Poly(1, x, domain='QQ')), [])


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


    # A tower with a symbolic constant parameter: t0 == exp(y*x).
    # exp(2*y*x) == t0**2 is a K-radical; the structure coefficients are
    # rational even though the constant field is QQ(y)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(y*t0, t0)],
        'exts': ['exp'], 'extargs': [y*x]})
    assert is_log_deriv_k_t_radical(Poly(2*y*x, t0), Poly(1, t0), DE) == \
        ([(t0, 2)], t0**2, 1, 0)

    # A constant is the logarithmic derivative of the radical 1
    DE = DifferentialExtension(extension={'D': [Poly(1, x)], 'exts': [],
        'extargs': []})
    assert is_log_deriv_k_t_radical(Poly(2, x), Poly(1, x), DE) == \
        ([], 1, 1, 2)

    # Real towers (Corollary 9.3.2 (ii)): hypertangent and arc-tangent
    # monomials never contribute.  t0 == exp(x), t1 == tan(x),
    # t2 == atan(x); exp(2*x) == t0**2
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t0, t0),
        Poly(1 + t1**2, t1), Poly(1/(x**2 + 1), t2)],
        'exts': ['exp', 'tan', 'atan'], 'extargs': [x, x, x]})
    assert is_log_deriv_k_t_radical(Poly(2*x, t2), Poly(1, t2), DE) == \
        ([(t0, 2)], t0**2, 1, 0)
    # ... but exp(atan(x)) is a new monomial over it
    assert is_log_deriv_k_t_radical(Poly(t2, t2), Poly(1, t2), DE) is None


def test_structure_theorem_guards():
    # The structure theorems need every monomial of the tower labeled
    # (len(exts) == len(D) - 1)
    for D in ([Poly(1, x), Poly(1/x, t1), Poly(t2, t2)],
              [Poly(1, x), Poly(1/x, t1), Poly(1/(x**2 + 1), t2)]):
        DE = DifferentialExtension(extension={'D': D, 'exts': ['log'],
            'extargs': [x]})
        for func in (is_deriv_k, is_log_deriv_k_t_radical, is_deriv_k_atan,
                is_log_deriv_k_t_radical_tan):
            try:
                func(Poly(t2, t2), Poly(1, t2), DE)
            except NotImplementedError as e:
                assert "labeled" in str(e)
            else:
                raise AssertionError("NotImplementedError was not raised")

    # The real version of the structure theorems, needed as soon as the
    # tower has a hypertangent or arc-tangent monomial, requires sqrt(-1)
    # not in the field: neither in the input ...
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t1),
        Poly(1 + t2**2, t2)], 'exts': ['log', 'tan'], 'extargs': [x, x]})
    for func in (is_deriv_k, is_log_deriv_k_t_radical, is_deriv_k_atan,
            is_log_deriv_k_t_radical_tan):
        try:
            func(Poly(I*x, t2), Poly(1, t2), DE)
        except NotImplementedError as e:
            assert "sqrt(-1)" in str(e)
        else:
            raise AssertionError("NotImplementedError was not raised")
    # ... nor in the tower (t2 == tan(I*x))
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t1),
        Poly(I*(1 + t2**2), t2)], 'exts': ['log', 'tan'], 'extargs': [x, I*x]})
    raises(NotImplementedError, lambda: is_deriv_k(Poly(x, t2), Poly(1, t2), DE))
    # whereas with only exponentials and logarithms, I is a fine constant
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(I*t1, t1)],
        'exts': ['exp'], 'extargs': [I*x]})
    assert is_deriv_k(Poly(t1, t1), Poly(1, t1), DE) == ([(I*x, 1)], I*x, 1)


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


    # Linearly dependent tower generators (log(x) and log(x**2)) make the
    # structure system underdetermined; the solution must still be a
    # correct full-length coefficient vector (the reduced system used to
    # be misread as a solution vector and silently zip-truncated)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t0),
        Poly(2/x, t1)], 'exts': ['log', 'log'],
        'extargs': [x, x**2]})
    assert is_deriv_k(Poly(x, t0), Poly(1, t0), DE) == \
        ([(t0, 1), (t1, 0)], t0, 1)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(2/x, t0),
        Poly(1/x, t1)], 'exts': ['log', 'log'],
        'extargs': [x**2, x]})
    assert is_deriv_k(Poly(x, t0), Poly(1, t0), DE) == \
        ([(t0, S.Half), (t1, 0)], t0/2, 1)

    # A constant is the derivative of 0
    DE = DifferentialExtension(extension={'D': [Poly(1, x)], 'exts': [],
        'extargs': []})
    assert is_deriv_k(Poly(2, x), Poly(1, x), DE) == ([], 0, 2)

    # Real towers (Corollary 9.3.2 (i)): hypertangent and arc-tangent
    # monomials never contribute.  t0 == exp(x), t1 == tan(x),
    # t2 == atan(x); log(exp(x)) == x
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t0, t0),
        Poly(1 + t1**2, t1), Poly(1/(x**2 + 1), t2)],
        'exts': ['exp', 'tan', 'atan'], 'extargs': [x, x, x]})
    assert is_deriv_k(Poly(t0, t2), Poly(1, t2), DE) == ([(x, 1)], x, 1)
    # ... but log(x) and log(tan(x)) are new monomials over it
    assert is_deriv_k(Poly(x, t2), Poly(1, t2), DE) is None
    assert is_deriv_k(Poly(t1, t2), Poly(1, t2), DE) is None


def test_is_deriv_k_atan():
    # atan(x) is a monomial over QQ(x) (Example 1 of Bronstein's 1989
    # paper, equation (4) with an empty sum)
    DE = DifferentialExtension(extension={'D': [Poly(1, x)], 'exts': [],
        'extargs': []})
    assert is_deriv_k_atan(Poly(x, x), Poly(1, x), DE) is None

    # t == atan(x): atan(x) == t and atan(2*x/(1 - x**2)) == 2*t, both
    # up to a constant
    DE = DifferentialExtension(extension={'D': [Poly(1, x),
        Poly(1/(x**2 + 1), t)], 'exts': ['atan'], 'extargs': [x]})
    assert is_deriv_k_atan(Poly(x, t), Poly(1, t), DE) == ([(t, 1)], t)
    assert is_deriv_k_atan(Poly(2*x, t), Poly(1 - x**2, t), DE) == \
        ([(t, 2)], 2*t)
    assert is_deriv_k_atan(Poly(x**2, t), Poly(1, t), DE) is None

    # t == tan(x): atan(tan(x)) == x up to a constant, and atan(x) is a
    # monomial over QQ(x, tan(x))
    DE = DifferentialExtension(extension={'D': [Poly(1, x),
        Poly(1 + t**2, t)], 'exts': ['tan'], 'extargs': [x]})
    assert is_deriv_k_atan(Poly(t, t), Poly(1, t), DE) == ([(x, 1)], x)
    assert is_deriv_k_atan(Poly(x, t), Poly(1, t), DE) is None

    # Exponentials and logarithms never contribute: atan(exp(x)) is a
    # monomial over QQ(x, exp(x))
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)],
        'exts': ['exp'], 'extargs': [x]})
    assert is_deriv_k_atan(Poly(t, t), Poly(1, t), DE) is None

    # Both kinds of monomials at once: t1 == atan(x), t2 == tan(x)
    # (the 1989 paper's Example 2): atan of the tangent addition
    # formula for tan(x + 2 - 3*atan(x)/4)
    DE = DifferentialExtension(extension={'D': [Poly(1, x),
        Poly(1/(x**2 + 1), t1), Poly(1 + t2**2, t2)],
        'exts': ['atan', 'tan'], 'extargs': [x, x]})
    assert is_deriv_k_atan(Poly(t2 + x, t2), Poly(1 - t2*x, t2), DE) == \
        ([(x, 1), (t1, 1)], t1 + x)


def test_is_log_deriv_k_t_radical_tan():
    # tan(x) is a monomial over QQ(x); tan(1) is a constant
    DE = DifferentialExtension(extension={'D': [Poly(1, x)], 'exts': [],
        'extargs': []})
    assert is_log_deriv_k_t_radical_tan(Poly(x, x), Poly(1, x), DE) is None
    assert is_log_deriv_k_t_radical_tan(Poly(1, x), Poly(1, x), DE) == \
        ([], 0, 1, 1)

    # Example 1 of Bronstein's 1989 paper: t == atan(x), and
    # tan(atan(x)/3) is algebraic of degree 3 over QQ(x, atan(x))
    DE = DifferentialExtension(extension={'D': [Poly(1, x),
        Poly(1/(x**2 + 1), t)], 'exts': ['atan'], 'extargs': [x]})
    assert is_log_deriv_k_t_radical_tan(Poly(t, t), Poly(3, t), DE) == \
        ([(t, 1)], t, 3, 0)
    # tan(2*atan(x)) == 2*x/(1 - x**2) is in the field (n == 1)
    assert is_log_deriv_k_t_radical_tan(Poly(2*t, t), Poly(1, t), DE) == \
        ([(t, 2)], 2*t, 1, 0)
    # tan(x) and tan(x*atan(x)) are monomials over QQ(x, atan(x))
    assert is_log_deriv_k_t_radical_tan(Poly(x, t), Poly(1, t), DE) is None
    assert is_log_deriv_k_t_radical_tan(Poly(x*t, t), Poly(1, t), DE) is None

    # Example 2 of the paper: t1 == atan(x), t2 == tan(x), and
    # tan(x + 2 - 3*atan(x)/4) is algebraic of degree 4 over
    # QQ(x, atan(x), tan(x), tan(2))
    DE = DifferentialExtension(extension={'D': [Poly(1, x),
        Poly(1/(x**2 + 1), t1), Poly(1 + t2**2, t2)],
        'exts': ['atan', 'tan'], 'extargs': [x, x]})
    assert is_log_deriv_k_t_radical_tan(Poly(x + 2 - Rational(3, 4)*t1, t2),
        Poly(1, t2), DE) == ([(x, 4), (t1, -3)], 4*x - 3*t1, 4, 2)

    # t == tan(x): tan(2*x + 1) is in QQ(x, tan(x), tan(1)) (n == 1)
    DE = DifferentialExtension(extension={'D': [Poly(1, x),
        Poly(1 + t**2, t)], 'exts': ['tan'], 'extargs': [x]})
    assert is_log_deriv_k_t_radical_tan(Poly(2*x + 1, t), Poly(1, t), DE) == \
        ([(x, 2)], 2*x, 1, 1)
    # but tan(x/2) is algebraic of degree 2 over it
    assert is_log_deriv_k_t_radical_tan(Poly(x, t), Poly(2, t), DE) == \
        ([(x, 1)], x, 2, 0)

    # Exponentials and logarithms never contribute: tan(log(x)) is a
    # monomial over QQ(x, log(x))
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)],
        'exts': ['log'], 'extargs': [x]})
    assert is_log_deriv_k_t_radical_tan(Poly(t, t), Poly(1, t), DE) is None


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

    # exp case: u must come back as an Expr, never as a malformed Poly
    # with the tower generator inside the coefficient domain (this used
    # to return Poly(t, x, domain='ZZ[t]')-shaped results, and issued
    # the deprecated Poly/Expr mixing warning when residue terms were
    # present)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    assert is_log_deriv_k_t_radical_in_field(Poly(1, t), Poly(1, t), DE) == \
        (1, t)
    # f == Dt/t + Dt/(t + 1): the log-derivative of t*(t + 1)
    assert is_log_deriv_k_t_radical_in_field(Poly(2*t + 1, t),
        Poly(t + 1, t), DE) == (1, t**2 + t)

    # Hypertangent case (Section 5.12), t == tan(x): after the residue
    # reduction p == a + b*t, and u == v*(t**2 + 1)**e with n*a == Dv/v
    # and e == n*b/(2*eta), eta == Dt/(t**2 + 1)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1 + t**2, t)]})
    # D(t**2 + 1)/(t**2 + 1) == 2*t
    assert is_log_deriv_k_t_radical_in_field(Poly(2*t, t), Poly(1, t), DE) == \
        (1, t**2 + 1)
    assert is_log_deriv_k_t_radical_in_field(Poly(t, t), Poly(1, t), DE) == \
        (2, t**2 + 1)
    # a != 0 recurses into k: D(x*(t**2 + 1))/(x*(t**2 + 1))
    assert is_log_deriv_k_t_radical_in_field(Poly(1 + 2*t*x, t), Poly(x, t),
        DE) == (1, t**2*x + x)
    # with a residue term: D(t - 1)/(t - 1) and D((t - 1)*(t**2 + 1))/...
    assert is_log_deriv_k_t_radical_in_field(Poly(1 + t**2, t), Poly(t - 1, t),
        DE) == (1, t - 1)
    assert is_log_deriv_k_t_radical_in_field(Poly(1 + t**2 + 2*t*(t - 1), t),
        Poly(t - 1, t), DE) == (1, t**3 - t**2 + t - 1)
    # 1 == Dv/v needs v == exp(x), and x*t needs e == n*x/2, neither in k(t)
    assert is_log_deriv_k_t_radical_in_field(Poly(1, t), Poly(1, t), DE) is None
    assert is_log_deriv_k_t_radical_in_field(Poly(x*t, t), Poly(1, t), DE) is None
    # the book's precondition sqrt(-1) not in k(t)
    raises(NotImplementedError, lambda: is_log_deriv_k_t_radical_in_field(
        Poly(I*t, t), Poly(1, t), DE))
    # nested hypertangents, t1 == tan(x), t2 == tan(tan(x)): the
    # recursion goes through the hypertangent case of the level below
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1 + t1**2, t1),
        Poly((1 + t1**2)*(1 + t2**2), t2)]})
    assert is_log_deriv_k_t_radical_in_field(Poly(2*t1 + 2*(1 + t1**2)*t2, t2),
        Poly(1, t2), DE) == (1, t1**2*t2**2 + t1**2 + t2**2 + 1)


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

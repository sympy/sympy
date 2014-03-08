"""
Algorithms for solving Coupled Differential System.

This method is used for solving Coupled Differential System.
Given a differential field K of characteristic 0 and f1, f2, g1, g2
in K, it decides whether the system of equations

     /     \        /          \    /    \       /    \
    |  Dy1  |      |  f1   af2  |  |  y1  |     |  g1  |
    |       |  +   |            |  |      |  =  |      |
    |  Dy2  |      |  f2   f1   |  |  y2  |     |  g2  |
     \     /        \          /    \    /       \    /

"""
from __future__ import print_function, division
from sympy import sqrt

from sympy.core import Dummy, oo

from sympy.polys import Poly, cancel

from sympy.integrals.risch import (NonElementaryIntegralException,
    frac_in, derivation, is_deriv, as_poly_1t, DecrementLevel)
from sympy.integrals.rde import (weak_normalizer,
    bound_degree, spde, normal_denom, special_denom)
from sympy.integrals.prde import (real_imag, is_log_deriv_k_t_radical_in_field)

def cds_cancel_primitive(a, b1, b2, c1, c2, DE, n):
    """
    Cancellation - primitive case

    Given a derivation D on k[t], n either an integer or positive
    infinity, a in Const(k), b1, b2 in k and c1, c2 in k[t] with
    Dt in k, sqrt(a) not in k(t) and b1 != 0 or b2 != 0, raises either
    "raise NonElementaryIntegralException" in which case the system as follows

         /     \       /        \    /  \        /    \
        |  Dq1  |     |  b1 ab2  |  | q1 |      |  c1  |
        |       |  +  |          |  |    |  =   |      |
        |  Dq2  |     |  b2  b1  |  | q2 |      |  c2  |
         \     /       \        /    \  /        \    /

    Equation 8.4 from Manuel Bronstien

    has no solution with both degrees at most n in k[t] or a solution
    q1, q2 in k[t] X k[t] of this system with deg(q1) <= n and deg(q2)
    <= n
    """
    t = DE.t
    k = Dummy('k')
    b1a, b1d = frac_in(b1, DE.t)
    b2a, b2d = frac_in(b2, DE.t)
    A1 = is_log_deriv_k_t_radical_in_field(b1a, b1d, DE)
    A2 = is_log_deriv_k_t_radical_in_field(b2a, b2d, DE)

    if A1 and A2:
        n1, u1 = A1
        n2, u2 = A2
        z1a, z1d = frac_in(cancel(u1*c1.as_expr() - u2*c2.as_expr()), DE.t)
        z2a, z2d = frac_in(cancel(u2*c1.as_expr() + u1*c2.as_expr()), DE.t)
        P1 = is_deriv(z1a, z1d, DE)
        P2 = is_deriv(z2a, z2d, DE)
        if P1 and P2:
            p1, _ = P1
            p2, _ = P2
            if Poly(p1).degree(t) <= n and Poly(p2).degree(t) <= n:
                num1 = Poly(cancel(u1*p1 + u2*p2), t)
                num2 = Poly(cancel(u1*p2 - u2*p1), t)
                denom = Poly(cancel(u1**2 + u2**2), t)
                return (cancel(num1/denom), cancel(num2/denom))
            else:
                raise NonElementaryIntegralException

    if c1 == 0 and c2 == 0:
        return (Poly(0, t), Poly(0, t))
    c1a, c1d = frac_in(c1, DE.t)
    c2a, c2d = frac_in(c2, DE.t)
    if n < max(c1a.degree(t) - c1d.degree(t), c2a.degree(t) - c2d.degree(t)):
        raise NonElementaryIntegralException
    q1, q2 = (Poly(0, t), Poly(0, t))

    while c1 or c2:
        c1a, c1d = frac_in(c1, DE.t)
        c2a, c2d = frac_in(c2, DE.t)
        m = max(c1a.degree(t) - c1d.degree(t), c2a.degree(t) - c2d.degree(t))
        if n < m:
            raise NonElementaryIntegralException

        c1k = as_poly_1t(c1, t, k).as_poly(t).nth(m)
        c2k = as_poly_1t(c2, t, k).as_poly(t).nth(m)
        with DecrementLevel(DE):
            (s1, s2) = coupled_DE_system(b1, b2, c1k, c2k, DE)
        (b1, b2) = (b1.as_poly(DE.t), b2.as_poly(DE.t))
        q1 = (q1 + s1.as_expr()*t**m).as_poly(DE.t)
        q2 = (q2 + s2.as_expr()*t**m).as_poly(DE.t)
        n = m - 1
        Ds1tm = derivation(s1*t**m, DE)
        Ds2tm = derivation(s2*t**m, DE)
        c1 = Poly(c1.as_expr() - Ds1tm.as_expr() - (b1*s1 - b2*s2).as_expr()*t**m, t)
        c2 = Poly(c2.as_expr() - Ds2tm.as_expr() - (b2*s1 + b1*s2).as_expr()*t**m, t)
    return (q1, q2)


def cds_cancel_exp(a, b1, b2, c1, c2, DE, n):
    """
    Cancellation - hyperexponential case

    Given a derivation D on k[t], n either an integer or +oo,
    a in Const(k), b1, b2 in k and c1, c2 in k[t] with Dt/t
    in k, sqrt(a) not in k(t) and b1 != 0 or b2 != 0, raises either
    "NonElementaryIntegralException", in which case the system as below

          /     \       /        \    /  \        /    \
         |  Dq1  |     |  b1 ab2  |  | q1 |      |  c1  |
         |       |  +  |          |  |    |  =   |      |
         |  Dq2  |     |  b2  b1  |  | q2 |      |  c2  |
          \     /       \        /    \  /        \    /

    has no solution with both degrees at most n in k[t], or a solution
    q1, q2 in k[t] X k[t] of this system with deg(q1)<=n & deg(g2)<=n
    """
    from sympy.integrals.prde import parametric_log_deriv
    k = Dummy('k')
    t = DE.t
    Dt = derivation(t, DE)
    wa, wd = frac_in(DE.d.quo(Poly(t, t)), t)

    with DecrementLevel(DE):
        wa, wd = frac_in(DE.d.quo(Poly(t, t)), DE.t)
        b1a, b1d = frac_in(b1, DE.t)
        b2a, b2d = frac_in(b2, DE.t)
        A1 = parametric_log_deriv(b1a, b1d, wa, wd, DE)
        A2 = parametric_log_deriv(b2a, b2d, wa, wd, DE)

    if A1 and A2:
        n1, m1, u1 = A1
        n2, m2, u2 = A2
        m = m1
        z1a, z1d = frac_in(cancel((u1.as_expr()*c1.as_expr() + a*u2.as_expr()*c2.as_expr())*t**m), t)
        z2a, z2d = frac_in(cancel((u2.as_expr()*c1.as_expr() + u1.as_expr()*c2.as_expr())*t**m), t)
        P1 = is_deriv(z1a, z1d, DE)
        P2 = is_deriv(z2a, z2d, DE)

        if P1 and P2:
            p1, i1 = P1
            p2, i2 = P2
            q1 = cancel((u1*p1 - a*u2*p2)*t**(-m)/(u1**2 - a*u2**2))
            q2 = cancel((u1*p2 - u2*p1)*t**(-m)/(u1**2 - a*u2**2))
            q1a, q1d = frac_in(q1, t)
            q2a, q2d = frac_in(q2, t)

            if q1d == Poly(1, DE.t) and q2d == Poly(2, t) and Poly(p1).degree(t) <= n and Poly(p2).degree(t) <= n:
                return (q1, q2)
            else:
                raise NonElementaryIntegralException

    if c1 == 0 and c2 == 0:
        return (Poly(0, t), Poly(0, t))
    c1a, c1d = frac_in(c1, DE.t)
    c2a, c2d = frac_in(c2, DE.t)

    if n < max(c1a.degree(t) - c1d.degree(t), c2a.degree(t) - c2d.degree(t)):
        raise NonElementaryIntegralException
    q1, q2 = (Poly(0, t), Poly(0, t))

    while c1 or c2:
        c1a, c1d = frac_in(c1, DE.t)
        c2a, c2d = frac_in(c2, DE.t)
        m = max(c1a.degree(t) - c1d.degree(t), c2a.degree(t) - c2d.degree(t))
        if n < m:
            raise NonElementaryIntegralException
        c1k = as_poly_1t(c1, t, k).as_poly(t).nth(m)
        c2k = as_poly_1t(c2, t, k).as_poly(t).nth(m)
        with DecrementLevel(DE):
            (s1, s2) = coupled_DE_system(b1 + m*Dt/t, b2, c1k, c2k, DE)
        (b1, b2) = (b1.as_poly(DE.t), b2.as_poly(DE.t))
        q1 = (q1 + s1.as_expr()*t**m).as_poly(DE.t)
        q2 = (q2 + s2.as_expr()*t**m).as_poly(DE.t)
        n = m - 1
        Ds1tm = derivation(s1*t**m, DE)
        Ds2tm = derivation(s2*t**m, DE)
        c1 = Poly(c1.as_expr() - Ds1tm.as_expr() - (b1*s1 + a*b2*s2).as_expr()*t**m, t)
        c2 = Poly(c2.as_expr() - Ds2tm.as_expr() - (b2*s1 + b1*s2).as_expr()*t**m, t)
    return (q1, q2)


def cds_cancel_tan(b0, b2, c1, c2, DE, n):
    """
    Cancellation - tangent case

    Given a derivation D on k[t], n either an integer or +oo,
    b0, b2 in k and c1, c2 in k[t] with Dt/(t**2 + 1) = eta in k, sqrt(-1)
    not in k(t) and b0 != 0 or b2 != 0, raises either "NonElementaryIntegralException",
    in which case the system
         /     \       /                           \    /    \     /  \
        |  Dq1  |     |  b0 - n*eta*t     -b2       |  |  q1  |   | c1 |
        |       |  +  |                             |  |      | = |    |
        |  Dq2  |     |      b2        b0 + n*eta*t |  |  q2  |   | c2 |
         \     /       \                           /    \    /     \  /
    has no solution with both degrees at most n in k[t], or a solution
    q1, q2 in k[t] X k[t] of this system with deg(q1)<=n and deg(q2)<=n
    """
    k = Dummy('k')
    t = DE.t
    if n == 0:
        if not c1.has(t) and not c2.has(t):
			A = coupled_DE_system(b0, b2, c1, c2)
        else:
            raise NonElementaryIntegralException

    a = sqrt(-1)
    p = t - a
    eta = DE.d.exquo(Poly(t**2 + 1, t))
    #u1 + u2*I = c1(I) + c2(I)*I
    c1_ = Poly(c1.eval(a), t)
    c2_ = Poly(c2.eval(a), t)
    ca, cd = frac_in(c1_ + c2_*a, t)
    u1a, u2a, u1d = real_imag(ca, cd, k)
    u1 = u1a.to_field().mul_ground(1/u1d)
    u2 = u2a.to_field().mul_ground(1/u1d)
    (s1, s2) = coupled_DE_system(b0, b2, u1, u2, DE)
    c = c1 - u1 + n*eta*(s1*t + s2) + (c2 - u2 + n*eta*(s2*t - s1))*a
    c = c.to_field().mul_ground(1/p)
    ca, cd = frac_in(c)
    d1a, d2a, d1d = real_imag(ca, cd, k)
    d1 = d1a.to_field().mul_ground(1/d1d)
    d2 = d2a.to_field().mul_ground(1/d1d)
    B = cds_cancel_tan(b0, b2 + eta, d1, d2, DE, n-1)
    h1, h2 = B
    return (h1*t + h2 + s1, h2*t - h1 + s2)


def non_cancellation_algo(alpha, beta, hn, hs, b, c, n, DE):
    from sympy.integrals.prde import real_imag
    from sympy.integrals.rde import (no_cancel_b_large,
        no_cancel_b_small, no_cancel_equal)
    k = Dummy('k')

    if not b.is_zero and (DE.case == 'base' or \
        b.degree(DE.t) > max(0, DE.d.degree(DE.t) - 1)):
        q = no_cancel_b_large(b, c, n, DE)
        qa , qd = frac_in(q, DE.t)
        qa_r, qa_i, qd = real_imag(alpha*qa + beta, hn*hs*qd, k)
        return (qa_r/qd, qa_i/qd)

    elif (b.is_zero or b.degree(DE.t) < DE.d.degree(DE.t) - 1) \
        and (DE.case == 'base' or DE.d.degree(DE.t) >= 2):
        q = no_cancel_b_small(b, c, n, DE)
        qa , qd = frac_in(q, DE.t)
        qa_r, qa_i, qd = real_imag(alpha*qa + beta, hn*hs*qd, k)
        return (qa_r/qd, qa_i/qd)

    elif DE.d.degree(DE.t) >= 2 and b.degree(DE.t) == DE.d.degree(DE.t)- 1 \
        and n > -b.as_poly(DE.t).LC()/DE.d.as_poly(DE.t).LC():
        q = no_cancel_equal(b, c, n, DE)
        qa , qd = frac_in(q, DE.t)
        qa_r, qa_i, qd = real_imag(alpha*qa + beta, hn*hs*qd, k)
        return (qa_r/qd, qa_i/qd)
    #spde passed but not a non-cancellation case
    return None


def cancellation_algo(a, b1, b2, c1, c2, DE, n):
    case = DE.case
    if case == 'base':
        raise NonElementaryIntegralException
    if case == 'primitive':
        (q1, q2) = cds_cancel_primitive(a, b1, b2, c1, c2, DE, n)
    if case == 'exp':
        (q1, q2) = cds_cancel_exp(a, b1, b2, c1, c2, DE, n)
    if case == 'tan':
        (q1, q2) = cds_cancel_tan(a, b1, b2, c1, c2, DE, n)
    return (q1, q2)


def coupled_DE_system(b1, b2, c1, c2, DE):
    """
    Algorithms for solving Coupled Differential System.

    This method is used for solving Coupled Differential System.
    Given a differential field K of characteristic 0 and f1, f2, g1, g2
    in K, it decides whether the system of equations
        /     \        /          \    /    \       /    \
       |  Dy1  |      |  b1   ab2  |  |  y1  |     |  c1  |
       |       |  +   |            |  |      |  =  |      |
       |  Dy2  |      |  b2   b1   |  |  y2  |     |  c2  |
        \     /        \          /    \    /       \    /

    Hence returning (y1, y2) if a solution exist, None otherwise
    """
    a = Poly(sqrt(-1) , DE.t)
    b1a, b1d = frac_in(b1, DE.t)
    b2a, b2d = frac_in(b2, DE.t)
    c1a, c1d = frac_in(c1, DE.t)
    c2a, c2d = frac_in(c2, DE.t)
    fa = b1a*b2d + b1d*b2a*sqrt(-1)
    fd = b1d*b2d
    ga = c1a*c2d + c1d*c2a*sqrt(-1)
    gd = c1d*c2d
    _, (fa, fd) = weak_normalizer(fa, fd, DE)
    a, (ba, bd), (ca, cd), hn = normal_denom(fa, fd, ga, gd, DE)
    A, B, C, hs = special_denom(a, ba, bd, ca, cd, DE)

    try:
       n = bound_degree(A, B, C, DE)
    except NotImplementedError:
       n = oo

    try:
       B, C, m, alpha, beta = spde(A, B, C, n, DE)
    except NonElementaryIntegralException:
       # Does not fall in non cancellation
       # Hence cancellation cases
       return cancellation_algo(a, b1, b2, c1, c2, DE, n)
    else:
       # non cancellation cases solve for q
       A = non_cancellation_algo(alpha, beta, hn, hs, B, C, m, DE)
       if A is None:
           return cancellation_algo(a, b1, b2, c1, c2, DE, n)
       return A

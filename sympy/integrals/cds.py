"""
Algorithms for solving Coupled Differential System.

This method is used for solving Coupled Differntial System.
Given a differntial field K of characterstic 0 and f1, f2, g1, g2
in K, it decides whether the system of equations
     /     \        /          \    /    \       /    \
    |  Dy1  |      |  f1   af2  |  |  y1  |     |  g1  |
    |       |  +   |            |  |      |  =  |      |
    |  Dy2  |      |  f2   f1   |  |  y2  |     |  g2  |
     \     /        \          /    \    /       \    /

Hence returning (y1, y2) if a solution exist, None otherwise
"""
from __future__ import with_statement

from sympy import sqrt

from sympy.core import Dummy, ilcm, Add, Mul, Pow, S

from sympy.polys import Poly, cancel, gcd

from sympy.integrals.risch import (frac_in, derivation, residue_reduce,
    residue_reduce_derivation, DecrementLevel)
from sympy.integrals.rde import (weak_normalizer,
    bound_degree, spde, solve_poly_rde, normal_denom, special_denom)
from sympy.integrals.prde import (real_imag)
def cds_cancel_prim(a, b1, b2, c1, c2, DE, n):
    """
    Cancellation - primitive case

    Given a derivation D on k[t], n either an integer or positive
    infinity, a in Const(k), b1, b2 in k and c1, c2 in k[t] with
    Dt in k, sqrt(a) not in k(t) and b1!=0 or b2=0, return either
    "no solution" in which case the system as follows

         /     \       /        \    /  \        /    \
        |  Dq1  |     |  b1 ab2  |  | q1 |      |  c1  |
        |       |  +  |          |  |    |  =   |      |
        |  Dq2  |     |  b2  b1  |  | q2 |      |  c1  |
         \     /       \        /    \  /        \    /

    Equation 8.4 from Manuel Bronstien

    has no solution with both degrees at most n in k[t] or a solution
    q1, q2 in k[t] X k[t] of this system with deg(q1) <= n and deg(q2)
    <= n
    """
    t = DE.t
    A1 = is_log_deriv_k_t_radical_in_field(b1)
    A2 = is_log_deriv_k_t_radical_in_field(b2)

    if A1 and A2:
        n1, u1 = A1
        n2, u2 = A2
        u = u1 + u2*sqrt(a)
        if is_deriv(u1*c1 + a*z2*c2) and is_deriv(z2*c1 + z1*c2):
            if p1.degree() <= n and p2.degree() <= n:
                denom = (z1**2 - a*z2**2).as_expr()
                return ((u1*p1 - a*z2*p2).as_expr()/denom, (z1*p2 - z2*p1).as_expr()/denom())
            else:
                return None
    if c1 == 0 and c2 == 0:
        return (0, 0)
    if n < max(c1.degree(), c2.degree()):
        return None
    q1, q2 = (0, 0)
    while c1 is not 0 or c2 is not 0:
        m = max(c1.degree(), c2.degree())
        if n < m:
            return None
        A = coupled_DE_System(b1, b2, Poly(c1.nth(m),DE.t), Poly(c2.nth(m),DE.t))
        if A is None:
            return None
        (s1, s2) = A
        q1 = q1 + s1*t**m
        q2 = q2 + s2*t**m
        n = m - 1
        Ds1tm = derivation(s1*t**m, DE.t)
        Ds2tm = derivation(s2*t**m, DE.t)
        c1 = c1 - Ds1tm - (b1*s1 + a*b2*s2)*t**m
        c2 = c2 - Ds2tm - (b2*s1 + b1*s2)*t**m
    return (q1, q2)


def cds_cancel_exp(a, b1, b2, c1, c2, DE, n):
    """
    Cancellation - hyperexponantial case

    Given a derivation D on k[t], n either an integer or positive
    infinity, a in Const(k), b1, b2 in k and c1, c2 in k[t] with Dt/t
    in k, sqrt(a) not in k(t) and b1 != 0  or b2 != 0, return either
    "no solution", in which case the system as below

          /     \       /        \    /  \        /    \
         |  Dq1  |     |  b1 ab2  |  | q1 |      |  c1  |
         |       |  +  |          |  |    |  =   |      |
         |  Dq2  |     |  b2  b1  |  | q2 |      |  c1  |
          \     /       \        /    \  /        \    /

    has no solution with both degrees at most n in k[t], or a solution
    q1, q2 in k[t] X k[t] of this system with deg(q1)<=n & deg(g2)<=n
    """
    from sympy.integrals.prde import parametric_log_deriv
    t = DE.t
    wa, wd = frac_in(DE.d.quo(Poly(DE.t, DE.t)), DE.t)
    b1a, b1d = frac_in(b1, DE.t)
    b2a, b2d = frac_in(b2, DE.t)
    A1 = parametric_log_deriv(b1a, b1d, wa, wd, DE)
    A2 = parametric_log_deriv(b2a, b2d, wa, wd, DE)
    print "printed"
    if A1 is not None and A2 is not None:
         n1, m1, u1 = A1
         n2, m2, u2 = A2
         m =  m1
         u = u1 + u2*sqrt(a)
         if is_deriv(u1*c1 + a*z2*c2) and is_deriv(z2*c1 + z1*c2):
             q1 = ((z1*p1 - a*z2*p2)*t**(-m)).as_expr()/(z1**2 - a*z2**2).as_expr()
             q2 = ((z1*p2 - z2*p1)*t**(-m)).as_expr()/(z1**2 - a*z2**2).as_expr()
             if q1 in k[t] and q2 in k[t] and q1.degree()<=n and p2.degree()<=n:
                 return (q1, q2)
             else:
                 return None
    if c1 == 0 and c2 == 0:
        return (0, 0)
    if n < max(c1.degree(), c2.degree()):
        return None
    q1, q2 = (0, 0)
    while c1 is not 0 and c2 is not 0:
         m = max(c1.degree(), c2.degree())
         if n <  m:
             return None
	 print b1, b2, Poly(c1.nth(m),DE.t), Poly(c2.nth(m),DE.t)
         A = coupled_DE_System(b1, b2, Poly(c1.nth(m),DE.t), Poly(c2.nth(m),DE.t), DE)
         if A is None:
             return None
         (s1, s2) = A
         q1 = q1 + s1*t**m
         q2 = q2 + s2*t**m
         n = m - 1
         Ds1tm = derivation(s1*t**m, DE.t)
         Ds2tm = derivation(s2*t**m, DE.t)
         c1 = c1 - Ds1tm - (b1*s1 + a*b2*s2)*t**m
         c2 = c2 - Ds2tm - (b2*s1 + b1*s2)*t**m
    return (q1, q2)


def cds_cancel_tan(b0, b2, c1, c2, DE, n):
    """
    Cancellation - tangent case

    Given a derivation D on k[t], n either an integer or positive infinity,
    b0, b2 in k and c1, c2 in k[t] with Dt/(t**2 + 1) = eta in k, sqrt(-1)
    not in k(t) and b0 != 0 or b2 != 0, return either "no solution", in
    which case the system
         /     \       /                           \    /    \     /  \
        |  Dq1  |     |  b0 - n*eta*t     -b2       |  |  q1  |   | c1 |
        |       |  +  |                             |  |      | = |    |
        |  Dq2  |     |      b2        b0 + n*eta*t |  |  q2  |   | c2 |
         \     /       \                           /    \    /     \  /
    has no solution with both degrees at most n in k[t], or a solution
    q1, q2 in k[t] X k[t] of this system with deg(q1)<=n and deg(q2)<=n
    """
    t = DE.t
    if n == 0:
      if c1 in k and c2 in k:
          A = coupled_DE_System(b0, b2, c1, c2)
      else:
          return None
    p = t - sqrt(-1)
    eta = DE.d.exquo(Poly(DE.t**2 + 1, DE.t))

    ca, cd = frac_in(c1 + c2*sqrt(-1), DE.t)

    ca_dict = ca.as_poly(sqrt(-1)).as_dict()
    ca_r = [value if key[0] % 2 == 0 else 0 for key, value in ca_dict.items()]
    ca_i = [value if key[0] % 2 == 1 else 0 for key, value in ca_dict.items()]
    ca_r = sum(r for r in  ca_r)
    ca_i = sum(r for r in  ca_i)
    cd_dict = cd.as_poly(sqrt(-1)).as_dict()
    cd_r = [value if key[0] % 2 == 0 else 0 for key, value in cd_dict.items()]
    cd_i = [value if key[0] % 2 == 1 else 0 for key, value in cd_dict.items()]
    cd_r = sum(r for r in  cd_r)
    cd_i = sum(r for r in  cd_i)

    u1a = ca_r*cd_r - ca_i*cd_i
    u1d = cd_r*cd_r + cd_i*cd_i
    u2a = ca_i*cd_r - ca_r*cd_i
    u2d = u1d

    u1 = (u1a.as_expr()/u1d.as_expr()).as_poly(DE.t)
    u2 = (u2a.as_expr()/u1d.as_expr()).as_poly(DE.t)
    A = coupled_DE_System(b0, b2, u1, u2, DE)
    if A is None:
        return None
    (s1, s2) = A
    c = c1 - z1 + n*eta*(s1*t + s2) + (c2 - z2 + n*eta*(s2*t - s1))*sqrt(-1)
    c = (c.as_expr()/p.as_expr()).as_poly(DE.t)
    d1a, d2a, d1d = real_imag(c)
    d1 = (d1a.as_expr()/d1d.as_expr()).as_poly(DE.t)
    d2 = (d2a.as_expr()/d1d.as_expr()).as_poly(DE.t)
    B = cds_cancel_tan(b0, b2 + eta, d1, d2, DE, n-1)
    if B is None:
        return None
    h1, h2 = B
    return (h1*t + h2 + s1, h2*t - h1 + s2)

def coupled_DE_System(b1, b2, c1, c2, DE):
    """
    """
<<<<<<< HEAD
    from sympy.integrals.prde import (prde_no_cancel_b_large,
        prde_no_cancel_b_small)
=======
    from sympy.integrals.rde import (no_cancel_b_large,
        no_cancel_b_small, no_cancel_equal)
>>>>>>> 2eed68e... added
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
    n = bound_degree(A, B, C, DE)
    B, C, m, alpha, beta = spde(A, B, C, n, DE)
    print "this is coupled"
    print B, C, m, alpha, beta
<<<<<<< HEAD
    if B !=0 and B.degree() >  max(0, DE.d.degree(DE.t) - 1):
        H, A = prde_no_cancel_b_large(B, [C], n, DE)
    elif Q.degree() > 0 and B.degree() < DE.d.degree(DE.t) - 1 and DE.d.degree(DE.t) >= 2:
        H, A = prde_no_cancel_b_small(B, [C], n, DE)
    print H, A
=======
    if B !=0 and B.degree() > max(0, DE.d.degree(DE.t) - 1):
        q = no_cancel_b_large(b, c, n, DE)
        return Poly(cancel((alpha*q + beta)/m) , DE.t)
    elif Q.degree() > 0 and B.degree() < DE.d.degree(DE.t) - 1 and DE.d.degree(DE.t) >= 2:
        q = no_cancel_b_small(b, c, n, DE)
        return Poly(cancel((alpha*q + beta)/m) , DE.t)
    elif DE.d.degree(DE.t) >=2 and B.degree() = DE.d.degree(DE.t) - 1:
        q = no_cancel_equal(b, c, n, DE)
        return Poly(cancel((alpha*q + beta)/m) , DE.t)
    cases = DE.case
    a = Poly(sqrt(-1) , DE.t)
    if case == 'primtive':
	(q1, q2) = cds_cancel_prim(a, b1, b2, c1, c2, DE, n)
    if case == 'exp':
	(q1, q2) = cds_cancel_exp(a, b1, b2, c1, c2, DE, n)
    if case == 'tan':
	(q1, q2) = cds_cancel_tan(a, b1, b2, c1, c2, DE, n)
    return Poly(cancel(((q1 + sqrt(-1)*q2)*alpha + beta)/m), DE.t)
>>>>>>> 2eed68e... added

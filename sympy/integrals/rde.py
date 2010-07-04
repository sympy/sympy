"""
Algorithms for solving the Risch differential equation.

Given a differential field K of characteristic 0 that is a simple monomial
extension of a base field k and f, g in K, the Risch Differential Equation
problem is to decide if there exist y in K such that Dy + f*y == g and to find
one if there are some.  If t is a monomial over k and the coefficients of f and
g are in k(t), then y is in k(t), and the outline of the algorithm here is
given as:

1. Compute the normal part n of the denominator of y.  The problem is then
reduced to finding y' in k<t>, where y == y'/n.
2. Compute the special part s of the denominator of y.   The problem is then
reduced to finding y'' in k[t], where y == y''/(n*s)
3. Bound the degree of y''.
4. Reduce the equation Dy + f*y == g to a similar equation with f, g in k[t].
5. Find the solutions in k[t] of bounded degree of the reduced equation.

See Chapter 6 of "Symbolic Integration I: Transcendental Functions" by Manuel
Bronstein.  See also the docstring of risch.py.
"""
from sympy.core.symbol import Symbol

from sympy.polys import Poly, gcd, ZZ, cancel

from sympy.integrals.risch import (gcdex_diophantine, derivation, splitfactor,
    NonElementaryIntegral, get_case)

from operator import mul
#    from pudb import set_trace; set_trace() # Debugging

# TODO: Add messages to NonElementaryIntegral errors
def order_at(a, p, t):
    """
    Computes the order of a at p, with respect to t.

    For a, p in k[t], the order of a at p is defined as nu_p(a) = max{n in Z+
    such that p**n|a}, where a != 0.  If a == 0, nu_p(a) = +oo, but this
    function will instead raise ValueError in that case.

    To compute the order at a rational function, a/b, use the fact that
    nu_p(a/b) == nu_p(a) - nu_p(b).
    """
    if a.is_zero:
        raise ValueError("a == 0, order_at(%s, %s, %s) == +oo" % (a, p, t))
    if p == Poly(t, t):
        return a.as_poly(t).ET()[0][0]

    n = -1
    p1 = Poly(1, t)
    r = Poly(0, t)
    while r.is_zero:
        n += 1
        p1 = p1*p
        r = a.rem(p1)

    return n

def weak_normalizer(a, d, D, x, T, z=None):
    """
    Weak normalization.

    Given a derivation D on k[t] and f == a/d in k(t), return q in k[t] such
    that f - Dq/q is weakly normalized with respect to t.

    f in k(t) is said to be "weakly normalized" with respect to t if
    residue_p(f) is not a positive integer for any normal irreducible p in k[t]
    such that f is in R_p (Definition 6.1.1).  If f has an elementary integral,
    this is equivalent to no logarithm of integral(f) whose argument depends on
    t has a positive integer coefficient, where the arguments of the logarithms
    not in k(t) are in k[t].

    Returns (q, f - Dq/q)
    """
    t = T[-1]
    z = z or Symbol('z', dummy=True)
    dn, ds = splitfactor(d, D, x, T)

    # Compute d1, where dn == d1*d2**2*...*dn**n is a square-free
    # factorization of d.
    g = gcd(dn, dn.diff(t))
    d_sqf_part = dn.quo(g)
    d1 = d_sqf_part.quo(gcd(d_sqf_part, g))

    a1, b = gcdex_diophantine(d.quo(d1), d1, a)
    r = (a - Poly(z, *T)*derivation(d1, D, x, T)).as_poly(t).resultant(d1.as_poly(t))
    r = Poly(r, z)

    if not r.has(z):
        return (Poly(1, *T), (a, d))

    N = [i for i in r.real_roots() if i in ZZ and i > 0]

    q = reduce(mul, [gcd(a - Poly(n, *T)*derivation(d1, D, x, T), d1) for n in N],
        Poly(1, *T))

    dq = derivation(q, D, x, T)
    sn = q*a - d*dq
    sd = q*d
    sn, sd = sn.cancel(sd, include=True)

    return (q, (sn, sd))

def normal_denominator(fa, fd, ga, gd, D, x, T):
    """
    Normal part of the denominator.

    Given a derivation D on k[t] and f, g in k(t) with f weakly normalized with
    respect to t, either raise NonElementaryIntegral, in which case the equation
    Dy + f*y == g has no solution in k(t), or the quadruplet (a, b, c, h) such
    that a, h in k[t], b, c in k<t>, and for any solution y in k(t) of
    Dy + f*y == g, q = y*h in k<t> satisfies a*Dq + b*q == c.

    This constitutes step 1 in the outline given in the rde.py docstring.
    """
    t = T[-1]
    dn, ds = splitfactor(fd, D, x, T)
    en, es = splitfactor(gd, D, x, T)

    p = gcd(dn, es)
    h = gcd(en, en.diff(t)).quo(gcd(p, p.diff(t)))

    a = dn*h
    c = a*h
    if c.div(en)[1]:
        # en does not divide dn*h**2
        raise NonElementaryIntegral
    ca = c*ga
    ca, cd = ca.cancel(gd, include=True)

    ba = a*fa - dn*derivation(h, D, x, T)*fd
    ba, bd = ba.cancel(fd, include=True)

    # (dn*h, dn*h*f - dn*Dh, dn*h**2*g, h)
    return (a, (ba, bd), (ca, cd), h)

def special_denom(a, ba, bd, ca, cd, D, x, T, case='auto'):
    """
    Special part of the denominator.

    case is one of {'exp', 'tan', 'primitive'} for the hyperexponential,
    hypertangent, and primitive cases, respectively.  For the hyperexponential
    (resp. hypertangent) case, given a derivation D on k[t] and a in k[t], b, c,
    in k<t> with Dt/t (resp. Dt/(t**2 + 1)) in k (sqrt(-1) not in k), a != 0,
    and gcd(a, t) == 1, return the quadruplet (A, B, C, 1/h) such that A, B, C,
    h in k[t] and for any solution q in k<t> of a*Dq + b*q == c, r = qh in k[t]
    satisfies A*Dr + B*r == C.

    For case == 'primitive', k<t> == k[t], so it returns (a, b, c, 1) in this
    case.

    This constitutes step 2 of the outline given in the rde.py docstring.
    """
    # TODO: finish writing this and write tests
    t = T[-1]
    d = D[-1]

    if case == 'auto':
        case = get_case(d, x, t)

    if case == 'exp':
        p = Poly(t, *T)
    elif case == 'tan':
        p = Poly(t**2 + 1, *T)
    elif case in ['primitive', 'base']:
        B = ba.quo(bd)
        C = cd.quo(cd)
        return (a, B, C, Poly(1, *T))
    else:
        raise ValueError("case must be one of {'exp', 'tan', 'primitive', " +
            "'base'}, not %s." % case)

    nb = order_at(ba, p, t) - order_at(bd, p, t)
    nc = order_at(ca, p, t) - order_at(cd, p, t)

    n = min(0, nc - min(0, nb))
    if not nb:
        # Possible cancelation.
        #
        # if case == 'exp':
        #     alpha = (-b/a).rem(p) == -b(0)/a(0)
        #     if alpha == m*Dt/t + Dz/z # parametric logarithmic derivative problem
        #         n = min(n, m)
        # elif case == 'tan':
        #     alpha*sqrt(-1) + beta = (-b/a).rem(p) == -b(sqrt(-1))/a(sqrt(-1))
        #     eta = derivation(t, D, x, t).quo(Poly(t**2 + 1, t)) # eta in k
        #     if 2*beta == Dv/v for some v in k* (see pg. 176) and \
        #     alpha*sqrt(-1) + beta == 2*m*eta*sqrt(-1) + Dz/z:
        #     # parametric logarithmic derivative problem
        #         n = min(n, m)
        raise NotImplementedError("The ability to solve the parametric " +
            "logarithmic derivative problem is required to solve this RDE.")

    N = max(0, -nb, n - nc)
    pN = p**N
    pn = p**n

    A = a*pN
    B = (ba*pN.quo(bd) + Poly(n, t)*a*derivation(p, D, x, T)*pN.quo(p))
    C = ca*pN.quo(pn).quo(cd)
    h = pn # This is 1/h

    # (ap**N, (b + n*a*Dp/p)*p**N, c*p**(N - n), p**-n)
    return (A, B, C, h)

def bound_degree(a, b, c, D, x, T, case='auto'):
    """
    Bound on polynomial solutions.

    Given a derivation D on k[t] and a, b, c in k[t] with a != 0, return n in
    ZZ such that deg(q) <= n for any solution q in k[t] of a*Dq + b*q == c.
    """
    # TODO: finish writing this and write tests
    t = T[-1]
    d = D[-1]

    if case == 'auto':
        case = get_case(d, x, t)

    da = a.degree(t)
    db = b.degree(t)
    dc = c.degree(t)

    alpha = -b.as_poly(t).LC().as_basic()/c.as_poly(t).LC().as_basic()
    alpha = cancel(alpha)

    if case == 'base':
        n = max(0, dc - max(db, da - 1))
        if db == da - 1 and alpha.is_Integer:
            n = max(0, alpha, dc - db)

    elif case == 'primitive':
        if db > da:
            n = max(0, dc - db)
        else:
            max(0, dc - da + 1)

        if db == da - 1:
            raise NotImplementedError("Possible cancellation cases are " +
                "not yet implemented for the primitive case.")
            # if alpha == m*Dt + Dz for z in k and m in ZZ:
                # n = max(n, m)

        if db == da:
            raise NotImplementedError("Possible cancellation cases are " +
            "not yet implemented for the primitive case.")
            # if alpha == Dz/z for z in k*:
                # beta = -lc(a*Dz + b*z)/(z*lc(a))
                # if beta == m*Dt + Dw for w in k and m in ZZ:
                    # n = max(n, m)

    elif case == 'exp':
        n = max(0, dc - max(db, da))
        if da == db:
            raise NotImplementedError("Possible cancellation cases are " +
                "not yet implemented hyperexponential case.")
            # if alpha == m*Dt/t + Dz/z for z in k* and m in ZZ:
                # n = max(n, m)

    elif case in ['tan', 'other_nonlinear']:
        delta = d.degree(t)
        lam = d.LC()
        n = max(0, dc - max(da + delta - 1, db))
        if db == da + delta - 1 and alpha.is_Integer:
            n = max(0, alpha, dc - db)

    else:
        raise ValueError("case must be one of {'exp', 'tan', 'primitive', " +
            "'other_nonlinear', 'base'}, not %s." % case)

    return n

def spde(a, b, c, D, n, x, T):
    """
    Rothstein's Special Polynomial Differential Equation algorithm.

    Given a derivation D on k[t], an integer n and a, b, c in k[t] with a != 0,
    either raise NonElementaryIntegral, in whih case the equation
    a*Dq + b*q == c has no solution of degree at most n in k[t], or return the
    tuple (B, C, m, alpha, beta) such that B, C, alpha, beta in k[t], m in ZZ,
    and any solution q in k[t] of degree at most n of a*Dq + b*q == c must be of
    the form q == alpha*h + beta, where h in k[t], deg(h) <= m, and
    Dh + B*h == C.
    """
    # TODO: Rewrite this non-recursively
    t = T[-1]
    zero = Poly(0, *T)
    if n < 0:
        if c.is_zero:
            return (zero, zero, 0, zero, zero)
        raise NonElementaryIntegral

    g = a.gcd(b)
    if not c.rem(g).is_zero: # g does not divide c
        raise NonElementaryIntegral

    a, b, c = a.quo(g), b.quo(g), c.quo(g)
    if a.degree(t) == 0:
        return (b.quo(a), c.quo(a), n, Poly(1, *T), zero)

    r, z = gcdex_diophantine(b.as_poly(t), a.as_poly(t), c.as_poly(t))
    r, z = Poly(r, *T), Poly(z, *T)
    u = (a, b + derivation(a, D, x, T), z - derivation(r, D, x, T), D,
        n - a.degree(t)) + (x, T)
    B, C, m, alpha, beta = spde(*u)

    return (B, C, m, a*alpha, a*beta + r)

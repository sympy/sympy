"""
Algorithms for solving the Risch differential equation.

Given a differential field K of characteristic 0 that is a simple
monomial extension of a base field k and f, g in K, the Risch
Differential Equation problem is to decide if there exist y in K such
that Dy + f*y == g and to find one if there are some.  If t is a
monomial over k and the coefficients of f and g are in k(t), then y is
in k(t), and the outline of the algorithm here is given as:

1. Compute the normal part n of the denominator of y.  The problem is
then reduced to finding y' in k<t>, where y == y'/n.
2. Compute the special part s of the denominator of y.   The problem is
then reduced to finding y'' in k[t], where y == y''/(n*s)
3. Bound the degree of y''.
4. Reduce the equation Dy + f*y == g to a similar equation with f, g in
k[t].
5. Find the solutions in k[t] of bounded degree of the reduced equation.

See Chapter 6 of "Symbolic Integration I: Transcendental Functions" by
Manuel Bronstein.  See also the docstring of risch.py.
"""
from __future__ import annotations

from operator import mul
from functools import reduce

from sympy.core import oo
from sympy.core.symbol import Dummy

from sympy.polys import Poly, gcd, ZZ, cancel

from sympy.integrals.risch import (gcdex_diophantine, frac_in, derivation,
    splitfactor, NonElementaryIntegralException, DecrementLevel, recognize_log_derivative)

def order_at(a, p, t):
    """
    Computes the order of a at p, with respect to t.

    Explanation
    ===========

    For a, p in k[t], the order of a at p is defined as nu_p(a) = max({n
    in Z+ such that p**n|a}), where a != 0.  If a == 0, nu_p(a) = +oo.

    To compute the order at a rational function, a/b, use the fact that
    nu_p(a/b) == nu_p(a) - nu_p(b).

    This is the order function nu_p from Section 4.1 of Bronstein's book.
    """
    if a.is_zero:
        return oo
    if p == Poly(t, t):
        return a.as_poly(t).ET()[0][0]

    # Uses binary search for calculating the power. power_list collects the tuples
    # (p^k,k) where each k is some power of 2. After deciding the largest k
    # such that k is power of 2 and p^k|a the loop iteratively calculates
    # the actual power.
    power_list = []
    p1 = p
    r = a.rem(p1)
    tracks_power = 1
    while r.is_zero:
        power_list.append((p1,tracks_power))
        p1 = p1*p1
        tracks_power *= 2
        r = a.rem(p1)
    n = 0
    product = Poly(1, t)
    while len(power_list) != 0:
        final = power_list.pop()
        productf = product*final[0]
        r = a.rem(productf)
        if r.is_zero:
            n += final[1]
            product = productf
    return n


def order_at_oo(a, d, t):
    """
    Computes the order of a/d at oo (infinity), with respect to t.

    For f in k(t), the order or f at oo is defined as deg(d) - deg(a), where
    f == a/d.

    This is the order at infinity nu_oo from Section 4.3 of Bronstein's
    book.
    """
    if a.is_zero:
        return oo
    return d.degree(t) - a.degree(t)


def weak_normalizer(a, d, DE, z=None):
    """
    Weak normalization.

    Explanation
    ===========

    Given a derivation D on k[t] and f == a/d in k(t), return q in k[t]
    such that f - Dq/q is weakly normalized with respect to t.

    f in k(t) is said to be "weakly normalized" with respect to t if
    residue_p(f) is not a positive integer for any normal irreducible p
    in k[t] such that f is in R_p (Definition 6.1.1).  If f has an
    elementary integral, this is equivalent to no logarithm of
    integral(f) whose argument depends on t has a positive integer
    coefficient, where the arguments of the logarithms not in k(t) are
    in k[t].

    Returns (q, f - Dq/q)

    This is ``WeakNormalizer`` from Section 6.1 of Bronstein's book.
    """
    z = z or Dummy('z')
    dn, ds = splitfactor(d, DE)

    # Compute d1, where dn == d1*d2**2*...*dn**n is a square-free
    # factorization of d.
    g = gcd(dn, dn.diff(DE.t))
    d_sqf_part = dn.quo(g)
    d1 = d_sqf_part.quo(gcd(d_sqf_part, g))

    a1, b = gcdex_diophantine(d.quo(d1).as_poly(DE.t), d1.as_poly(DE.t),
        a.as_poly(DE.t))
    r = (a - Poly(z, DE.t)*derivation(d1, DE)).as_poly(DE.t).resultant(
        d1.as_poly(DE.t))
    r = Poly(r, z)

    if not r.expr.has(z):
        return (Poly(1, DE.t), (a, d))

    N = [i for i in r.real_roots() if i in ZZ and i > 0]

    q = reduce(mul, [gcd(a - Poly(n, DE.t)*derivation(d1, DE), d1) for n in N],
        Poly(1, DE.t))

    dq = derivation(q, DE)
    sn = q*a - d*dq
    sd = q*d
    sn, sd = sn.cancel(sd, include=True)

    return (q, (sn, sd))


def normal_denom(fa, fd, ga, gd, DE):
    """
    Normal part of the denominator.

    Explanation
    ===========

    Given a derivation D on k[t] and f, g in k(t) with f weakly
    normalized with respect to t, either raise NonElementaryIntegralException,
    in which case the equation Dy + f*y == g has no solution in k(t), or the
    quadruplet (a, b, c, h) such that a, h in k[t], b, c in k<t>, and for any
    solution y in k(t) of Dy + f*y == g, q = y*h in k<t> satisfies
    a*Dq + b*q == c.

    This constitutes step 1 in the outline given in the rde.py docstring.

    This is ``RdeNormalDenominator`` from Section 6.1 of Bronstein's book.
    """
    dn, ds = splitfactor(fd, DE)
    en, es = splitfactor(gd, DE)

    p = dn.gcd(en)
    h = en.gcd(en.diff(DE.t)).quo(p.gcd(p.diff(DE.t)))

    a = dn*h
    c = a*h
    if c.div(en)[1]:
        raise NonElementaryIntegralException("en does not divide dn*h**2 in "
            "normal_denom().")
    ca = c*ga
    ca, cd = ca.cancel(gd, include=True)

    ba = a*fa - dn*derivation(h, DE)*fd
    ba, bd = ba.cancel(fd, include=True)

    # (dn*h, dn*h*f - dn*Dh, dn*h**2*g, h)
    return (a, (ba, bd), (ca, cd), h)


def _special_denom_cancel_bound(a, ba, bd, n, DE, case):
    """
    Sharpen the bound n on the special part of the denominator in the
    possible-cancellation case nu_p(b) == 0.

    This is the part of the special denominator algorithms (Sections 6.2 and
    7.1 of Bronstein) shared by special_denom() and prde_special_denom().
    ``case`` is 'exp' or 'tan'.

    This is the possible-cancellation part of
    ``RdeSpecialDenomExp``/``RdeSpecialDenomTan`` (Section 6.2) and
    ``ParamRdeSpecialDenomExp``/``ParamRdeSpecialDenomTan`` (Section 7.1)
    of Bronstein's book.
    """
    # Delayed imports, since prde imports this module at the module level.
    from .prde import parametric_log_deriv, real_imag

    if case == 'exp':
        dcoeff = DE.d.quo(Poly(DE.t, DE.t))
        with DecrementLevel(DE):  # We are guaranteed to not have problems,
                                  # because case != 'base'.
            alphaa, alphad = frac_in(-ba.eval(0)/bd.eval(0)/a.eval(0), DE.t)
            etaa, etad = frac_in(dcoeff, DE.t)
            A = parametric_log_deriv(alphaa, alphad, etaa, etad, DE)
            if A is not None:
                Q, m, z = A
                if Q == 1:
                    n = min(n, m)

    elif case == 'tan':
        dcoeff = DE.d.quo(Poly(DE.t**2 + 1, DE.t))
        # alpha*sqrt(-1) + beta == Remainder(-b/a, t**2 + 1), with alpha and
        # beta in k.  real_imag() computes this without introducing sqrt(-1)
        # into the coefficient domain.  This must be done at the current
        # level, where t is the hypertangent monomial.
        betaa, alphaa, alphad = real_imag(-ba, bd*a, DE.t)
        betad = alphad
        alpha = alphaa.as_expr()/alphad.as_expr()
        beta = betaa.as_expr()/betad.as_expr()

        with DecrementLevel(DE):  # We are guaranteed to not have problems,
                                  # because case != 'base'.
            alphaa, alphad = frac_in(alpha, DE.t)
            betaa, betad = frac_in(beta, DE.t)
            etaa, etad = frac_in(dcoeff, DE.t)

            if recognize_log_derivative(Poly(2, DE.t)*betaa, betad, DE):
                A = parametric_log_deriv(alphaa, alphad, etaa, etad, DE)
                if A is not None:
                    Q, m, z = A
                    # The condition from the book is
                    # alpha*sqrt(-1) + beta == 2*m*eta*sqrt(-1) + Dz/z for
                    # z in k(sqrt(-1))* and m in ZZ, in which case
                    # n = min(n, m).  parametric_log_deriv() solves
                    # n*alpha == Dv/v + m*eta over k, so its m corresponds
                    # to 2*m in the book's condition.  This is only an
                    # approximation of the real condition (z ranges over
                    # k(sqrt(-1))*, not k*); the complete version needs the
                    # structure theorems.
                    if Q == 1 and m % 2 == 0:
                        n = min(n, m//2)

    return n


def special_denom(a, ba, bd, ca, cd, DE, case='auto'):
    """
    Special part of the denominator.

    Explanation
    ===========

    case is one of {'exp', 'tan', 'primitive'} for the hyperexponential,
    hypertangent, and primitive cases, respectively.  For the
    hyperexponential (resp. hypertangent) case, given a derivation D on
    k[t] and a in k[t], b, c, in k<t> with Dt/t in k (resp. Dt/(t**2 + 1) in
    k, sqrt(-1) not in k), a != 0, and gcd(a, t) == 1 (resp.
    gcd(a, t**2 + 1) == 1), return the quadruplet (A, B, C, 1/h) such that
    A, B, C, h in k[t] and for any solution q in k<t> of a*Dq + b*q == c,
    r = qh in k[t] satisfies A*Dr + B*r == C.

    For ``case == 'primitive'``, k<t> == k[t], so it returns (a, b, c, 1) in
    this case.

    This constitutes step 2 of the outline given in the rde.py docstring.

    This is ``RdeSpecialDenomExp`` and ``RdeSpecialDenomTan`` from Section
    6.2 of Bronstein's book.
    """
    # TODO: finish writing this and write tests

    if case == 'auto':
        case = DE.case

    if case == 'exp':
        p = Poly(DE.t, DE.t)
    elif case == 'tan':
        p = Poly(DE.t**2 + 1, DE.t)
    elif case in ('primitive', 'base'):
        B = ba.to_field().quo(bd)
        C = ca.to_field().quo(cd)
        return (a, B, C, Poly(1, DE.t))
    elif case in ('other_linear', 'other_nonlinear'):
        raise NotImplementedError("The %s case is not implemented in "
            "special_denom()." % case)
    else:
        raise ValueError("case must be one of {'exp', 'tan', 'primitive', "
            "'base', 'other_linear', 'other_nonlinear'}, not %s." % case)

    nb = order_at(ba, p, DE.t) - order_at(bd, p, DE.t)
    nc = order_at(ca, p, DE.t) - order_at(cd, p, DE.t)

    n = min(0, nc - min(0, nb))
    if not nb:
        # Possible cancellation.
        n = _special_denom_cancel_bound(a, ba, bd, n, DE, case)

    # Note that this N, from Section 6.2, intentionally differs from the
    # N == max(0, -nb) of the parametric version in prde_special_denom()
    # (Section 7.1).
    N = max(0, -nb, n - nc)
    pN = p**N
    pn = p**-n

    A = a*pN
    B = ba*pN.quo(bd) + Poly(n, DE.t)*a*derivation(p, DE).quo(p)*pN
    C = (ca*pN*pn).quo(cd)
    h = pn

    # (a*p**N, (b + n*a*Dp/p)*p**N, c*p**(N - n), p**-n)
    return (A, B, C, h)


def bound_degree(a, b, cQ, DE, case='auto', parametric=False):
    """
    Bound on polynomial solutions.

    Explanation
    ===========

    Given a derivation D on k[t] and ``a``, ``b``, ``c`` in k[t] with ``a != 0``, return
    n in ZZ such that deg(q) <= n for any solution q in k[t] of
    a*Dq + b*q == c, when parametric=False, or deg(q) <= n for any solution
    c1, ..., cm in Const(k) and q in k[t] of a*Dq + b*q == Sum(ci*gi, (i, 1, m))
    when parametric=True.

    For ``parametric=False``, ``cQ`` is ``c``, a ``Poly``; for ``parametric=True``, ``cQ`` is Q ==
    [q1, ..., qm], a list of Polys.

    This constitutes step 3 of the outline given in the rde.py docstring.

    This combines ``RdeBoundDegreeBase``, ``RdeBoundDegreePrim``,
    ``RdeBoundDegreeExp``, and ``RdeBoundDegreeNonLinear`` from Section 6.3
    of Bronstein's book (and their parametric analogues from Section 7.1).
    """
    # TODO: finish writing this and write tests

    if case == 'auto':
        case = DE.case

    da = a.degree(DE.t)
    db = b.degree(DE.t)

    # The parametric and regular cases are identical, except for this part
    if parametric:
        dc = max(i.degree(DE.t) for i in cQ)
    else:
        dc = cQ.degree(DE.t)

    alpha = cancel(-b.as_poly(DE.t).LC().as_expr()/
        a.as_poly(DE.t).LC().as_expr())

    if case == 'base':
        n = max(0, dc - max(db, da - 1))
        if db == da - 1 and alpha.is_Integer:
            n = max(0, alpha, dc - db)

    elif case == 'primitive':
        if db > da:
            n = max(0, dc - db)
        else:
            n = max(0, dc - da + 1)

        etaa, etad = frac_in(DE.d, DE.T[DE.level - 1])

        t1 = DE.t
        with DecrementLevel(DE):
            alphaa, alphad = frac_in(alpha, DE.t)
            if db == da - 1:
                from .prde import limited_integrate
                # if alpha == m*Dt + Dz for z in k and m in ZZ:
                try:
                    A = limited_integrate(alphaa, alphad, [(etaa, etad)], DE)
                except NonElementaryIntegralException:
                    A = None
                if A is not None:
                    (za, zd), m = A
                    if len(m) != 1:
                        raise ValueError("Length of m should be 1")
                    m = m[0].as_expr()
                    if m.is_Integer:
                        n = max(n, m)

            elif db == da:
                # if alpha == Dz/z for z in k*:
                    # beta = -lc(a*Dz + b*z)/(z*lc(a))
                    # if beta == m*Dt + Dw for w in k and m in ZZ:
                        # n = max(n, m)
                from .prde import is_log_deriv_k_t_radical_in_field
                A = is_log_deriv_k_t_radical_in_field(alphaa, alphad, DE)
                if A is not None:
                    aa, z = A
                    if aa == 1:
                        beta = -(a*derivation(z, DE, basic=True).as_poly(t1) +
                            b*z.as_poly(t1)).LC()/(z.as_expr()*a.LC())
                        betaa, betad = frac_in(beta, DE.t)
                        from .prde import limited_integrate
                        try:
                            A = limited_integrate(betaa, betad,
                                [(etaa, etad)], DE)
                        except NonElementaryIntegralException:
                            A = None
                        if A is not None:
                            (za, zd), m = A
                            if len(m) != 1:
                                raise ValueError("Length of m should be 1")
                            m = m[0].as_expr()
                            if m.is_Integer:
                                n = max(n, m)

    elif case == 'exp':
        from .prde import parametric_log_deriv

        n = max(0, dc - max(db, da))
        if da == db:
            etaa, etad = frac_in(DE.d.quo(Poly(DE.t, DE.t)), DE.T[DE.level - 1])
            with DecrementLevel(DE):
                alphaa, alphad = frac_in(alpha, DE.t)
                A = parametric_log_deriv(alphaa, alphad, etaa, etad, DE)
                if A is not None:
                    # if alpha == m*Dt/t + Dz/z for z in k* and m in ZZ:
                        # n = max(n, m)
                    a, m, z = A
                    if a == 1:
                        n = max(n, m)

    elif case in ('tan', 'other_nonlinear'):
        delta = DE.d.degree(DE.t)
        lam = DE.d.LC()
        alpha = cancel(alpha/lam)
        n = max(0, dc - max(da + delta - 1, db))
        if db == da + delta - 1 and alpha.is_Integer:
            n = max(0, alpha, dc - db)

    elif case == 'other_linear':
        raise NotImplementedError("The other_linear case is not implemented "
            "in bound_degree().")
    else:
        raise ValueError("case must be one of {'exp', 'tan', 'primitive', "
            "'base', 'other_linear', 'other_nonlinear'}, not %s." % case)

    return n


def spde(a, b, c, n, DE):
    """
    Rothstein's Special Polynomial Differential Equation algorithm.

    Explanation
    ===========

    Given a derivation D on k[t], an integer n and ``a``,``b``,``c`` in k[t] with
    ``a != 0``, either raise NonElementaryIntegralException, in which case the
    equation a*Dq + b*q == c has no solution of degree at most ``n`` in
    k[t], or return the tuple (B, C, m, alpha, beta) such that B, C,
    alpha, beta in k[t], m in ZZ, and any solution q in k[t] of degree
    at most n of a*Dq + b*q == c must be of the form
    q == alpha*h + beta, where h in k[t], deg(h) <= m, and Dh + B*h == C.

    This constitutes step 4 of the outline given in the rde.py docstring.

    This is ``SPDE`` from Section 6.4 of Bronstein's book.
    """
    zero = Poly(0, DE.t)

    alpha = Poly(1, DE.t)
    beta = Poly(0, DE.t)

    # Generous pass allowance for n == oo (see below); any solution of
    # degree d is found within about d/deg(a) passes, so this only
    # defers solutions of very high degree to an honest error.
    oo_passes = 4*(a.degree(DE.t) + b.degree(DE.t) + c.degree(DE.t) + 4)

    while True:
        if c.is_zero:
            return (zero, zero, 0, zero, beta)  # -1 is more to the point
        if (n < 0) is True:
            raise NonElementaryIntegralException("n became negative in "
                "spde().")

        g = a.gcd(b)
        if not c.rem(g).is_zero:  # g does not divide c
            raise NonElementaryIntegralException("gcd(a, b) does not divide "
                "c in spde().")

        a, b, c = a.quo(g), b.quo(g), c.quo(g)

        if a.degree(DE.t) == 0:
            b = b.to_field().quo(a)
            c = c.to_field().quo(a)
            return (b, c, n, alpha, beta)

        if n is oo:
            # Each pass through the reduction below decreases n by
            # deg(a) > 0 and "no solution" is detected by n < 0, so with
            # no degree bound the loop can run forever when there is no
            # solution (e.g. a == t, b == 1, c == x over t == exp(x):
            # gcd(a, b) == 1 stays trivial, c cycles through nonzero
            # constants, and no other exit is ever taken).  Allow a
            # generous number of passes (solvable cases terminate
            # quickly, with c becoming zero), then give up with an
            # honest error -- NotImplementedError, not a nonexistence
            # claim.
            oo_passes -= 1
            if oo_passes < 0:
                raise NotImplementedError("spde() ran out of passes with "
                    "an infinite degree bound; cannot decide solvability.")

        r, z = gcdex_diophantine(b, a, c)
        b += derivation(a, DE)
        c = z - derivation(r, DE)
        n -= a.degree(DE.t)

        beta += alpha * r
        alpha *= a

def no_cancel_b_large(b, c, n, DE):
    """
    Poly Risch Differential Equation - No cancellation: deg(b) large enough.

    Explanation
    ===========

    Given a derivation D on k[t], ``n`` either an integer or +oo, and ``b``,``c``
    in k[t] with ``b != 0`` and either D == d/dt or
    deg(b) > max(0, deg(D) - 1), either raise NonElementaryIntegralException, in
    which case the equation ``Dq + b*q == c`` has no solution of degree at
    most n in k[t], or a solution q in k[t] of this equation with
    ``deg(q) < n``.

    This is ``PolyRischDENoCancel1`` from Section 6.5 of Bronstein's book.
    """
    q = Poly(0, DE.t)

    while not c.is_zero:
        m = c.degree(DE.t) - b.degree(DE.t)
        if not 0 <= m <= n:  # n < 0 or m < 0 or m > n
            raise NonElementaryIntegralException("deg(c) - deg(b) is not "
                "between 0 and n in no_cancel_b_large().")

        p = Poly(c.as_poly(DE.t).LC()/b.as_poly(DE.t).LC()*DE.t**m, DE.t,
            expand=False)
        q = q + p
        n = m - 1
        c = c - derivation(p, DE) - b*p

    return q


def no_cancel_b_small(b, c, n, DE):
    """
    Poly Risch Differential Equation - No cancellation: deg(b) small enough.

    Explanation
    ===========

    Given a derivation D on k[t], ``n`` either an integer or +oo, and ``b``,``c``
    in k[t] with deg(b) < deg(D) - 1 and either D == d/dt or
    deg(D) >= 2, either raise NonElementaryIntegralException, in which case the
    equation Dq + b*q == c has no solution of degree at most n in k[t],
    or a solution q in k[t] of this equation with deg(q) <= n, or the
    tuple (h, b0, c0) such that h in k[t], b0, c0, in k, and for any
    solution q in k[t] of degree at most n of Dq + bq == c, y == q - h
    is a solution in k of Dy + b0*y == c0.

    This is ``PolyRischDENoCancel2`` from Section 6.5 of Bronstein's book.
    """
    q = Poly(0, DE.t)

    while not c.is_zero:
        if n == 0:
            m = 0
        else:
            m = c.degree(DE.t) - DE.d.degree(DE.t) + 1

        if not 0 <= m <= n:  # n < 0 or m < 0 or m > n
            raise NonElementaryIntegralException("m is not between 0 and n "
                "in no_cancel_b_small().")

        if m > 0:
            p = Poly(c.as_poly(DE.t).LC()/(m*DE.d.as_poly(DE.t).LC())*DE.t**m,
                DE.t, expand=False)
        else:
            if b.degree(DE.t) != c.degree(DE.t):
                raise NonElementaryIntegralException("deg(b) != deg(c) in "
                    "the m == 0 case of no_cancel_b_small().")
            if b.degree(DE.t) == 0:
                return (q, b.as_poly(DE.T[DE.level - 1]),
                    c.as_poly(DE.T[DE.level - 1]))
            p = Poly(c.as_poly(DE.t).LC()/b.as_poly(DE.t).LC(), DE.t,
                expand=False)

        q = q + p
        n = m - 1
        c = c - derivation(p, DE) - b*p

    return q


# TODO: better name for this function
def no_cancel_equal(b, c, n, DE):
    """
    Poly Risch Differential Equation - No cancellation: deg(b) == deg(D) - 1

    Explanation
    ===========

    Given a derivation D on k[t] with deg(D) >= 2, n either an integer
    or +oo, and b, c in k[t] with deg(b) == deg(D) - 1, either raise
    NonElementaryIntegralException, in which case the equation Dq + b*q == c has
    no solution of degree at most n in k[t], or a solution q in k[t] of
    this equation with deg(q) <= n, or the tuple (h, m, C) such that h
    in k[t], m in ZZ, and C in k[t], and for any solution q in k[t] of
    degree at most n of Dq + b*q == c, y == q - h is a solution in k[t]
    of degree at most m of Dy + b*y == C.

    This is ``PolyRischDENoCancel3`` from Section 6.5 of Bronstein's book.
    """
    q = Poly(0, DE.t)
    lc = cancel(-b.as_poly(DE.t).LC()/DE.d.as_poly(DE.t).LC())
    if lc.is_Integer and lc.is_positive:
        M = lc
    else:
        M = -1

    while not c.is_zero:
        m = max(M, c.degree(DE.t) - DE.d.degree(DE.t) + 1)

        if not 0 <= m <= n:  # n < 0 or m < 0 or m > n
            raise NonElementaryIntegralException("m is not between 0 and n "
                "in no_cancel_equal().")

        u = cancel(m*DE.d.as_poly(DE.t).LC() + b.as_poly(DE.t).LC())
        if u.is_zero:
            return (q, m, c)
        if m > 0:
            p = Poly(c.as_poly(DE.t).LC()/u*DE.t**m, DE.t, expand=False)
        else:
            if c.degree(DE.t) != DE.d.degree(DE.t) - 1:
                raise NonElementaryIntegralException("deg(c) != deg(D) - 1 "
                    "in the m == 0 case of no_cancel_equal().")
            else:
                p = Poly(c.as_poly(DE.t).LC()/b.as_poly(DE.t).LC(), DE.t,
                    expand=False)

        q = q + p
        n = m - 1
        c = c - derivation(p, DE) - b*p

    return q


def cancel_primitive(b, c, n, DE):
    """
    Poly Risch Differential Equation - Cancellation: Primitive case.

    Explanation
    ===========

    Given a derivation D on k[t], n either an integer or +oo, ``b`` in k, and
    ``c`` in k[t] with Dt in k and ``b != 0``, either raise
    NonElementaryIntegralException, in which case the equation Dq + b*q == c
    has no solution of degree at most n in k[t], or a solution q in k[t] of
    this equation with deg(q) <= n.

    This is ``PolyRischDECancelPrim`` from Section 6.6 of Bronstein's book.
    """
    # Delayed imports
    from .prde import is_log_deriv_k_t_radical_in_field, is_deriv_in_field
    with DecrementLevel(DE):
        ba, bd = frac_in(b, DE.t)
        A = is_log_deriv_k_t_radical_in_field(ba, bd, DE)

    if A is not None:
        m, z = A
        if m == 1:  # b == Dz/z
            # D(z*q) == z*(Dq + b*q), so q solves Dq + b*q == c iff
            # p == z*q satisfies Dp == z*c.  So find an antiderivative
            # p in k[t] of z*c with deg(p) <= n, and return q == p/z.
            zca, zcd = frac_in(z*c.as_expr(), DE.t, cancel=True)
            B = is_deriv_in_field(zca, zcd, DE)
            if B is None:
                raise NonElementaryIntegralException("z*c is not the "
                    "derivative of an element of k(t) in cancel_primitive().")
            pa, pd = B
            # For a primitive t, an antiderivative in k(t) of an element
            # of k[t] is itself in k[t] (it can have no poles), and any
            # additive constant keeps q == p/z in k[t], so any choice of
            # p will do.
            q = cancel(pa.as_expr()/(pd.as_expr()*z)).as_poly(DE.t)
            if q is None:
                raise NonElementaryIntegralException("p/z is not in k[t] "
                    "in cancel_primitive().")
            if q.degree(DE.t) > n:
                raise NonElementaryIntegralException("deg(p/z) > n in "
                    "cancel_primitive().")
            return q

    if c.is_zero:
        return c  # return 0

    if n < c.degree(DE.t):
        raise NonElementaryIntegralException("n < deg(c) in "
            "cancel_primitive().")

    q = Poly(0, DE.t)
    while not c.is_zero:
        m = c.degree(DE.t)
        if n < m:
            raise NonElementaryIntegralException("n < deg(c) in "
                "cancel_primitive().")
        with DecrementLevel(DE):
            a2a, a2d = frac_in(c.LC(), DE.t)
            sa, sd = rischDE(ba, bd, a2a, a2d, DE)
        stm = Poly(sa.as_expr()/sd.as_expr()*DE.t**m, DE.t, expand=False)
        q += stm
        n = m - 1
        c -= b*stm + derivation(stm, DE)

    return q


def cancel_exp(b, c, n, DE):
    """
    Poly Risch Differential Equation - Cancellation: Hyperexponential case.

    Explanation
    ===========

    Given a derivation D on k[t], n either an integer or +oo, ``b`` in k, and
    ``c`` in k[t] with Dt/t in k and ``b != 0``, either raise
    NonElementaryIntegralException, in which case the equation Dq + b*q == c
    has no solution of degree at most n in k[t], or a solution q in k[t] of
    this equation with deg(q) <= n.

    This is ``PolyRischDECancelExp`` from Section 6.6 of Bronstein's book.
    """
    from .prde import parametric_log_deriv, is_deriv_in_field
    eta = DE.d.quo(Poly(DE.t, DE.t)).as_expr()

    with DecrementLevel(DE):
        etaa, etad = frac_in(eta, DE.t)
        ba, bd = frac_in(b, DE.t)
        A = parametric_log_deriv(ba, bd, etaa, etad, DE)

    if A is not None:
        a, m, z = A
        if a == 1:
            # b == Dz/z + m*Dt/t, so D(z*t**m*q) == z*t**m*(Dq + b*q).
            # Hence q solves Dq + b*q == c iff p == z*t**m*q satisfies
            # Dp == z*t**m*c.  So find an antiderivative p in k<t> of
            # z*t**m*c; then q == p/(z*t**m) must be in k[t] with
            # deg(q) <= n.
            pa, pd = frac_in(z*DE.t**m*c.as_expr(), DE.t, cancel=True)
            B = is_deriv_in_field(pa, pd, DE)
            if B is None:
                raise NonElementaryIntegralException("z*t**m*c is not the "
                    "derivative of an element of k(t) in cancel_exp().")
            va, vd = B
            p = cancel(va.as_expr()/vd.as_expr())
            if m > 0:
                # p is determined only up to an additive constant.  For
                # m > 0 the true p == z*t**m*q has no t-free term (all its
                # monomials have degree >= m), while a constant term would
                # make p/(z*t**m) non-polynomial, so remove it.
                p -= p.subs(DE.t, 0)
            # For m <= 0 any additive constant C leaves q in k[t]
            # (C/(z*t**m) == C*t**(-m)/z), and q + C*t**(-m)/z is still a
            # solution, since t**(-m)/z solves the homogeneous equation.
            q = cancel(p/(z*DE.t**m)).as_poly(DE.t)
            if q is None:
                raise NonElementaryIntegralException("p/(z*t**m) is not in "
                    "k[t] in cancel_exp().")
            if q.degree(DE.t) > n:
                if m < 0 and q.degree(DE.t) == -m and \
                        not cancel(z*q.nth(-m)).has(*DE.T):
                    # The homogeneous solutions are C*t**(-m)/z for
                    # constants C, so the t**(-m) coefficient of q is
                    # adjustable by the constant of integration exactly
                    # when it is C/z for a constant C; in that case,
                    # remove that term.  (A non-constant coefficient
                    # cannot be removed by any choice of the constant,
                    # and stripping it anyway would return a q that does
                    # not solve the equation.)
                    q = q - Poly(q.nth(-m)*DE.t**(-m), DE.t)
                if q.degree(DE.t) > n:
                    raise NonElementaryIntegralException("deg(p/(z*t**m)) "
                        "> n in cancel_exp().")
            return q

    if c.is_zero:
        return c  # return 0

    if n < c.degree(DE.t):
        raise NonElementaryIntegralException("n < deg(c) in cancel_exp().")

    q = Poly(0, DE.t)
    while not c.is_zero:
        m = c.degree(DE.t)
        if n < m:
            raise NonElementaryIntegralException("n < deg(c) in "
                "cancel_exp().")
        # a1 = b + m*Dt/t
        a1 = b.as_expr()
        with DecrementLevel(DE):
            # TODO: Write a dummy function that does this idiom
            a1a, a1d = frac_in(a1, DE.t)
            a1a = a1a*etad + etaa*a1d*Poly(m, DE.t)
            a1d = a1d*etad

            a2a, a2d = frac_in(c.LC(), DE.t)

            sa, sd = rischDE(a1a, a1d, a2a, a2d, DE)
        stm = Poly(sa.as_expr()/sd.as_expr()*DE.t**m, DE.t, expand=False)
        q += stm
        n = m - 1
        c -= b*stm + derivation(stm, DE)  # deg(c) becomes smaller
    return q


def solve_poly_rde(b, c, n, DE):
    """
    Solve a Polynomial Risch Differential Equation with degree bound ``n``.

    This constitutes step 4 of the outline given in the rde.py docstring.

    This dispatches among the algorithms of Sections 6.5 and 6.6 of
    Bronstein's book; it has no single named counterpart.  The parametric
    analogue of this dispatch is param_poly_rischDE() in prde.py.
    """
    # No cancellation
    if not b.is_zero and (DE.case == 'base' or
            b.degree(DE.t) > max(0, DE.d.degree(DE.t) - 1)):

        return no_cancel_b_large(b, c, n, DE)

    elif (b.is_zero or b.degree(DE.t) < DE.d.degree(DE.t) - 1) and \
            (DE.case == 'base' or DE.d.degree(DE.t) >= 2):

        R = no_cancel_b_small(b, c, n, DE)

        if isinstance(R, Poly):
            return R
        else:
            # XXX: Might k be a field? (pg. 209)
            h, b0, c0 = R
            with DecrementLevel(DE):
                b0, c0 = b0.as_poly(DE.t), c0.as_poly(DE.t)
                if b0 is None:  # See above comment
                    raise ValueError("b0 should be a non-Null value")
                if c0 is  None:
                    raise ValueError("c0 should be a non-Null value")
                y = solve_poly_rde(b0, c0, n, DE).as_poly(DE.t)
            return h + y

    elif DE.d.degree(DE.t) >= 2 and b.degree(DE.t) == DE.d.degree(DE.t) - 1 and \
            n > -b.as_poly(DE.t).LC()/DE.d.as_poly(DE.t).LC():

        # TODO: Is this check necessary, and if so, what should it do if it fails?
        # b comes from the first element returned from spde()
        if not b.as_poly(DE.t).LC().is_number:
            raise TypeError("Result should be a number")

        R = no_cancel_equal(b, c, n, DE)

        if isinstance(R, Poly):
            return R
        else:
            h, m, C = R
            # XXX: Or should it be rischDE()?
            y = solve_poly_rde(b, C, m, DE)
            return h + y

    else:
        # Cancellation
        if b.is_zero:
            # Dq == c: in-field integration (Section 6.6)
            from .prde import is_deriv_in_field
            B = is_deriv_in_field(c, Poly(1, DE.t), DE)
            if B is None:
                raise NonElementaryIntegralException("c is not the "
                    "derivative of an element of k(t) in solve_poly_rde().")
            va, vd = B
            # For a primitive or hyperexponential t, an antiderivative in
            # k(t) of an element of k[t] is itself in k[t], and any
            # additive constant is also in k[t], so any choice of the
            # antiderivative will do.
            q = cancel(va.as_expr()/vd.as_expr()).as_poly(DE.t)
            if q is None:
                raise NonElementaryIntegralException("The antiderivative "
                    "of c is not in k[t] in solve_poly_rde().")
            if q.degree(DE.t) > n:
                raise NonElementaryIntegralException("The antiderivative "
                    "of c has degree > n in solve_poly_rde().")
            return q
        else:
            if DE.case == 'exp':
                return cancel_exp(b, c, n, DE)

            elif DE.case == 'primitive':
                return cancel_primitive(b, c, n, DE)

            else:
                raise NotImplementedError("Other Poly (P)RDE cancellation "
                    "cases are not yet implemented (%s)." % DE.case)


def rischDE(fa, fd, ga, gd, DE):
    """
    Solve a Risch Differential Equation: Dy + f*y == g.

    Explanation
    ===========

    See the outline in the docstring of rde.py for more information
    about the procedure used.  Either raise NonElementaryIntegralException, in
    which case there is no solution y in the given differential field,
    or return y in k(t) satisfying Dy + f*y == g, or raise
    NotImplementedError, in which case, the algorithms necessary to
    solve the given Risch Differential Equation have not yet been
    implemented.

    This chains together the algorithms of Sections 6.1-6.6 of Bronstein's
    book; the book does not give it as a single named pseudocode function.
    """
    # Substitute y == z/q, where q is the weak normalizer of f.  By Theorem
    # 6.1.2, y in k(t) solves Dy + f*y == g iff z == q*y solves
    # Dz + (f - Dq/q)*z == q*g, and f - Dq/q is weakly normalized, as
    # required by normal_denom().
    q, (fa, fd) = weak_normalizer(fa, fd, DE)
    ga, gd = (ga*q).cancel(gd, include=True)
    a, (ba, bd), (ca, cd), hn = normal_denom(fa, fd, ga, gd, DE)
    A, B, C, hs = special_denom(a, ba, bd, ca, cd, DE)
    try:
        n = bound_degree(A, B, C, DE)
    except NotImplementedError:
        # bound_degree() could not decide one of its structure-theorem
        # queries.  Proceeding with n == oo is sound: no bound is needed
        # to accept a solution, and the no-cancellation and cancellation
        # algorithms all terminate by descending on deg(c).  The only
        # step that needs a finite bound to terminate is the deg(a) > 0
        # reduction in spde(), which raises NotImplementedError instead
        # of looping forever when it would need one.
        n = oo

    B, C, m, alpha, beta = spde(A, B, C, n, DE)
    if C.is_zero:
        y = C
    else:
        y = solve_poly_rde(B, C, m, DE)

    # The solution found so far is z == (alpha*y + beta)/(hn*hs); undo the
    # weak normalizer substitution y == z/q.
    return (alpha*y + beta, q*hn*hs)

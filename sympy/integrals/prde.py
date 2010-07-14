"""
Algorithms for solving Parametric Risch Differential Equations.

The methods used for solving Parametric Risch Differential Equations parallel
those for solving Risch Differential Equations.  See the outline in the
docstring of rde.py for more information.
"""
from sympy.core import Symbol

from sympy.solvers import solve

from sympy.polys import Poly, PolynomialError, lcm, cancel

from sympy.integrals.risch import (derivation, get_case, NonElementaryIntegral,
    residue_reduce, splitfactor, residue_reduce_derivation)
from sympy.integrals.rde import order_at, weak_normalizer

#    from pudb import set_trace; set_trace() # Debugging

def prde_normal_denom(fa, fd, G, D, T):
    """
    Parametric Risch Differential Equation - Normal part of the denominator.

    Given a derivation D on k[t] and f, g1, ..., gm in k(t) with f weakly
    normalized with respect to t, return the tuple (a, b, G, h) such that
    a, h in k[t], b in k<t>, G = [g1, ..., gm] in k(t)^m, and for any solution
    c1, ..., cm in Const(k) and y in k(t) of Dy + f*y == Sum(ci*gi, (i, 1, m)),
    q == y*h in k<t> satisfies a*Dq + b*q == Sum(ci*Gi, (i, 1, m)).
    """
    t = T[-1]
    dn, ds = splitfactor(fd, D, T)
    gd = reduce(lambda i, j: i.lcm(j), zip(*G)[1])
    en, es = splitfactor(gd, D, T)

    p = dn.gcd(en)
    h = en.gcd(en.diff(t)).quo(p.gcd(p.diff(t)))

    a = dn*h
    c = a*h

    ba = a*fa - dn*derivation(h, D, T)*fd
    ba, bd = ba.cancel(fd, include=True)

    G = [(c*A).cancel(D, include=True) for A, D in G]

    return (a, (ba, bd), G, h)

def prde_special_denom(a, ba, bd, G, D, T, case='auto'):
    """
    Parametric Risch Differential Equation - Special part of the denominator.

    case is on of {'exp', 'tan', 'primitive'} for the hyperexponential,
    hypertangent, and primitive cases, respectively.  For the hyperexponential
    (resp. hypertangent) case, given a derivation D on k[t] and a in k[t],
    b in k<t>, and g1, ..., gm in k(t) with Dt/t in k (resp. Dt/(t**2 + 1) in
    k, sqrt(-1) not in k), a != 0, and gcd(a, t) == 1 (resp.
    gcd(a, t**2 + 1) == 1), return the tuple (A, B, GG, h) such that A, B, h in
    k[t], GG = [gg1, ..., ggm] in k(t)^m, and for any solution c1, ..., cm in
    Const(k) and q in k<t> of a*Dq + b*q == Sum(ci*gi, (i, 1, m)), r == q*h in
    k[t] satisfies A*Dr + B*r == Sum(ci*ggi, (i, 1, m)).

    For case == 'primitive', k<t> == k[t], so it returns (a, b, G, 1) in this
    case.
    """
    t = T[-1]
    d = D[-1]

    if case == 'auto':
        case = get_case(d, t)

    if case == 'exp':
        p = Poly(t, t)
    elif case == 'tan':
        p = Poly(t**2 + 1, t)
    elif case in ['primitive', 'base']:
        B = ba.quo(bd)
        return (a, B, G, Poly(1, t))
    else:
        raise ValueError("case must be one of {'exp', 'tan', 'primitive', " +
            "'base'}, not %s." % case)

    nb = order_at(ba, p, t) - order_at(bd, p, t)
    nc = min([order_at(Ga, p, t) - order_at(Gd, p, t) for Ga, Gd in G])

    n = min(0, nc - min(0, nb))
    if not nb:
        # Possible cancelation
        #
        # if case == 'exp':
        #     alpha = (-b/a).rem(p) == -b(0)/a(0)
        #     if alpha == m*Dt/t + Dz/z # parametric logarithmic derivative problem
        #         n = min(n, m)
        # elif case == 'tan':
        #     alpha*sqrt(-1) + beta = (-b/a)/rem(p) == -b(sqrt(-1))/a(sqrt(-1))
        #     eta = derivation(t, D, T).quo(Poly(t**2 + 1, t)) # eta in k
        #     if 2*beta == Db/b for some v in k* (see pg. 176) and \
        #     alpha*sqrt(-1) + beta == 2*b*eta*sqrt(-1) + Dz/z:
        #     # parametric logarithmic derivative problem
        #         n = min(n, m)
        raise NotImplementedError("The ability to solve the parametric " +
            "logarithmic derivative problem is required to solve this PRDE.")

    N = max(0, -nb)
    pN = p**N
    pn = p**-n # This is 1/h

    A = a*pN
    B = ba*pN.quo(bd) + Poly(n, t)*a*derivation(p, D, T).quo(p)*pN
    G = [(Ga*pN*pn).quo(Gd) for Ga, Gd in G]
    h = pn

    # (a*p**N, (b + n*a*Dp/p)*p**N, g1*p**(N - n), ..., gm*p**(N - n), p**-n)
    return (A, B, G, h)

def param_rischDE(fa, fd, G, D, T):
    """
    Solve a Parametric Risch Differential Equation: Dy + f*y == Sum(ci*Gi, (i, 1, m)).
    """
    _, (fa, fd) = weak_normalizer(fa, fd, D, T)
    a, (ba, bd), G, hn = prde_normal_denom(a, ga, gd, G, D, T)
    A, B, G, hs = prde_special_denom(a, ba, bd, G, D, T)


def parametric_log_deriv_heu(fa, fd, wa, wd, D, T):
    """
    Parametric logarithmic derivative heuristic.

    Given a derivation D on k[t], f in k(t), and a hyperexponential monomial
    theta over k(t), raises either NotImplementedError, in which case the
    heuristic failed, or NonElementaryIntegral, in which case it has proven
    that no solution exists, or returns a solution (n, m, v) of the equation
    n*f == Dv/v + m*Dtheta/theta, with v in k(t)* and n, m in ZZ with n != 0.

    If this heuristic fails, the structure theorem approach will need to be
    used.

    The argument w == Dtheta/theta
    """
    # TODO: finish writing this and write tests
    from pudb import set_trace; set_trace() # Debugging
    t = T[-1]
    c1 = Symbol('c', dummy=True)

    p, a = fa.div(fd)
    q, b = wa.div(wd)

    B = max(0, derivation(t, D, T).degree(t) - 1)
    C = max(p.degree(t), q.degree(t))

    if q.degree(t) > B:
        eqs = [p.nth(i) - c1*q.nth(i) for i in range(B + 1, C + 1)]
        s = solve(eqs, c1)
        if not s or not s[c1].is_Rational:
            raise NonElementaryIntegral("parametric_log_deriv_heu(): deg(q) " +
            "> B, no solution for c.")

        N, M = s[c1].as_numer_denom() # N, M are integers
        N, M = Poly(N, t), Poly(M, t)

        nfmwa = N*fa*wd - M*wa*fd
        nfmwd = fd*wd
        Qv = is_log_deriv_k_t_radical(N*fa*wd - M*wa*fd, fd*wd, D, T, 'auto')
        if Qv is None:
            raise NonElementaryIntegral("parametric_log_deriv_heu(): %s/%s " +
                "(N*f - M*w) is not the logarithmic derivaitive of a k(t)-radical."
                % (nfmwa, nfmwd))
        Q, v = Qv

        if Q.is_zero or v.is_zero:
            raise NonElementaryIntegral("parametric_log_deriv_heu(): Q == 0 " +
                "or v == 0.")

        return (Q*N, Q*M, v)

    if p.degree(t) > B:
        raise NonElementaryIntegral("parametric_log_deriv_heu(): p.degree() " +
            "B.")

    c = lcm(fd.as_poly(t).LC(),wd.as_poly(t).LC())
    l = fd.monic().lcm(wd.monic())*Poly(c, t)
    ln, ls = splitfactor(l, D, T)
    z = ls*ln.gcd(ln.diff(t))

    if not z.has(t):
        raise NotImplementedError("parametric_log_deriv_heu() " +
            "heuristic failed: z in k.")

    u1, r1 = (fa*l.quo(fd)).div(z) # (l*f).div(z)
    u2, r2 = (wa*l.quo(wd)).div(z) # (l*w).div(z)

    eqs = [r1.nth(i) - c1*r2.nth(i) for i in range(z.degree(t))]
    s = solve(eqs, c1)
    if not s or not s[c1].is_Rational:
        raise NonElementaryIntegral("parametric_log_deriv_heu(): deg(q) " +
            "<= B, no solution for c.")

    M, N = s[c1].as_numer_denom()
    M, N = Poly(M, t), Poly(N, t)

    nfmwa = N*fa*wd - M*wa*fd
    nfmwd = fd*wd
    Qv = is_log_deriv_k_t_radical(nfmwa, nfmwd, D, T, 'auto')
    if Qv is None:
        raise NonElementaryIntegral("parametric_log_deriv_heu(): %s/%s " +
            "(N*f - M*w) is not the logarithmic derivaitive of a k(t)-radical."
            % (nfmwa, nfmwd))
    Q, v = Qv

    if Q.is_zero or v.is_zero:
        raise NonElementaryIntegral("parametric_log_deriv_heu(): Q == 0 " +
            "or v == 0.")

    return (Q*N, Q*M, v)

def parametric_log_deriv(fa, fd, wa, wd, D, T):
    # TODO: Write the full algorithm using the structure theorems.
    return parametric_log_deriv_heu(fa, fd, wa, wd, D, T)

def is_log_deriv_k_t_radical(fa, fd, D, T, case='auto'):
    """
    Checks if f can be written as the logarithmic derivative of a k(t)-radical.

    f in k(t) can be written as the logarithmic derivative of a k(t) radical if
    there exist n in ZZ and u in k(t) with n, u != 0 such that n*f == Du/u.
    Either returns (n, u) or None, which means that f cannot be written as the
    logarithmic derivative of a k(t)-radical.

    case is one of {'primitive', 'exp', 'tan', 'auto'} for the primitive,
    hyperexponential, and hypertangent cases, respectively.  If case it 'auto',
    it will attempt to determine the type of the derivation automatically.
    """
    # TODO: finish writing this and write tests
    from pudb import set_trace; set_trace() # Debugging
    fa, fd = fa.cancel(fd, include=True)

    # TODO: Fix for x in T
    if not T:
        # Base case.
        # These had better be True.
        assert case in ['auto', 'primitive', 'base']
        assert not D
        t = x
        d = Poly(1, x)
        case = 'base'
    else:
        t = T[-1]
        d = D[-1]

    if case == 'auto':
        case = get_case(d, t)

    # f must be simple
    n, s = splitfactor(fd, D, T)
    if not s.is_one:
        return None

    z = Symbol('z', dummy=True)
    H, b = residue_reduce(fa, fd, D, T, z=z)
    if not b:
        # Note, according to the note on page 255 of Bronstein's book, this
        # should never happen with the parametric logarithmic derivative
        # problems.  I don't know, however, if it will happen when building up
        # the differential field for an integrand.
        raise NotImplementedError("f has a non-elementary integral, cannot " +
            "determine if it is the logarithmic derivative of a k(t)-radical.")

    p = cancel(fa.as_basic()/fd.as_basic() - residue_reduce_derivation(H, D, T, z))
    try:
        p = Poly(p, t)
    except PolynomialError:
        # f - Dg will be in k[t] if f is the logarithmic derivaitve of a k(t)-radical
        return None

    if p.degree(t) >= max(1, d.degree(t)):
        return None

    if case == 'exp':
        wa, wd = derivation(t, D, T).cancel(Poly(t, t), include=True)
        try:
            n, e, u = parametric_log_deriv(p, Poly(1, t), wa, wd, D, T)
        except NonElementaryIntegral:
            return None
        u **= e
        raise NotImplementedError("The hyperexponential case is " +
        "not yet completely implemented for is_log_deriv_k_t_radical().")

    elif case == 'primitive':
        n, u = is_log_deriv_k_t_radical(fa, fd, D[:-1], T[:-1], case='auto')
        raise NotImplementedError("The primitive case is " +
        "not yet completely implemented for is_log_deriv_k_t_radical()")

    elif case == 'base':
        raise NotImplementedError("The base case is " +
        "not yet implemented for is_log_deriv_k_t_radical().")

    elif case == 'tan':
        raise NotImplementedError("The hypertangent case is " +
        "not yet implemented for is_log_deriv_k_t_radical()")

    elif case in ['other_linear', 'other_nonlinear']:
        # XXX: If these are supported by the structure theorems, change to NotImplementedError.
        raise ValueError("The %s case is not supported in this function." % case)

    else:
        raise ValueError("case must be one of {'primitive', 'exp', 'tan', "+
        "'base', 'auto'}, not %s""" % case)

    return (n, u)

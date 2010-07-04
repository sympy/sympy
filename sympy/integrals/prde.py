"""
Algorithms for solving Parametric Risch Differential Equations.
"""
from sympy.core import Symbol

from sympy.solvers import solve

from sympy.polys import Poly, PolynomialError, lcm, cancel

from sympy.integrals.risch import (derivation, get_case, NonElementaryIntegral,
    residue_reduce, splitfactor, residue_reduce_derivation)

def parametric_log_deriv_heu(fa, fd, wa, wd, D, x, t):
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
    from pudb import set_trace; set_trace() # Debugging
    c1 = Symbol('c', dummy=True)

    p, a = fa.div(fd)
    q, b = wa.div(wd)

    B = max(0, derivation(t, D, x, t).degree(t) - 1)
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
        Qv = is_log_deriv_k_t_radical(N*fa*wd - M*wa*fd, fd*wd, D, x, t, 'auto')
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

    c = lcm(fd.LC(),wd.LC())
    l = fd.monic().lcm(wd.monic())*Poly(c, t)
    ln, ls = splitfactor(l, D, x, t)
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
    Qv = is_log_deriv_k_t_radical(nfmwa, nfmwd, D, x, t, 'auto')
    if Qv is None:
        raise NonElementaryIntegral("parametric_log_deriv_heu(): %s/%s " +
        "(N*f - M*w) is not the logarithmic derivaitive of a k(t)-radical."
        % (nfmwa, nfmwd))
    Q, v = Qv

    if Q.is_zero or v.is_zero:
        raise NonElementaryIntegral("parametric_log_deriv_heu(): Q == 0 " +
            "or v == 0.")

    return (Q*N, Q*M, v)

def parametric_log_deriv(fa, fd, wa, wd, D, x, t):
    # TODO: Write the full algorithm using the structure theorems.
    return parametric_log_deriv_heu(fa, fd, wa, wd, D, x, t)

def is_log_deriv_k_t_radical(fa, fd, D, x, t, case='auto'):
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
    from pudb import set_trace; set_trace() # Debugging
    if case == 'auto':
        case = get_case(D, x, t)

    # f must be simple
    n, s = splitfactor(fd, D, x, t)
    if not s.is_one:
        return None

    z = Symbol('z', dummy=True)
    H, b = residue_reduce(fa, fd, D, x, t, z=z)
    if not b:
        # Note, according to the note on page 255 of Bronstein's book, this
        # should never happen with the parametric logarithmic derivative
        # problems.  I don't know, however, if it will happen when building up
        # the differential field for an integrand.
        raise NotImplementedError("f has a non-elementary integral, cannot " +
            "determine if it is the logarithmic derivative of a k(t)-radical.")

    p = cancel(fa.as_basic()/fd.as_basic() - residue_reduce_derivation(H, D, x, t, z))
    try:
        p = Poly(p, t)
    except PolynomialError:
        # f - Dg will be in k[t] if f is the logarithmic derivaitve of a k(t)-radical
        return None

    if p.degree(t) >= max(1, D.degree(t)):
        return None

    if case == 'exp':
        wa, wd = derivation(t, D, x, t).cancel(Poly(t, t), include=True)
        try:
            n, e, u = parametric_log_deriv(p, Poly(1, t), wa, wd, D, x, t)
        except NonElementaryIntegral:
            return None

    elif case in ['tan', 'primitive']:
        raise NotImplementedError(r"The hypertangent and primitive cases are " +
        "not yet implemented for is_log_deriv_k_t_radical()")
    elif case in ['other_linear', 'other_nonlinear']:
        raise ValueError("The %s case is not supported in this function." % case)
    else:
        raise ValueError("case must be one of {'primitive', 'exp', 'tan', "+
        "'auto'}, not %s""" % case)

    return (n, u)

"""
Algorithms for solving Parametric Risch Differential Equations.

The methods used for solving Parametric Risch Differential Equations parallel
those for solving Risch Differential Equations.  See the outline in the
docstring of rde.py for more information.

The Parametric Risch Differential Equation problem is, given f, g1, ..., gm in
K(t), to determine if there exist y in K(t) and c1, ..., cm in Const(K) such
that Dy + f*y == Sum(ci*gi, (i, 1, m)), and to find such y and ci if they exist.

For the algorithms here G is a list of tuples of factions of the terms on the
right hand side of the equation (i.e., gi in k(t)), and Q is a list of terms on
the right hand side of the equation (i.e., qi in k[t]).  See the docstring of
each function for more information.
"""
from __future__ import print_function, division

from sympy.core import Dummy, ilcm, Add, Mul, Pow, S

from sympy.matrices import Matrix, zeros, eye

from sympy.solvers import solve

from sympy.polys import Poly, lcm, cancel, sqf_list

from sympy.integrals.risch import (gcdex_diophantine, frac_in, derivation,
    NonElementaryIntegralException, residue_reduce, splitfactor,
    residue_reduce_derivation, DecrementLevel, recognize_log_derivative)
from sympy.integrals.rde import (order_at, order_at_oo, weak_normalizer,
    bound_degree, spde, solve_poly_rde)
from sympy.core.compatibility import reduce, range
from sympy.utilities.misc import debug


def prde_normal_denom(fa, fd, G, DE):
    """
    Parametric Risch Differential Equation - Normal part of the denominator.

    Given a derivation D on k[t] and f, g1, ..., gm in k(t) with f weakly
    normalized with respect to t, return the tuple (a, b, G, h) such that
    a, h in k[t], b in k<t>, G = [g1, ..., gm] in k(t)^m, and for any solution
    c1, ..., cm in Const(k) and y in k(t) of Dy + f*y == Sum(ci*gi, (i, 1, m)),
    q == y*h in k<t> satisfies a*Dq + b*q == Sum(ci*Gi, (i, 1, m)).
    """
    dn, ds = splitfactor(fd, DE)
    Gas, Gds = list(zip(*G))
    gd = reduce(lambda i, j: i.lcm(j), Gds, Poly(1, DE.t))
    en, es = splitfactor(gd, DE)

    p = dn.gcd(en)
    h = en.gcd(en.diff(DE.t)).quo(p.gcd(p.diff(DE.t)))

    a = dn*h
    c = a*h

    ba = a*fa - dn*derivation(h, DE)*fd
    ba, bd = ba.cancel(fd, include=True)

    G = [(c*A).cancel(D, include=True) for A, D in G]

    return (a, (ba, bd), G, h)

def real_imag(ba, bd, gen):
    """
    Helper function, to get the real and imaginary part of a rational function
    evaluated at sqrt(-1) without actually evaluating it at sqrt(-1)

    Separates the even and odd power terms by checking the degree of terms wrt
    mod 4. Returns a tuple (ba[0], ba[1], bd) where ba[0] is real part
    of the numerator ba[1] is the imaginary part and bd is the denominator
    of the rational function.
    """
    bd = bd.as_poly(gen).as_dict()
    ba = ba.as_poly(gen).as_dict()
    denom_real = [value if key[0] % 4 == 0 else -value if key[0] % 4 == 2 else 0 for key, value in bd.items()]
    denom_imag = [value if key[0] % 4 == 1 else -value if key[0] % 4 == 3 else 0 for key, value in bd.items()]
    bd_real = sum(r for r in denom_real)
    bd_imag = sum(r for r in denom_imag)
    num_real = [value if key[0] % 4 == 0 else -value if key[0] % 4 == 2 else 0 for key, value in ba.items()]
    num_imag = [value if key[0] % 4 == 1 else -value if key[0] % 4 == 3 else 0 for key, value in ba.items()]
    ba_real = sum(r for r in num_real)
    ba_imag = sum(r for r in num_imag)
    ba = ((ba_real*bd_real + ba_imag*bd_imag).as_poly(gen), (ba_imag*bd_real - ba_real*bd_imag).as_poly(gen))
    bd = (bd_real*bd_real + bd_imag*bd_imag).as_poly(gen)
    return (ba[0], ba[1], bd)


def prde_special_denom(a, ba, bd, G, DE, case='auto'):
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
    # TODO: Merge this with the very similar special_denom() in rde.py
    if case == 'auto':
        case = DE.case

    if case == 'exp':
        p = Poly(DE.t, DE.t)
    elif case == 'tan':
        p = Poly(DE.t**2 + 1, DE.t)
    elif case in ['primitive', 'base']:
        B = ba.quo(bd)
        return (a, B, G, Poly(1, DE.t))
    else:
        raise ValueError("case must be one of {'exp', 'tan', 'primitive', "
            "'base'}, not %s." % case)

    nb = order_at(ba, p, DE.t) - order_at(bd, p, DE.t)
    nc = min([order_at(Ga, p, DE.t) - order_at(Gd, p, DE.t) for Ga, Gd in G])
    n = min(0, nc - min(0, nb))
    if not nb:
        # Possible cancellation.
        if case == 'exp':
            dcoeff = DE.d.quo(Poly(DE.t, DE.t))
            with DecrementLevel(DE):  # We are guaranteed to not have problems,
                                      # because case != 'base'.
                alphaa, alphad = frac_in(-ba.eval(0)/bd.eval(0)/a.eval(0), DE.t)
                etaa, etad = frac_in(dcoeff, DE.t)
                A = parametric_log_deriv(alphaa, alphad, etaa, etad, DE)
                if A is not None:
                    a, m, z = A
                    if a == 1:
                        n = min(n, m)

        elif case == 'tan':
            dcoeff = DE.d.quo(Poly(DE.t**2 + 1, DE.t))
            with DecrementLevel(DE):  # We are guaranteed to not have problems,
                                      # because case != 'base'.
                betaa, alphaa, alphad =  real_imag(ba, bd*a, DE.t)
                betad = alphad
                etaa, etad = frac_in(dcoeff, DE.t)
                if recognize_log_derivative(2*betaa, betad, DE):
                    A = parametric_log_deriv(alphaa, alphad, etaa, etad, DE)
                    B = parametric_log_deriv(betaa, betad, etaa, etad, DE)
                    if A is not None and B is not None:
                        a, s, z = A
                        if a == 1:
                             n = min(n, s/2)

    N = max(0, -nb)
    pN = p**N
    pn = p**-n  # This is 1/h

    A = a*pN
    B = ba*pN.quo(bd) + Poly(n, DE.t)*a*derivation(p, DE).quo(p)*pN
    G = [(Ga*pN*pn).cancel(Gd, include=True) for Ga, Gd in G]
    h = pn

    # (a*p**N, (b + n*a*Dp/p)*p**N, g1*p**(N - n), ..., gm*p**(N - n), p**-n)
    return (A, B, G, h)


def prde_linear_constraints(a, b, G, DE):
    """
    Parametric Risch Differential Equation - Generate linear constraints on the constants.

    Given a derivation D on k[t], a, b, in k[t] with gcd(a, b) == 1, and
    G = [g1, ..., gm] in k(t)^m, return Q = [q1, ..., qm] in k[t]^m and a
    matrix M with entries in k(t) such that for any solution c1, ..., cm in
    Const(k) and p in k[t] of a*Dp + b*p == Sum(ci*gi, (i, 1, m)),
    (c1, ..., cm) is a solution of Mx == 0, and p and the ci satisfy
    a*Dp + b*p == Sum(ci*qi, (i, 1, m)).

    Because M has entries in k(t), and because Matrix doesn't play well with
    Poly, M will be a Matrix of Basic expressions.
    """
    m = len(G)

    Gns, Gds = list(zip(*G))
    d = reduce(lambda i, j: i.lcm(j), Gds)
    d = Poly(d, field=True)
    Q = [(ga*(d).quo(gd)).div(d) for ga, gd in G]

    if not all([ri.is_zero for _, ri in Q]):
        N = max([ri.degree(DE.t) for _, ri in Q])
        M = Matrix(N + 1, m, lambda i, j: Q[j][1].nth(i))
    else:
        M = Matrix()  # No constraints, return the empty matrix.

    qs, _ = list(zip(*Q))
    return (qs, M)


def constant_system(A, u, DE):
    """
    Generate a system for the constant solutions.

    Given a differential field (K, D) with constant field C = Const(K), a Matrix
    A, and a vector (Matrix) u with coefficients in K, returns the tuple
    (B, v, s), where B is a Matrix with coefficients in C and v is a vector
    (Matrix) such that either v has coefficients in C, in which case s is True
    and the solutions in C of Ax == u are exactly all the solutions of Bx == v,
    or v has a non-constant coefficient, in which case s is False Ax == u has no
    constant solution.

    This algorithm is used both in solving parametric problems and in
    determining if an element a of K is a derivative of an element of K or the
    logarithmic derivative of a K-radical using the structure theorem approach.

    Because Poly does not play well with Matrix yet, this algorithm assumes that
    all matrix entries are Basic expressions.
    """
    if not A:
        return A, u
    Au = A.row_join(u)
    Au = Au.rref(simplify=cancel)[0]
    # Warning: This will NOT return correct results if cancel() cannot reduce
    # an identically zero expression to 0.  The danger is that we might
    # incorrectly prove that an integral is nonelementary (such as
    # risch_integrate(exp((sin(x)**2 + cos(x)**2 - 1)*x**2), x).
    # But this is a limitation in computer algebra in general, and implicit
    # in the correctness of the Risch Algorithm is the computability of the
    # constant field (actually, this same correctness problem exists in any
    # algorithm that uses rref()).
    #
    # We therefore limit ourselves to constant fields that are computable
    # via the cancel() function, in order to prevent a speed bottleneck from
    # calling some more complex simplification function (rational function
    # coefficients will fall into this class).  Furthermore, (I believe) this
    # problem will only crop up if the integral explicitly contains an
    # expression in the constant field that is identically zero, but cannot
    # be reduced to such by cancel().  Therefore, a careful user can avoid this
    # problem entirely by being careful with the sorts of expressions that
    # appear in his integrand in the variables other than the integration
    # variable (the structure theorems should be able to completely decide these
    # problems in the integration variable).

    Au = Au.applyfunc(cancel)
    A, u = Au[:, :-1], Au[:, -1]

    for j in range(A.cols):
        for i in range(A.rows):
            if A[i, j].has(*DE.T):
                # This assumes that const(F(t0, ..., tn) == const(K) == F
                Ri = A[i, :]
                # Rm+1; m = A.rows
                Rm1 = Ri.applyfunc(lambda x: derivation(x, DE, basic=True)/
                    derivation(A[i, j], DE, basic=True))
                Rm1 = Rm1.applyfunc(cancel)
                um1 = cancel(derivation(u[i], DE, basic=True)/
                    derivation(A[i, j], DE, basic=True))

                for s in range(A.rows):
                    # A[s, :] = A[s, :] - A[s, i]*A[:, m+1]
                    Asj = A[s, j]
                    A.row_op(s, lambda r, jj: cancel(r - Asj*Rm1[jj]))
                    # u[s] = u[s] - A[s, j]*u[m+1
                    u.row_op(s, lambda r, jj: cancel(r - Asj*um1))

                A = A.col_join(Rm1)
                u = u.col_join(Matrix([um1]))

    return (A, u)


def prde_spde(a, b, Q, n, DE):
    """
    Special Polynomial Differential Equation algorithm: Parametric Version.

    Given a derivation D on k[t], an integer n, and a, b, q1, ..., qm in k[t]
    with deg(a) > 0 and gcd(a, b) == 1, return (A, B, Q, R, n1), with
    Qq = [q1, ..., qm] and R = [r1, ..., rm], such that for any solution
    c1, ..., cm in Const(k) and q in k[t] of degree at most n of
    a*Dq + b*q == Sum(ci*gi, (i, 1, m)), p = (q - Sum(ci*ri, (i, 1, m)))/a has
    degree at most n1 and satisfies A*Dp + B*p == Sum(ci*qi, (i, 1, m))
    """
    R, Z = list(zip(*[gcdex_diophantine(b, a, qi) for qi in Q]))

    A = a
    B = b + derivation(a, DE)
    Qq = [zi - derivation(ri, DE) for ri, zi in zip(R, Z)]
    R = list(R)
    n1 = n - a.degree(DE.t)

    return (A, B, Qq, R, n1)


def prde_no_cancel_b_large(b, Q, n, DE):
    """
    Parametric Poly Risch Differential Equation - No cancellation: deg(b) large enough.

    Given a derivation D on k[t], n in ZZ, and b, q1, ..., qm in k[t] with
    b != 0 and either D == d/dt or deg(b) > max(0, deg(D) - 1), returns
    h1, ..., hr in k[r] and a matrix A with coefficients in Const(k) such that
    if c1, ..., cm in Const(k) and q in k[t] satisfy deg(q) <= n and
    Dq + b*Q == Sum(ci*qi, (i, 1, m)), then q = Sum(dj*hj, (j, 1, r)), where
    d1, ..., dr in Const(k) and A*Matrix([[c1, ..., cm, d1, ..., dr]]).T == 0.
    """
    db = b.degree(DE.t)
    m = len(Q)
    H = [Poly(0, DE.t)]*m

    for N in range(n, -1, -1):  # [n, ..., 0]
        for i in range(m):
            si = Q[i].nth(N + db)/b.LC()
            sitn = Poly(si*DE.t**N, DE.t)
            H[i] = H[i] + sitn
            Q[i] = Q[i] - derivation(sitn, DE) - b*sitn

    if all(qi.is_zero for qi in Q):
        dc = -1
        M = zeros(0, 2)
    else:
        dc = max([qi.degree(t) for qi in Q])
        M = Matrix(dc + 1, m, lambda i, j: Q[j].nth(i))
    A, u = constant_system(M, zeros(dc + 1, 1), DE)
    c = eye(m)
    A = A.row_join(zeros(A.rows, m)).col_join(c.row_join(-c))

    return (H, A)


def prde_no_cancel_b_small(b, Q, n, DE):
    """
    Parametric Poly Risch Differential Equation - No cancellation: deg(b) small enough.

    Given a derivation D on k[t], n in ZZ, and b, q1, ..., qm in k[t] with
    deg(b) < deg(D) - 1 and either D == d/dt or deg(D) >= 2, returns
    h1, ..., hr in k[t] and a matrix A with coefficients in Const(k) such that
    if c1, ..., cm in Const(k) and q in k[t] satisfy deg(q) <= n and
    Dq + b*q == Sum(ci*qi, (i, 1, m)) then q = Sum(dj*hj, (j, 1, r)) where
    d1, ..., dr in Const(k) and A*Matrix([[c1, ..., cm, d1, ..., dr]]).T == 0.
    """
    m = len(Q)
    H = [Poly(0, DE.t)]*m

    for N in range(n, 0, -1):  # [n, ..., 1]
        for i in range(m):
            si = Q[i].nth(N + DE.d.degree(DE.t) - 1)/(N*DE.d.LC())
            sitn = Poly(si*DE.t**N, DE.t)
            H[i] = H[i] + sitn
            Q[i] = Q[i] - derivation(sitn, DE) - b*sitn

    if b.degree(DE.t) > 0:
        for i in range(m):
            si = Poly(Q[i].nth(b.degree(DE.t))/b.LC(), DE.t)
            H[i] = H[i] + si
            Q[i] = Q[i] - derivation(si, DE) - b*si
        if all(qi.is_zero for qi in Q):
            dc = -1
            M = Matrix()
        else:
            dc = max([qi.degree(DE.t) for qi in Q])
            M = Matrix(dc + 1, m, lambda i, j: Q[j].nth(i))
        A, u = constant_system(M, zeros(dc + 1, 1), DE)
        c = eye(m)
        A = A.row_join(zeros(A.rows, m)).col_join(c.row_join(-c))
        return (H, A)
    else:
        # TODO: implement this (requires recursive param_rischDE() call)
        raise NotImplementedError


def param_rischDE(fa, fd, G, DE):
    """
    Solve a Parametric Risch Differential Equation: Dy + f*y == Sum(ci*Gi, (i, 1, m)).
    """
    _, (fa, fd) = weak_normalizer(fa, fd, DE)
    a, (ba, bd), G, hn = prde_normal_denom(ga, gd, G, DE)
    A, B, G, hs = prde_special_denom(a, ba, bd, G, DE)
    g = gcd(A, B)
    A, B, G = A.quo(g), B.quo(g), [gia.cancel(gid*g, include=True) for
        gia, gid in G]
    Q, M = prde_linear_constraints(A, B, G, DE)
    M, _ = constant_system(M, zeros(M.rows, 1), DE)
    # Reduce number of constants at this point
    try:
        # Similar to rischDE(), we try oo, even though it might lead to
        # non-termination when there is no solution.  At least for prde_spde,
        # it will always terminate no matter what n is.
        n = bound_degree(A, B, G, DE, parametric=True)
    except NotImplementedError:
        debug("param_rischDE: Proceeding with n = oo; may cause "
              "non-termination.")
        n = oo

    A, B, Q, R, n1 = prde_spde(A, B, Q, n, DE)


def limited_integrate_reduce(fa, fd, G, DE):
    """
    Simpler version of step 1 & 2 for the limited integration problem.

    Given a derivation D on k(t) and f, g1, ..., gn in k(t), return
    (a, b, h, N, g, V) such that a, b, h in k[t], N is a non-negative integer,
    g in k(t), V == [v1, ..., vm] in k(t)^m, and for any solution v in k(t),
    c1, ..., cm in C of f == Dv + Sum(ci*wi, (i, 1, m)), p = v*h is in k<t>, and
    p and the ci satisfy a*Dp + b*p == g + Sum(ci*vi, (i, 1, m)).  Furthermore,
    if S1irr == Sirr, then p is in k[t], and if t is nonlinear or Liouvillian
    over k, then deg(p) <= N.

    So that the special part is always computed, this function calls the more
    general prde_special_denom() automatically if it cannot determine that
    S1irr == Sirr.  Furthermore, it will automatically call bound_degree() when
    t is linear and non-Liouvillian, which for the transcendental case, implies
    that Dt == a*t + b with for some a, b in k*.
    """
    dn, ds = splitfactor(fd, DE)
    E = [splitfactor(gd, DE) for _, gd in G]
    En, Es = list(zip(*E))
    c = reduce(lambda i, j: i.lcm(j), (dn,) + En)  # lcm(dn, en1, ..., enm)
    hn = c.gcd(c.diff(DE.t))
    a = hn
    b = -derivation(hn, DE)
    N = 0

    # These are the cases where we know that S1irr = Sirr, but there could be
    # others, and this algorithm will need to be extended to handle them.
    if DE.case in ['base', 'primitive', 'exp', 'tan']:
        hs = reduce(lambda i, j: i.lcm(j), (ds,) + Es)  # lcm(ds, es1, ..., esm)
        a = hn*hs
        b = -derivation(hn, DE) - (hn*derivation(hs, DE)).quo(hs)
        mu = min(order_at_oo(fa, fd, DE.t), min([order_at_oo(ga, gd, DE.t) for
            ga, gd in G]))
        # So far, all the above are also nonlinear or Liouvillian, but if this
        # changes, then this will need to be updated to call bound_degree()
        # as per the docstring of this function (DE.case == 'other_linear').
        N = hn.degree(DE.t) + hs.degree(DE.t) + max(0, 1 - DE.d.degree(DE.t) - mu)
    else:
        # TODO: implement this
        raise NotImplementedError

    V = [(-a*hn*ga).cancel(gd, include=True) for ga, gd in G]
    return (a, b, a, N, (a*hn*fa).cancel(fd, include=True), V)


def limited_integrate(fa, fd, G, DE):
    """
    Solves the limited integration problem:  f = Dv + Sum(ci*wi, (i, 1, n))
    """
    fa, fd = fa*Poly(1/fd.LC(), DE.t), fd.monic()
    A, B, h, N, g, V = limited_integrate_reduce(fa, fd, G, DE)
    V = [g] + V
    g = A.gcd(B)
    A, B, V = A.quo(g), B.quo(g), [via.cancel(vid*g, include=True) for
        via, vid in V]
    Q, M = prde_linear_constraints(A, B, V, DE)
    M, _ = constant_system(M, zeros(M.rows, 1), DE)
    l = M.nullspace()
    if M == Matrix() or len(l) > 1:
        # Continue with param_rischDE()
        raise NotImplementedError("param_rischDE() is required to solve this "
            "integral.")
    elif len(l) == 0:
        raise NonElementaryIntegralException
    elif len(l) == 1:
        # The c1 == 1.  In this case, we can assume a normal Risch DE
        if l[0][0].is_zero:
            raise NonElementaryIntegralException
        else:
            l[0] *= 1/l[0][0]
            C = sum([Poly(i, DE.t)*q for (i, q) in zip(l[0], Q)])
            # Custom version of rischDE() that uses the already computed
            # denominator and degree bound from above.
            B, C, m, alpha, beta = spde(A, B, C, N, DE)
            y = solve_poly_rde(B, C, m, DE)

            return ((alpha*y + beta, h), list(l[0][1:]))
    else:
        raise NotImplementedError


def parametric_log_deriv_heu(fa, fd, wa, wd, DE, c1=None):
    """
    Parametric logarithmic derivative heuristic.

    Given a derivation D on k[t], f in k(t), and a hyperexponential monomial
    theta over k(t), raises either NotImplementedError, in which case the
    heuristic failed, or returns None, in which case it has proven that no
    solution exists, or returns a solution (n, m, v) of the equation
    n*f == Dv/v + m*Dtheta/theta, with v in k(t)* and n, m in ZZ with n != 0.

    If this heuristic fails, the structure theorem approach will need to be
    used.

    The argument w == Dtheta/theta
    """
    # TODO: finish writing this and write tests
    c1 = c1 or Dummy('c1')

    p, a = fa.div(fd)
    q, b = wa.div(wd)

    B = max(0, derivation(DE.t, DE).degree(DE.t) - 1)
    C = max(p.degree(DE.t), q.degree(DE.t))

    if q.degree(DE.t) > B:
        eqs = [p.nth(i) - c1*q.nth(i) for i in range(B + 1, C + 1)]
        s = solve(eqs, c1)
        if not s or not s[c1].is_Rational:
            # deg(q) > B, no solution for c.
            return None

        N, M = s[c1].as_numer_denom()  # N and M are integers
        N, M = Poly(N, DE.t), Poly(M, DE.t)

        nfmwa = N*fa*wd - M*wa*fd
        nfmwd = fd*wd
        Qv = is_log_deriv_k_t_radical_in_field(N*fa*wd - M*wa*fd, fd*wd, DE,
            'auto')
        if Qv is None:
            # (N*f - M*w) is not the logarithmic derivative of a k(t)-radical.
            return None

        Q, e, v = Qv
        if e != 1:
            return None

        if Q.is_zero or v.is_zero:
            return None

        return (Q*N, Q*M, v)

    if p.degree(DE.t) > B:
        return None

    c = lcm(fd.as_poly(DE.t).LC(), wd.as_poly(DE.t).LC())
    l = fd.monic().lcm(wd.monic())*Poly(c, DE.t)
    ln, ls = splitfactor(l, DE)
    z = ls*ln.gcd(ln.diff(DE.t))

    if not z.has(DE.t):
        raise NotImplementedError("parametric_log_deriv_heu() "
            "heuristic failed: z in k.")

    u1, r1 = (fa*l.quo(fd)).div(z)  # (l*f).div(z)
    u2, r2 = (wa*l.quo(wd)).div(z)  # (l*w).div(z)

    eqs = [r1.nth(i) - c1*r2.nth(i) for i in range(z.degree(DE.t))]
    s = solve(eqs, c1)
    if not s or not s[c1].is_Rational:
        # deg(q) <= B, no solution for c.
        return None

    M, N = s[c1].as_numer_denom()

    nfmwa = N.as_poly(DE.t)*fa*wd - M.as_poly(DE.t)*wa*fd
    nfmwd = fd*wd
    Qv = is_log_deriv_k_t_radical_in_field(nfmwa, nfmwd, DE)
    if Qv is None:
        # (N*f - M*w) is not the logarithmic derivative of a k(t)-radical.
        return None

    Q, v = Qv

    if Q.is_zero or v.is_zero:
        return None

    return (Q*N, Q*M, v)


def parametric_log_deriv(fa, fd, wa, wd, DE):
    # TODO: Write the full algorithm using the structure theorems.
#    try:
    A = parametric_log_deriv_heu(fa, fd, wa, wd, DE)
#    except NotImplementedError:
        # Heuristic failed, we have to use the full method.
        # TODO: This could be implemented more efficiently.
        # It isn't too worrisome, because the heuristic handles most difficult
        # cases.
    return A


def is_deriv_k(fa, fd, DE):
    """
    Checks if Df/f is the derivative of an element of k(t).

    a in k(t) is the derivative of an element of k(t) if there exists b in k(t)
    such that a = Db.  Either returns (ans, u), such that Df/f == Du, or None,
    which means that Df/f is not the derivative of an element of k(t).  ans is
    a list of tuples such that Add(*[i*j for i, j in ans]) == u.  This is useful
    for seeing exactly which elements of k(t) produce u.

    This function uses the structure theorem approach, which says that for any
    f in K, Df/f is the derivative of a element of K if and only if there are ri
    in QQ such that::

            ---               ---       Dt
            \    r  * Dt   +  \    r  *   i      Df
            /     i     i     /     i   ---   =  --.
            ---               ---        t        f
         i in L            i in E         i
               K/C(x)            K/C(x)


    Where C = Const(K), L_K/C(x) = { i in {1, ..., n} such that t_i is
    transcendental over C(x)(t_1, ..., t_i-1) and Dt_i = Da_i/a_i, for some a_i
    in C(x)(t_1, ..., t_i-1)* } (i.e., the set of all indices of logarithmic
    monomials of K over C(x)), and E_K/C(x) = { i in {1, ..., n} such that t_i
    is transcendental over C(x)(t_1, ..., t_i-1) and Dt_i/t_i = Da_i, for some
    a_i in C(x)(t_1, ..., t_i-1) } (i.e., the set of all indices of
    hyperexponential monomials of K over C(x)).  If K is an elementary extension
    over C(x), then the cardinality of L_K/C(x) U E_K/C(x) is exactly the
    transcendence degree of K over C(x).  Furthermore, because Const_D(K) ==
    Const_D(C(x)) == C, deg(Dt_i) == 1 when t_i is in E_K/C(x) and
    deg(Dt_i) == 0 when t_i is in L_K/C(x), implying in particular that E_K/C(x)
    and L_K/C(x) are disjoint.

    The sets L_K/C(x) and E_K/C(x) must, by their nature, be computed
    recursively using this same function.  Therefore, it is required to pass
    them as indices to D (or T).  E_args are the arguments of the
    hyperexponentials indexed by E_K (i.e., if i is in E_K, then T[i] ==
    exp(E_args[i])).  This is needed to compute the final answer u such that
    Df/f == Du.

    log(f) will be the same as u up to a additive constant.  This is because
    they will both behave the same as monomials. For example, both log(x) and
    log(2*x) == log(x) + log(2) satisfy Dt == 1/x, because log(2) is constant.
    Therefore, the term const is returned.  const is such that
    log(const) + f == u.  This is calculated by dividing the arguments of one
    logarithm from the other.  Therefore, it is necessary to pass the arguments
    of the logarithmic terms in L_args.

    To handle the case where we are given Df/f, not f, use is_deriv_k_in_field().
    """
    # Compute Df/f
    dfa, dfd = fd*(fd*derivation(fa, DE) - fa*derivation(fd, DE)), fd**2*fa
    dfa, dfd = dfa.cancel(dfd, include=True)

    # Our assumption here is that each monomial is recursively transcendental
    if len(DE.L_K) + len(DE.E_K) != len(DE.D) - 1:
        if [i for i in DE.cases if i == 'tan'] or \
                set([i for i in DE.cases if i == 'primitive']) - set(DE.L_K):
            raise NotImplementedError("Real version of the structure "
                "theorems with hypertangent support is not yet implemented.")

        # TODO: What should really be done in this case?
        raise NotImplementedError("Nonelementary extensions not supported "
            "in the structure theorems.")

    E_part = [DE.D[i].quo(Poly(DE.T[i], DE.T[i])).as_expr() for i in DE.E_K]
    L_part = [DE.D[i].as_expr() for i in DE.L_K]

    lhs = Matrix([E_part + L_part])
    rhs = Matrix([dfa.as_expr()/dfd.as_expr()])

    A, u = constant_system(lhs, rhs, DE)

    if not all(derivation(i, DE, basic=True).is_zero for i in u) or not A:
        # If the elements of u are not all constant
        # Note: See comment in constant_system

        # Also note: derivation(basic=True) calls cancel()
        return None
    else:
        if not all(i.is_Rational for i in u):
            raise NotImplementedError("Cannot work with non-rational "
                "coefficients in this case.")
        else:
            terms = DE.E_args + [DE.T[i] for i in DE.L_K]
            ans = list(zip(terms, u))
            result = Add(*[Mul(i, j) for i, j in ans])
            argterms = [DE.T[i] for i in DE.E_K] + DE.L_args
            l = []
            ld = []
            for i, j in zip(argterms, u):
                # We need to get around things like sqrt(x**2) != x
                # and also sqrt(x**2 + 2*x + 1) != x + 1
                # Issue 10798: i need not be a polynomial
                i, d = i.as_numer_denom()
                icoeff, iterms = sqf_list(i)
                l.append(Mul(*([Pow(icoeff, j)] + [Pow(b, e*j) for b, e in iterms])))
                dcoeff, dterms = sqf_list(d)
                ld.append(Mul(*([Pow(dcoeff, j)] + [Pow(b, e*j) for b, e in dterms])))
            const = cancel(fa.as_expr()/fd.as_expr()/Mul(*l)*Mul(*ld))

            return (ans, result, const)


def is_log_deriv_k_t_radical(fa, fd, DE, Df=True):
    """
    Checks if Df is the logarithmic derivative of a k(t)-radical.

    b in k(t) can be written as the logarithmic derivative of a k(t) radical if
    there exist n in ZZ and u in k(t) with n, u != 0 such that n*b == Du/u.
    Either returns (ans, u, n, const) or None, which means that Df cannot be
    written as the logarithmic derivative of a k(t)-radical.  ans is a list of
    tuples such that Mul(*[i**j for i, j in ans]) == u.  This is useful for
    seeing exactly what elements of k(t) produce u.

    This function uses the structure theorem approach, which says that for any
    f in K, Df is the logarithmic derivative of a K-radical if and only if there
    are ri in QQ such that::

            ---               ---       Dt
            \    r  * Dt   +  \    r  *   i
            /     i     i     /     i   ---   =  Df.
            ---               ---        t
         i in L            i in E         i
               K/C(x)            K/C(x)


    Where C = Const(K), L_K/C(x) = { i in {1, ..., n} such that t_i is
    transcendental over C(x)(t_1, ..., t_i-1) and Dt_i = Da_i/a_i, for some a_i
    in C(x)(t_1, ..., t_i-1)* } (i.e., the set of all indices of logarithmic
    monomials of K over C(x)), and E_K/C(x) = { i in {1, ..., n} such that t_i
    is transcendental over C(x)(t_1, ..., t_i-1) and Dt_i/t_i = Da_i, for some
    a_i in C(x)(t_1, ..., t_i-1) } (i.e., the set of all indices of
    hyperexponential monomials of K over C(x)).  If K is an elementary extension
    over C(x), then the cardinality of L_K/C(x) U E_K/C(x) is exactly the
    transcendence degree of K over C(x).  Furthermore, because Const_D(K) ==
    Const_D(C(x)) == C, deg(Dt_i) == 1 when t_i is in E_K/C(x) and
    deg(Dt_i) == 0 when t_i is in L_K/C(x), implying in particular that E_K/C(x)
    and L_K/C(x) are disjoint.

    The sets L_K/C(x) and E_K/C(x) must, by their nature, be computed
    recursively using this same function.  Therefore, it is required to pass
    them as indices to D (or T).  L_args are the arguments of the logarithms
    indexed by L_K (i.e., if i is in L_K, then T[i] == log(L_args[i])).  This is
    needed to compute the final answer u such that n*f == Du/u.

    exp(f) will be the same as u up to a multiplicative constant.  This is
    because they will both behave the same as monomials.  For example, both
    exp(x) and exp(x + 1) == E*exp(x) satisfy Dt == t. Therefore, the term const
    is returned.  const is such that exp(const)*f == u.  This is calculated by
    subtracting the arguments of one exponential from the other.  Therefore, it
    is necessary to pass the arguments of the exponential terms in E_args.

    To handle the case where we are given Df, not f, use
    is_log_deriv_k_t_radical_in_field().
    """
    H = []
    if Df:
        dfa, dfd = (fd*derivation(fa, DE) - fa*derivation(fd, DE)).cancel(fd**2,
            include=True)
    else:
        dfa, dfd = fa, fd

    # Our assumption here is that each monomial is recursively transcendental
    if len(DE.L_K) + len(DE.E_K) != len(DE.D) - 1:
        if [i for i in DE.cases if i == 'tan'] or \
                set([i for i in DE.cases if i == 'primitive']) - set(DE.L_K):
            raise NotImplementedError("Real version of the structure "
                "theorems with hypertangent support is not yet implemented.")

        # TODO: What should really be done in this case?
        raise NotImplementedError("Nonelementary extensions not supported "
            "in the structure theorems.")

    E_part = [DE.D[i].quo(Poly(DE.T[i], DE.T[i])).as_expr() for i in DE.E_K]
    L_part = [DE.D[i].as_expr() for i in DE.L_K]

    lhs = Matrix([E_part + L_part])
    rhs = Matrix([dfa.as_expr()/dfd.as_expr()])

    A, u = constant_system(lhs, rhs, DE)
    if not all(derivation(i, DE, basic=True).is_zero for i in u) or not A:
        # If the elements of u are not all constant
        # Note: See comment in constant_system

        # Also note: derivation(basic=True) calls cancel()
        return None
    else:
        if not all(i.is_Rational for i in u):
            # TODO: But maybe we can tell if they're not rational, like
            # log(2)/log(3). Also, there should be an option to continue
            # anyway, even if the result might potentially be wrong.
            raise NotImplementedError("Cannot work with non-rational "
                "coefficients in this case.")
        else:
            n = reduce(ilcm, [i.as_numer_denom()[1] for i in u])
            u *= n
            terms = [DE.T[i] for i in DE.E_K] + DE.L_args
            ans = list(zip(terms, u))
            result = Mul(*[Pow(i, j) for i, j in ans])

            # exp(f) will be the same as result up to a multiplicative
            # constant.  We now find the log of that constant.
            argterms = DE.E_args + [DE.T[i] for i in DE.L_K]
            const = cancel(fa.as_expr()/fd.as_expr() -
                Add(*[Mul(i, j/n) for i, j in zip(argterms, u)]))

            return (ans, result, n, const)


def is_log_deriv_k_t_radical_in_field(fa, fd, DE, case='auto', z=None):
    """
    Checks if f can be written as the logarithmic derivative of a k(t)-radical.

    f in k(t) can be written as the logarithmic derivative of a k(t) radical if
    there exist n in ZZ and u in k(t) with n, u != 0 such that n*f == Du/u.
    Either returns (n, u) or None, which means that f cannot be written as the
    logarithmic derivative of a k(t)-radical.

    case is one of {'primitive', 'exp', 'tan', 'auto'} for the primitive,
    hyperexponential, and hypertangent cases, respectively.  If case is 'auto',
    it will attempt to determine the type of the derivation automatically.
    """
    fa, fd = fa.cancel(fd, include=True)

    # f must be simple
    n, s = splitfactor(fd, DE)
    if not s.is_one:
        pass

    z = z or Dummy('z')
    H, b = residue_reduce(fa, fd, DE, z=z)
    if not b:
        # I will have to verify, but I believe that the answer should be
        # None in this case. This should never happen for the
        # functions given when solving the parametric logarithmic
        # derivative problem when integration elementary functions (see
        # Bronstein's book, page 255), so most likely this indicates a bug.
        return None

    roots = [(i, i.real_roots()) for i, _ in H]
    if not all(len(j) == i.degree() and all(k.is_Rational for k in j) for
               i, j in roots):
        # If f is the logarithmic derivative of a k(t)-radical, then all the
        # roots of the resultant must be rational numbers.
        return None

    # [(a, i), ...], where i*log(a) is a term in the log-part of the integral
    # of f
    respolys, residues = list(zip(*roots)) or [[], []]
    # Note: this might be empty, but everything below should work find in that
    # case (it should be the same as if it were [[1, 1]])
    residueterms = [(H[j][1].subs(z, i), i) for j in range(len(H)) for
        i in residues[j]]

    # TODO: finish writing this and write tests

    p = cancel(fa.as_expr()/fd.as_expr() - residue_reduce_derivation(H, DE, z))

    p = p.as_poly(DE.t)
    if p is None:
        # f - Dg will be in k[t] if f is the logarithmic derivative of a k(t)-radical
        return None

    if p.degree(DE.t) >= max(1, DE.d.degree(DE.t)):
        return None

    if case == 'auto':
        case = DE.case

    if case == 'exp':
        wa, wd = derivation(DE.t, DE).cancel(Poly(DE.t, DE.t), include=True)
        with DecrementLevel(DE):
            pa, pd = frac_in(p, DE.t, cancel=True)
            wa, wd = frac_in((wa, wd), DE.t)
            A = parametric_log_deriv(pa, pd, wa, wd, DE)
        if A is None:
            return None
        n, e, u = A
        u *= DE.t**e

    elif case == 'primitive':
        with DecrementLevel(DE):
            pa, pd = frac_in(p, DE.t)
            A = is_log_deriv_k_t_radical_in_field(pa, pd, DE, case='auto')
        if A is None:
            return None
        n, u = A

    elif case == 'base':
        # TODO: we can use more efficient residue reduction from ratint()
        if not fd.is_sqf or fa.degree() >= fd.degree():
            # f is the logarithmic derivative in the base case if and only if
            # f = fa/fd, fd is square-free, deg(fa) < deg(fd), and
            # gcd(fa, fd) == 1.  The last condition is handled by cancel() above.
            return None
        # Note: if residueterms = [], returns (1, 1)
        # f had better be 0 in that case.
        n = reduce(ilcm, [i.as_numer_denom()[1] for _, i in residueterms], S(1))
        u = Mul(*[Pow(i, j*n) for i, j in residueterms])
        return (n, u)

    elif case == 'tan':
        raise NotImplementedError("The hypertangent case is "
        "not yet implemented for is_log_deriv_k_t_radical_in_field()")

    elif case in ['other_linear', 'other_nonlinear']:
        # XXX: If these are supported by the structure theorems, change to NotImplementedError.
        raise ValueError("The %s case is not supported in this function." % case)

    else:
        raise ValueError("case must be one of {'primitive', 'exp', 'tan', "
        "'base', 'auto'}, not %s" % case)

    common_denom = reduce(ilcm, [i.as_numer_denom()[1] for i in [j for _, j in
        residueterms]] + [n], S(1))
    residueterms = [(i, j*common_denom) for i, j in residueterms]
    m = common_denom//n
    if common_denom != n*m:  # Verify exact division
        raise ValueError("Inexact division")
    u = cancel(u**m*Mul(*[Pow(i, j) for i, j in residueterms]))

    return (common_denom, u)

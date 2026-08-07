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
from __future__ import annotations
import itertools
from functools import reduce

from sympy.core.intfunc import ilcm, igcd
from sympy.core import Dummy, Add, Mul, Pow, S
from sympy.integrals.rde import (order_at, order_at_oo, weak_normalizer,
    bound_degree, _special_denom_cancel_bound)
from sympy.integrals.risch import (gcdex_diophantine, frac_in, derivation,
    residue_reduce, splitfactor, residue_reduce_derivation, DecrementLevel)
from sympy.polys import Poly, lcm, cancel, sqf_list
from sympy.polys.polymatrix import PolyMatrix as Matrix
from sympy.solvers import solve

zeros = Matrix.zeros
eye = Matrix.eye


def prde_normal_denom(fa, fd, G, DE):
    """
    Parametric Risch Differential Equation - Normal part of the denominator.

    Explanation
    ===========

    Given a derivation D on k[t] and f, g1, ..., gm in k(t) with f weakly
    normalized with respect to t, return the tuple (a, b, G, h) such that
    a, h in k[t], b in k<t>, G = [g1, ..., gm] in k(t)^m, and for any solution
    c1, ..., cm in Const(k) and y in k(t) of Dy + f*y == Sum(ci*gi, (i, 1, m)),
    q == y*h in k<t> satisfies a*Dq + b*q == Sum(ci*Gi, (i, 1, m)).

    This is ``ParamRdeNormalDenominator`` from Section 7.1 of Bronstein's
    book.
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
    evaluated at sqrt(-1) without actually evaluating it at sqrt(-1).

    Explanation
    ===========

    Separates the even and odd power terms by checking the degree of terms wrt
    mod 4. Returns a tuple (ba[0], ba[1], bd) where ba[0] is real part
    of the numerator ba[1] is the imaginary part and bd is the denominator
    of the rational function.
    """
    bd = bd.as_poly(gen).as_dict()
    ba = ba.as_poly(gen).as_dict()
    denom_real = [value if key[0] % 4 == 0 else -value if key[0] % 4 == 2 else 0 for key, value in bd.items()]
    denom_imag = [value if key[0] % 4 == 1 else -value if key[0] % 4 == 3 else 0 for key, value in bd.items()]
    bd_real = sum(denom_real, S.Zero)
    bd_imag = sum(denom_imag, S.Zero)
    num_real = [value if key[0] % 4 == 0 else -value if key[0] % 4 == 2 else 0 for key, value in ba.items()]
    num_imag = [value if key[0] % 4 == 1 else -value if key[0] % 4 == 3 else 0 for key, value in ba.items()]
    ba_real = sum(num_real, S.Zero)
    ba_imag = sum(num_imag, S.Zero)
    ba = (Poly(ba_real*bd_real + ba_imag*bd_imag, gen), Poly(ba_imag*bd_real - ba_real*bd_imag, gen))
    bd = Poly(bd_real*bd_real + bd_imag*bd_imag, gen)
    return (ba[0], ba[1], bd)


def prde_special_denom(a, ba, bd, G, DE, case='auto'):
    """
    Parametric Risch Differential Equation - Special part of the denominator.

    Explanation
    ===========

    Case is one of {'exp', 'tan', 'primitive'} for the hyperexponential,
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

    This is ``ParamRdeSpecialDenomExp`` and ``ParamRdeSpecialDenomTan``
    from Section 7.1 of Bronstein's book.
    """
    # The cancellation-case bound is shared with special_denom() in rde.py
    # via _special_denom_cancel_bound().  Note that N below is
    # max(0, -nb), following the parametric algorithms of Section 7.1,
    # while the non-parametric special_denom() uses max(0, -nb, n - nc)
    # per Section 6.2; this difference is intentional (it matches the
    # book), not a divergence.
    if case == 'auto':
        case = DE.case

    if case == 'exp':
        p = Poly(DE.t, DE.t)
    elif case == 'tan':
        p = Poly(DE.t**2 + 1, DE.t)
    elif case in ('primitive', 'base'):
        B = ba.quo(bd)
        return (a, B, G, Poly(1, DE.t))
    elif case in ('other_linear', 'other_nonlinear'):
        raise NotImplementedError("The %s case is not implemented in "
            "prde_special_denom()." % case)
    else:
        raise ValueError("case must be one of {'exp', 'tan', 'primitive', "
            "'base', 'other_linear', 'other_nonlinear'}, not %s." % case)

    nb = order_at(ba, p, DE.t) - order_at(bd, p, DE.t)
    nc = min(order_at(Ga, p, DE.t) - order_at(Gd, p, DE.t) for Ga, Gd in G)
    n = min(0, nc - min(0, nb))
    if not nb:
        # Possible cancellation.
        n = _special_denom_cancel_bound(a, ba, bd, n, DE, case)

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

    Explanation
    ===========

    Given a derivation D on k[t], a, b, in k[t] with gcd(a, b) == 1, and
    G = [g1, ..., gm] in k(t)^m, return Q = [q1, ..., qm] in k[t]^m and a
    matrix M with entries in k(t) such that for any solution c1, ..., cm in
    Const(k) and p in k[t] of a*Dp + b*p == Sum(ci*gi, (i, 1, m)),
    (c1, ..., cm) is a solution of Mx == 0, and p and the ci satisfy
    a*Dp + b*p == Sum(ci*qi, (i, 1, m)).

    Because M has entries in k(t), and because Matrix does not play well with
    Poly, M will be a Matrix of Basic expressions.

    This is ``LinearConstraints`` from Section 7.1 of Bronstein's book.
    """
    m = len(G)

    Gns, Gds = list(zip(*G))
    d = reduce(lambda i, j: i.lcm(j), Gds)
    d = Poly(d, field=True)
    Q = [(ga*(d).quo(gd)).div(d) for ga, gd in G]

    if not all(ri.is_zero for _, ri in Q):
        N = max(ri.degree(DE.t) for _, ri in Q)
        M = Matrix(N + 1, m, lambda i, j: Q[j][1].nth(i), DE.t)
    else:
        M = Matrix(0, m, [], DE.t)  # No constraints, return the empty matrix.

    qs, _ = list(zip(*Q))
    return (qs, M)

def poly_linear_constraints(p, d):
    """
    Given p = [p1, ..., pm] in k[t]^m and d in k[t], return
    q = [q1, ..., qm] in k[t]^m and a matrix M with entries in k such
    that Sum(ci*pi, (i, 1, m)), for c1, ..., cm in k, is divisible
    by d if and only if (c1, ..., cm) is a solution of Mx = 0, in
    which case the quotient is Sum(ci*qi, (i, 1, m)).

    There is no named counterpart in Bronstein's book; this is a subroutine
    used by the Section 7.1 algorithms.
    """
    m = len(p)
    q, r = zip(*[pi.div(d) for pi in p])

    if not all(ri.is_zero for ri in r):
        n = max(ri.degree() for ri in r)
        M = Matrix(n + 1, m, lambda i, j: r[j].nth(i), d.gens)
    else:
        M = Matrix(0, m, [], d.gens)  # No constraints.

    return q, M

def constant_system(A, u, DE):
    """
    Generate a system for the constant solutions.

    Explanation
    ===========

    Given a differential field (K, D) with constant field C = Const(K), a Matrix
    A, and a vector (Matrix) u with coefficients in K, returns the tuple
    (B, v), where B is a Matrix with coefficients in C and v is a vector
    (Matrix) such that either v has coefficients in C, in which case the
    solutions in C of Ax == u are exactly all the solutions of Bx == v, or v
    has a non-constant coefficient, in which case Ax == u has no constant
    solution.  Note that B and v are a reduced *system*, not a solution:
    callers must check that v is constant and then solve Bx == v themselves.

    This algorithm is used both in solving parametric problems and in
    determining if an element a of K is a derivative of an element of K or the
    logarithmic derivative of a K-radical using the structure theorem approach.

    Because Poly does not play well with Matrix yet, this algorithm assumes that
    all matrix entries are Basic expressions.

    This is ``ConstantSystem`` from Section 7.1 of Bronstein's book.
    """
    if not A:
        return A, u
    Au = A.row_join(u)
    Au, _ = Au.rref()
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

    A, u = Au[:, :-1], Au[:, -1]

    D = lambda x: derivation(x, DE, basic=True)

    # Repeat until no entry of A is nonconstant, processing rows appended
    # in earlier iterations as well ("while A is not constant do" in the
    # book's ConstantSystem; iterating over a snapshot of the original
    # rows would let nonconstant entries in the appended rows escape).
    # The safety cap guards against nontermination in case cancel() fails
    # to recognize an identically zero derivation (see the comment above).
    for _ in range((Au.rows + 1)*(Au.cols + 1)*(len(DE.T) + 1) + 10):
        for j, i in itertools.product(range(A.cols), range(A.rows)):
            if A[i, j].expr.has(*DE.T):
                break
        else:
            break
        # This assumes that const(F(t0, ..., tn)) == const(K) == F
        # TODO: If A[i, j] contains symbolic constants, D(A[i, j]) can
        # vanish for special values of them (e.g. an entry y*t is
        # nonconstant except at y == 0), and dividing by it makes the
        # reduced system valid only generically.  The correct result
        # would be a Piecewise over the (solve()-generated) vanishing
        # conditions of each pivot.
        Ri = A[i, :]
        # Rm+1; m = A.rows
        DAij = D(A[i, j])
        Rm1 = Ri.applyfunc(lambda x: D(x) / DAij)
        um1 = D(u[i]) / DAij

        Aj = A[:, j]
        A = A - Aj * Rm1
        u = u - Aj * um1

        A = A.col_join(Rm1)
        u = u.col_join(Matrix([um1], u.gens))
    else:
        raise RuntimeError("constant_system() did not reach a constant "
            "system; the constant field may not be computable by cancel().")

    return (A, u)


def prde_spde(a, b, Q, n, DE):
    """
    Special Polynomial Differential Equation algorithm: Parametric Version.

    Explanation
    ===========

    Given a derivation D on k[t], an integer n, and a, b, q1, ..., qm in k[t]
    with deg(a) > 0 and gcd(a, b) == 1, return (A, B, Q, R, n1), with
    Qq = [q1, ..., qm] and R = [r1, ..., rm], such that for any solution
    c1, ..., cm in Const(k) and q in k[t] of degree at most n of
    a*Dq + b*q == Sum(ci*gi, (i, 1, m)), p = (q - Sum(ci*ri, (i, 1, m)))/a has
    degree at most n1 and satisfies A*Dp + B*p == Sum(ci*qi, (i, 1, m))

    This is ``ParSPDE`` from Section 7.1 of Bronstein's book.
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

    Explanation
    ===========

    Given a derivation D on k[t], n in ZZ, and b, q1, ..., qm in k[t] with
    b != 0 and either D == d/dt or deg(b) > max(0, deg(D) - 1), returns
    h1, ..., hr in k[t] and a matrix A with coefficients in Const(k) such that
    if c1, ..., cm in Const(k) and q in k[t] satisfy deg(q) <= n and
    Dq + b*q == Sum(ci*qi, (i, 1, m)), then q = Sum(dj*hj, (j, 1, r)), where
    d1, ..., dr in Const(k) and A*Matrix([[c1, ..., cm, d1, ..., dr]]).T == 0.

    This is ``ParamPolyRischDENoCancel1`` from Section 7.1 of Bronstein's
    book.
    """
    db = b.degree(DE.t)
    m = len(Q)
    H = [Poly(0, DE.t)]*m

    for N, i in itertools.product(range(n, -1, -1), range(m)):  # [n, ..., 0]
        si = Q[i].nth(N + db)/b.LC()
        sitn = Poly(si*DE.t**N, DE.t)
        H[i] = H[i] + sitn
        Q[i] = Q[i] - derivation(sitn, DE) - b*sitn

    if all(qi.is_zero for qi in Q):
        dc = -1
    else:
        dc = max(qi.degree(DE.t) for qi in Q)
    M = Matrix(dc + 1, m, lambda i, j: Q[j].nth(i), DE.t)
    A, u = constant_system(M, zeros(dc + 1, 1, DE.t), DE)
    c = eye(m, DE.t)
    A = A.row_join(zeros(A.rows, m, DE.t)).col_join(c.row_join(-c))

    return (H, A)


def prde_no_cancel_b_small(b, Q, n, DE):
    """
    Parametric Poly Risch Differential Equation - No cancellation: deg(b) small enough.

    Explanation
    ===========

    Given a derivation D on k[t], n in ZZ, and b, q1, ..., qm in k[t] with
    deg(b) < deg(D) - 1 and either D == d/dt or deg(D) >= 2, returns
    h1, ..., hr in k[t] and a matrix A with coefficients in Const(k) such that
    if c1, ..., cm in Const(k) and q in k[t] satisfy deg(q) <= n and
    Dq + b*q == Sum(ci*qi, (i, 1, m)) then q = Sum(dj*hj, (j, 1, r)) where
    d1, ..., dr in Const(k) and A*Matrix([[c1, ..., cm, d1, ..., dr]]).T == 0.

    This is ``ParamPolyRischDENoCancel2`` from Section 7.1 of Bronstein's
    book.
    """
    m = len(Q)
    H = [Poly(0, DE.t)]*m

    for N, i in itertools.product(range(n, 0, -1), range(m)):  # [n, ..., 1]
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
        else:
            dc = max(qi.degree(DE.t) for qi in Q)
        M = Matrix(dc + 1, m, lambda i, j: Q[j].nth(i), DE.t)
        A, u = constant_system(M, zeros(dc + 1, 1, DE.t), DE)
        c = eye(m, DE.t)
        A = A.row_join(zeros(A.rows, m, DE.t)).col_join(c.row_join(-c))
        return (H, A)

    # else: b is in k, deg(qi) < deg(Dt)

    t = DE.t
    if DE.case != 'base':
        with DecrementLevel(DE):
            t0 = DE.t  # k = k0(t0)
            ba, bd = frac_in(b, t0, field=True)
            Q0 = [frac_in(qi.TC(), t0, field=True) for qi in Q]
            f, B = param_rischDE(ba, bd, Q0, DE)

            # f = [f1, ..., fr] in k^r and B is a matrix with
            # m + r columns and entries in Const(k) = Const(k0)
            # such that Dy0 + b*y0 = Sum(ci*qi, (i, 1, m)) has
            # a solution y0 in k with c1, ..., cm in Const(k)
            # if and only y0 = Sum(dj*fj, (j, 1, r)) where
            # d1, ..., dr ar in Const(k) and
            # B*Matrix([c1, ..., cm, d1, ..., dr]) == 0.

        # Transform fractions (fa, fd) in f into constant
        # polynomials fa/fd in k[t].
        # (Is there a better way?)
        f = [Poly(fa.as_expr()/fd.as_expr(), t, field=True)
             for fa, fd in f]
        B = Matrix.from_Matrix(B.to_Matrix(), t)
    else:
        # Base case. Dy == 0 for all y in k and b == 0.
        # Dy + b*y = Sum(ci*qi) is solvable if and only if
        # Sum(ci*qi) == 0 in which case the solutions are
        # y = d1*f1 for f1 = 1 and any d1 in Const(k) = k.

        f = [Poly(1, t, field=True)]  # r = 1
        B = Matrix([[qi.TC() for qi in Q] + [S.Zero]], DE.t)
        # The condition for solvability is
        # B*Matrix([c1, ..., cm, d1]) == 0
        # There are no constraints on d1.

    # Coefficients of t^j (j > 0) in Sum(ci*qi) must be zero.
    d = max(qi.degree(DE.t) for qi in Q)
    if d > 0:
        M = Matrix(d, m, lambda i, j: Q[j].nth(i + 1), DE.t)
        A, _ = constant_system(M, zeros(d, 1, DE.t), DE)
    else:
        # No constraints on the hj.
        A = Matrix(0, m, [], DE.t)

    # Solutions of the original equation are
    #    y = Sum(dj*fj, (j, 1, r) + Sum(ei*hi, (i, 1, m)),
    # where  ei == ci  (i = 1, ..., m),  when
    # A*Matrix([c1, ..., cm]) == 0 and
    # B*Matrix([c1, ..., cm, d1, ..., dr]) == 0

    # Build combined constraint matrix with m + r + m columns.

    r = len(f)
    I = eye(m, DE.t)
    A = A.row_join(zeros(A.rows, r + m, DE.t))
    B = B.row_join(zeros(B.rows, m, DE.t))
    C = I.row_join(zeros(m, r, DE.t)).row_join(-I)

    return f + H, A.col_join(B).col_join(C)


def prde_no_cancel_b_equal(b, Q, n, DE):
    """
    Parametric Poly Risch Differential Equation - No cancellation: deg(b) == delta(t) - 1

    Explanation
    ===========

    Given a derivation D on k[t] with delta(t) >= 2, n in ZZ, and
    b, q1, ..., qm in k[t] with deg(b) == delta(t) - 1 and
    n > -lc(b)/lc(Dt), returns h1, ..., hr in k[t] and a matrix A with
    coefficients in Const(k) such that if c1, ..., cm in Const(k) and
    q in k[t] satisfy deg(q) <= n and Dq + b*q == Sum(ci*qi, (i, 1, m)),
    then q == Sum(dj*hj, (j, 1, r)), where d1, ..., dr in Const(k) and
    A*Matrix([c1, ..., cm, d1, ..., dr]).T == 0.

    This implements the "When delta(t) >= 2 and deg(b) == delta(t) - 1"
    discussion of Section 7.1 of Bronstein's book (neither edition gives
    it as named pseudocode).  If the possible-cancellation degree
    -lc(b)/lc(Dt) is a positive integer at most n, the remaining problem
    is delegated to the cancellation algorithms via param_poly_rischDE(),
    which for delta(t) >= 2 are not yet implemented, so this case
    currently raises NotImplementedError.
    """
    m = len(Q)
    delta = DE.d.degree(DE.t)
    lam = DE.d.LC()

    # The degree at which cancellation could occur.
    # TODO: If b or Dt contains symbolic constants, Mc may be an integer
    # only for special values of them (e.g. b == -y*t gives Mc == y), and
    # the divisions by u below then produce a result that is valid only
    # generically: e.g. H == [-t**3/(y - 3) - 3*t/(y**2 - 4*y + 3)] is
    # undefined at y == 3 and y == 1.  The correct result would be a
    # Piecewise with one condition per loop iteration where
    # solve(N*lam + lc(b), y) is nonempty, with the cancellation
    # algorithms handling those special values.
    Mc = cancel(-b.LC()/lam)
    if Mc.is_Integer and Mc > 0:
        M = int(Mc)
    else:
        M = -1

    H = [Poly(0, DE.t)]*m

    for N in range(n, M, -1):  # [n, ..., M + 1]
        u = cancel(N*lam + b.LC())  # nonzero, since N != -lc(b)/lam
        for i in range(m):
            si = Q[i].nth(N + delta - 1)/u
            sitn = Poly(si*DE.t**N, DE.t)
            H[i] = H[i] + sitn
            Q[i] = Q[i] - derivation(sitn, DE) - b*sitn

    if M == -1:
        # The loop ran through N == 0, so as in prde_no_cancel_b_large()
        # any solution is q == Sum(ci*hi) and the (updated) right hand
        # sides must satisfy Sum(ci*qi) == 0.
        if all(qi.is_zero for qi in Q):
            dc = -1
        else:
            dc = max(qi.degree(DE.t) for qi in Q)
        Mmat = Matrix(dc + 1, m, lambda i, j: Q[j].nth(i), DE.t)
        A, _ = constant_system(Mmat, zeros(dc + 1, 1, DE.t), DE)
        c = eye(m, DE.t)
        A = A.row_join(zeros(A.rows, m, DE.t)).col_join(c.row_join(-c))
        return (H, A)

    # We reached the possible-cancellation degree M > 0.  Solve the
    # remaining problem, with the updated right hand sides and degree
    # bound M, by the cancellation algorithms.
    f, B = param_poly_rischDE(Poly(1, DE.t), b, Q, M, DE)
    # Any solution is q == Sum(ci*hi) + Sum(dj*fj), where the fj and the
    # matrix B constrain (c1, ..., cm, d1, ..., dr).  Combine into a
    # relation matrix over the columns
    # (c1, ..., cm, d_h1, ..., d_hm, d_f1, ..., d_fr), where the first
    # block of rows enforces d_hi == ci.
    r = len(f)
    Bc = Matrix(B.rows, m, lambda i, j: B[i, j], DE.t)
    Bd = Matrix(B.rows, r, lambda i, j: B[i, m + j], DE.t)
    top = eye(m, DE.t).row_join(-eye(m, DE.t)).row_join(zeros(m, r, DE.t))
    bottom = Bc.row_join(zeros(B.rows, m, DE.t)).row_join(Bd)
    return (H + f, top.col_join(bottom))


def prde_cancel_liouvillian(b, Q, n, DE):
    """
    Parametric Poly Risch Differential Equation - Cancellation: Liouvillian case.

    This implements the primitive and hyperexponential cancellation cases
    from the discussion in Section 7.1 of Bronstein's book (neither edition
    gives it as named pseudocode).
    """
    if DE.case not in ('primitive', 'exp'):
        raise ValueError("case must be 'primitive' or 'exp', not %s." %
            DE.case)

    H = []

    if DE.case == 'exp':
        # The coefficient of t**i satisfies the equation
        # D(q_i) + (b + i*eta)*q_i == rhs_i over k, where eta == Dt/t.
        # eta must be computed at the current level, before DecrementLevel
        # is entered below.
        eta = DE.d.quo(Poly(DE.t, DE.t)).as_expr()

    # Why use DecrementLevel? Below line answers that:
    # Assuming that we can solve such problems over 'k' (not k[t])
    if DE.case == 'primitive':
        with DecrementLevel(DE):
            ba, bd = frac_in(b, DE.t, field=True)

    for i in range(n, -1, -1):
        if DE.case == 'exp': # this re-checking can be avoided
            with DecrementLevel(DE):
                ba, bd = frac_in(b.as_expr() + i*eta, DE.t, field=True)
        with DecrementLevel(DE):
            Qy = [frac_in(q.nth(i), DE.t, field=True) for q in Q]
            fi, Ai = param_rischDE(ba, bd, Qy, DE)
        fi = [Poly(fa.as_expr()/fd.as_expr(), DE.t, field=True)
                for fa, fd in fi]
        Ai = Ai.set_gens(DE.t)

        ri = len(fi)

        if i == n:
            M = Ai
        else:
            M = Ai.col_join(M.row_join(zeros(M.rows, ri, DE.t)))

        Fi, hi = [None]*ri, [None]*ri

        # Substituting q == d*h + q_rest into Dq + b*q == Sum(ci*qi)
        # leaves the residual equation
        # D(q_rest) + b*q_rest == Sum(ci*qi) - d*(D(h) + b*h)
        # (see the equation following (7.17) in Section 7.1).
        for j in range(ri):
            hji = fi[j] * (DE.t**i).as_poly(fi[j].gens)
            hi[j] = hji
            # building up Sum(dji*(D(fji*t^i) + b*fji*t^i))
            Fi[j] = -(derivation(hji, DE) + b*hji)

        H += hi
        # in the next loop instead of Q it has
        # to be Q + Fi taking its place
        Q = Q + Fi

    return (H, M)


def param_poly_rischDE(a, b, q, n, DE):
    """Polynomial solutions of a parametric Risch differential equation.

    Explanation
    ===========

    Given a derivation D in k[t], a, b in k[t] relatively prime, and q
    = [q1, ..., qm] in k[t]^m, return h = [h1, ..., hr] in k[t]^r and
    a matrix A with m + r columns and entries in Const(k) such that
    a*Dp + b*p = Sum(ci*qi, (i, 1, m)) has a solution p of degree <= n
    in k[t] with c1, ..., cm in Const(k) if and only if p = Sum(dj*hj,
    (j, 1, r)) where d1, ..., dr are in Const(k) and (c1, ..., cm,
    d1, ..., dr) is a solution of Ax == 0.

    This implements the polynomial parametric RDE procedure of Section 7.1
    of Bronstein's book; there is no single corresponding pseudocode
    function.
    """
    m = len(q)
    if n < 0:
        # Only the trivial zero solution is possible.
        # Find relations between the qi.
        if all(qi.is_zero for qi in q):
            return [], zeros(1, m, DE.t)  # No constraints.

        N = max(qi.degree(DE.t) for qi in q)
        M = Matrix(N + 1, m, lambda i, j: q[j].nth(i), DE.t)
        A, _ = constant_system(M, zeros(M.rows, 1, DE.t), DE)

        return [], A

    if a.is_ground:
        # Normalization: a = 1.
        a = a.LC()
        b, q = b.to_field().exquo_ground(a), [qi.to_field().exquo_ground(a) for qi in q]

        if not b.is_zero and (DE.case == 'base' or
                b.degree() > max(0, DE.d.degree() - 1)):
            return prde_no_cancel_b_large(b, q, n, DE)

        elif ((b.is_zero or b.degree() < DE.d.degree() - 1)
                and (DE.case == 'base' or DE.d.degree() >= 2)):
            return prde_no_cancel_b_small(b, q, n, DE)

        elif (DE.d.degree() >= 2 and
              b.degree() == DE.d.degree() - 1 and
              n > -b.as_poly().LC()/DE.d.as_poly().LC()):
            return prde_no_cancel_b_equal(b, q, n, DE)

        else:
            # Liouvillian cases
            if DE.case in ('primitive', 'exp'):
                return prde_cancel_liouvillian(b, q, n, DE)
            else:
                raise NotImplementedError("non-linear and hypertangent "
                        "cases have not yet been implemented")

    # else: deg(a) > 0

    # Iterate SPDE as long as possible cumulating coefficient
    # and terms for the recovery of original solutions.
    alpha, beta = a.one, [a.zero]*m
    while n >= 0:  # and a, b relatively prime
        a, b, q, r, n = prde_spde(a, b, q, n, DE)
        beta = [betai + alpha*ri for betai, ri in zip(beta, r)]
        alpha *= a
        # Solutions p of a*Dp + b*p = Sum(ci*qi) correspond to
        # solutions alpha*p + Sum(ci*betai) of the initial equation.
        d = a.gcd(b)
        if not d.is_ground:
            break

    # a*Dp + b*p = Sum(ci*qi) may have a polynomial solution
    # only if the sum is divisible by d.

    qq, M = poly_linear_constraints(q, d)
    # qq = [qq1, ..., qqm] where qqi = qi.quo(d).
    # M is a matrix with m columns an entries in k.
    # Sum(fi*qi, (i, 1, m)), where f1, ..., fm are elements of k, is
    # divisible by d if and only if M*Matrix([f1, ..., fm]) == 0,
    # in which case the quotient is Sum(fi*qqi).

    A, _ = constant_system(M, zeros(M.rows, 1, DE.t), DE)
    # A is a matrix with m columns and entries in Const(k).
    # Sum(ci*qqi) is Sum(ci*qi).quo(d), and the remainder is zero
    # for c1, ..., cm in Const(k) if and only if
    # A*Matrix([c1, ...,cm]) == 0.

    V = A.nullspace()
    # V = [v1, ..., vu] where each vj is a column matrix with
    # entries aj1, ..., ajm in Const(k).
    # Sum(aji*qi) is divisible by d with exact quotient Sum(aji*qqi).
    # Sum(ci*qi) is divisible by d if and only if ci = Sum(dj*aji)
    # (i = 1, ..., m) for some d1, ..., du in Const(k).
    # In that case, solutions of
    #     a*Dp + b*p = Sum(ci*qi) = Sum(dj*Sum(aji*qi))
    # are the same as those of
    #     (a/d)*Dp + (b/d)*p = Sum(dj*rj)
    # where rj = Sum(aji*qqi).

    if not V:  # No non-trivial solution.
        return [], eye(m, DE.t)  # Could return A, but this has
                                 # the minimum number of rows.

    Mqq = Matrix([qq])  # A single row.
    r = [(Mqq*vj)[0] for vj in V]  # [r1, ..., ru]

    # Solutions of (a/d)*Dp + (b/d)*p = Sum(dj*rj) correspond to
    # solutions alpha*p + Sum(Sum(dj*aji)*betai) of the initial
    # equation. These are equal to alpha*p + Sum(dj*fj) where
    # fj = Sum(aji*betai).
    Mbeta = Matrix([beta])
    f = [(Mbeta*vj)[0] for vj in V]  # [f1, ..., fu]

    #
    # Solve the reduced equation recursively.
    #
    g, B = param_poly_rischDE(a.quo(d), b.quo(d), r, n, DE)

    # g = [g1, ..., gv] in k[t]^v and and B is a matrix with u + v
    # columns and entries in Const(k) such that
    # (a/d)*Dp + (b/d)*p = Sum(dj*rj) has a solution p of degree <= n
    # in k[t] if and only if p = Sum(ek*gk) where e1, ..., ev are in
    # Const(k) and B*Matrix([d1, ..., du, e1, ..., ev]) == 0.
    # The solutions of the original equation are then
    # Sum(dj*fj, (j, 1, u)) + alpha*Sum(ek*gk, (k, 1, v)).

    # Collect solution components.
    h = f + [alpha*gk for gk in g]

    # Build combined relation matrix.
    A = -eye(m, DE.t)
    for vj in V:
        A = A.row_join(vj)
    A = A.row_join(zeros(m, len(g), DE.t))
    A = A.col_join(zeros(B.rows, m, DE.t).row_join(B))

    return h, A


def _prde_normalized_solve(A, B, G, gamma, DE, n=None):
    """
    Common tail of param_rischDE() and limited_integrate().

    Explanation
    ===========

    Given A, B in k[t] and G = [G1, ..., Gm] in k(t)^m with the equation
    A*Dp + B*p == Sum(ci*Gi, (i, 1, m)) already reduced (weakly
    normalized, denominators cleared), return h = [h1, ..., hr] in
    k(t)^r and a matrix C with entries in Const(k) such that solutions
    of the original problem are y == Sum(dj*hj, (j, 1, r))/1 with
    (c1, ..., cm, d1, ..., dr) in the nullspace of C, where the
    returned hj have already been divided by gamma.  If ``n`` is None,
    the degree bound is computed with bound_degree(); otherwise ``n``
    is used as a proven degree bound, avoiding bound_degree() entirely.
    """
    m = len(G)
    g = A.gcd(B)
    a, b, g = A.quo(g), B.quo(g), [gia.cancel(gid*g, include=True) for
        gia, gid in G]

    # a*Dp + b*p = Sum(ci*gi)  may have a polynomial solution
    # only if the sum is in k[t].

    q, M = prde_linear_constraints(a, b, g, DE)

    # q = [q1, ..., qm] where qi in k[t] is the polynomial component
    # of the partial fraction expansion of gi.
    # M is a matrix with m columns and entries in k.
    # Sum(fi*gi, (i, 1, m)), where f1, ..., fm are elements of k,
    # is a polynomial if and only if M*Matrix([f1, ..., fm]) == 0,
    # in which case the sum is equal to Sum(fi*qi).

    M, _ = constant_system(M, zeros(M.rows, 1, DE.t), DE)
    # M is a matrix with m columns and entries in Const(k).
    # Sum(ci*gi) is in k[t] for c1, ..., cm in Const(k)
    # if and only if M*Matrix([c1, ..., cm]) == 0,
    # in which case the sum is Sum(ci*qi).

    ## Reduce number of constants at this point

    V = M.nullspace()
    # V = [v1, ..., vu] where each vj is a column matrix with
    # entries aj1, ..., ajm in Const(k).
    # Sum(aji*gi) is in k[t] and equal to Sum(aji*qi) (j = 1, ..., u).
    # Sum(ci*gi) is in k[t] if and only is ci = Sum(dj*aji)
    # (i = 1, ..., m) for some d1, ..., du in Const(k).
    # In that case,
    #     Sum(ci*gi) = Sum(ci*qi) = Sum(dj*Sum(aji*qi)) = Sum(dj*rj)
    # where rj = Sum(aji*qi) (j = 1, ..., u) in k[t].

    if not V:  # No non-trivial solution
        return [], eye(m, DE.t)

    Mq = Matrix([q])  # A single row.
    r = [(Mq*vj)[0] for vj in V]  # [r1, ..., ru]

    # Solutions of a*Dp + b*p = Sum(dj*rj) correspond to solutions
    # y = p/gamma of the initial equation with ci = Sum(dj*aji).

    if n is None:
        # bound_degree() raises NotImplementedError when it cannot decide
        # one of its structure-theorem queries, and we let that propagate.
        # Guessing a finite bound here instead (n = 5 from 2017 to 2026)
        # is unsound: solutions of degree above the guess are silently
        # lost, so a caller like is_deriv_in_field() can be tricked into
        # a false proof that an elementary integral is nonelementary.
        # Falling back to n == oo (as rischDE() does) is not an option
        # either, since the parametric no-cancellation and cancellation
        # algorithms iterate over range(n, ...), which requires a finite
        # bound.
        n = bound_degree(a, b, r, DE, parametric=True)

    h, B = param_poly_rischDE(a, b, r, n, DE)

    # h = [h1, ..., hv] in k[t]^v and and B is a matrix with u + v
    # columns and entries in Const(k) such that
    # a*Dp + b*p = Sum(dj*rj) has a solution p of degree <= n
    # in k[t] if and only if p = Sum(ek*hk) where e1, ..., ev are in
    # Const(k) and B*Matrix([d1, ..., du, e1, ..., ev]) == 0.
    # The solutions of the original equation for ci = Sum(dj*aji)
    # (i = 1, ..., m) are then y = Sum(ek*hk, (k, 1, v))/gamma.

    ## Build combined relation matrix with m + u + v columns.

    A = -eye(m, DE.t)
    for vj in V:
        A = A.row_join(vj)
    A = A.row_join(zeros(m, len(h), DE.t))
    A = A.col_join(zeros(B.rows, m, DE.t).row_join(B))

    ## Eliminate d1, ..., du.

    W = A.nullspace()

    # W = [w1, ..., wt] where each wl is a column matrix with
    # entries blk (k = 1, ..., m + u + v) in Const(k).
    # The vectors (bl1, ..., blm) generate the space of those
    # constant families (c1, ..., cm) for which a solution of
    # the equation Dy + f*y == Sum(ci*Gi) exists. They generate
    # the space and form a basis except possibly when Dy + f*y == 0
    # is solvable in k(t}. The corresponding solutions are
    # y = Sum(blk'*hk, (k, 1, v))/gamma, where k' = k + m + u.

    v = len(h)
    shape = (len(W), m+v)
    elements = [wl[:m] + wl[-v:] for wl in W] # excise dj's.
    items = [e for row in elements for e in row]

    # Need to set the shape in case W is empty
    M = Matrix(*shape, items, DE.t)
    N = M.nullspace()

    # N = [n1, ..., ns] where the ni in Const(k)^(m + v) are column
    # vectors generating the space of linear relations between
    # c1, ..., cm, e1, ..., ev.

    C = Matrix([ni[:] for ni in N], DE.t)  # rows n1, ..., ns.

    return [hk.cancel(gamma, include=True) for hk in h], C


def param_rischDE(fa, fd, G, DE):
    """
    Solve a Parametric Risch Differential Equation: Dy + f*y == Sum(ci*Gi, (i, 1, m)).

    Explanation
    ===========

    Given a derivation D in k(t), f in k(t), and G
    = [G1, ..., Gm] in k(t)^m, return h = [h1, ..., hr] in k(t)^r and
    a matrix A with m + r columns and entries in Const(k) such that
    Dy + f*y = Sum(ci*Gi, (i, 1, m)) has a solution y
    in k(t) with c1, ..., cm in Const(k) if and only if y = Sum(dj*hj,
    (j, 1, r)) where d1, ..., dr are in Const(k) and (c1, ..., cm,
    d1, ..., dr) is a solution of Ax == 0.

    Elements of k(t) are tuples (a, d) with a and d in k[t].

    This chains together the parametric algorithms of Section 7.1 of
    Bronstein's book; the book does not give it as a single named
    pseudocode function.
    """
    q, (fa, fd) = weak_normalizer(fa, fd, DE)
    # Solutions of the weakly normalized equation Dz + f*z = q*Sum(ci*Gi)
    # correspond to solutions y = z/q of the original equation.
    gamma = q
    G = [(q*ga).cancel(gd, include=True) for ga, gd in G]

    a, (ba, bd), G, hn = prde_normal_denom(fa, fd, G, DE)
    # Solutions q in k<t> of  a*Dq + b*q = Sum(ci*Gi) correspond
    # to solutions z = q/hn of the weakly normalized equation.
    gamma *= hn

    A, B, G, hs = prde_special_denom(a, ba, bd, G, DE)
    # Solutions p in k[t] of  A*Dp + B*p = Sum(ci*Gi) correspond
    # to solutions q = p/hs of the previous equation.
    gamma *= hs

    return _prde_normalized_solve(A, B, G, gamma, DE)


def limited_integrate_reduce(fa, fd, G, DE):
    """
    Simpler version of step 1 & 2 for the limited integration problem.

    Explanation
    ===========

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

    This is ``LimitedIntegrateReduce`` from Section 7.2 of Bronstein's
    book.
    """
    dn, ds = splitfactor(fd, DE)
    E = [splitfactor(gd, DE) for _, gd in G]
    En, Es = list(zip(*E)) if E else ((), ())
    c = reduce(lambda i, j: i.lcm(j), (dn,) + En)  # lcm(dn, en1, ..., enm)
    hn = c.gcd(c.diff(DE.t))
    a = hn
    b = -derivation(hn, DE)
    N = 0

    # These are the cases where we know that S1irr = Sirr, but there could be
    # others, and this algorithm will need to be extended to handle them.
    if DE.case in ('base', 'primitive', 'exp', 'tan'):
        hs = reduce(lambda i, j: i.lcm(j), (ds,) + Es)  # lcm(ds, es1, ..., esm)
        a = hn*hs
        b -= (hn*derivation(hs, DE)).quo(hs)
        mu = min([order_at_oo(fa, fd, DE.t)] + [order_at_oo(ga, gd, DE.t) for
            ga, gd in G])
        # So far, all the above are also nonlinear or Liouvillian, but if this
        # changes, then this will need to be updated to call bound_degree()
        # as per the docstring of this function (DE.case == 'other_linear').
        N = hn.degree(DE.t) + hs.degree(DE.t) + max(0, 1 - DE.d.degree(DE.t) - mu)
    else:
        # TODO: implement this
        raise NotImplementedError

    V = [(-a*hn*ga).cancel(gd, include=True) for ga, gd in G]
    # Note: the first component is hn, not a.  Both editions of the book
    # return a here, but that does not satisfy the stated contract when
    # hs != 1: from v == p/a, the original equation multiplied by a**2
    # gives a*Dp - Da*p == a**2*f - Sum(ci*a**2*wi), and dividing through
    # by hs (which divides Da exactly, since hs is special) gives
    # hn*Dp + b*p == a*hn*f - Sum(ci*a*hn*wi) with the b and right hand
    # sides constructed above.
    return (hn, b, a, N, (a*hn*fa).cancel(fd, include=True), V)


def limited_integrate(fa, fd, G, DE):
    """
    Solves the limited integration problem:  f = Dv + Sum(ci*wi, (i, 1, n))

    This implements the limited integration procedure of Section 7.2 of
    Bronstein's book (there is no single corresponding pseudocode
    function).
    """
    fa, fd = fa*Poly(1/fd.LC(), DE.t), fd.monic()
    # Reduction to a polynomial problem (LimitedIntegrateReduce, Section
    # 7.2): for any solution v in k(t), c1, ..., cm in Const(k) of
    # f == Dv + Sum(ci*wi), p == v*h is in k[t] and satisfies
    # A*Dp + b*p == g + Sum(ci*vi) with deg(p) <= N.
    A, b, h, N, g, V = limited_integrate_reduce(fa, fd, G, DE)
    # Solve the parametric problem with right hand sides [g] + V.
    # Solutions of the original problem correspond to parametric
    # solutions whose g-coefficient is 1.  N is a proven degree bound,
    # so bound_degree() is not needed.
    Q = [g] + V
    hs, C = _prde_normalized_solve(A, b, Q, h, DE, n=N)
    W = C.nullspace()
    W = [w for w in W if w[0] != 0]
    if not W:
        return None
    # we can take any vector from W; take W[0], scaled so that the
    # g-coefficient is 1
    # TODO: If W[0][0] contains symbolic constants, it can vanish for
    # special values of them, making the result valid only generically;
    # a different nullspace vector may serve those special values, and
    # the correct result would be a Piecewise over the corresponding
    # conditions.
    w = W[0]/W[0][0]
    r = len(hs)
    m = len(w) - r - 1
    C = list(w[1: m + 1])
    y = sum((w[m + 1 + i]*hs[i][0].as_expr()/hs[i][1].as_expr()
            for i in range(r)), S.Zero)
    y_num, y_den = y.as_numer_denom()
    Ya, Yd = Poly(y_num, DE.t), Poly(y_den, DE.t)
    Y = Ya*Poly(1/Yd.LC(), DE.t), Yd.monic()
    return Y, C


def is_deriv_in_field(fa, fd, DE):
    """
    Checks if f can be written as the derivative of an element of k(t).

    Explanation
    ===========

    f in k(t) is the derivative of an element of k(t) if there exists v in
    k(t) such that f == Dv.  Either returns (va, vd), with va, vd in k[t]
    such that f == D(va/vd), or None, which means that f is not the
    derivative of an element of k(t).  Note that v is determined only up
    to an additive constant.

    This is the in-field integration problem from Section 5.12 of
    Bronstein's book (the "Recognizing Derivatives" subsection; neither
    edition gives it as named pseudocode).  It is solved here as the
    special case of the limited integration problem (Section 7.2) with an
    empty list of special elements.

    See also
    ========
    limited_integrate, is_log_deriv_k_t_radical_in_field
    """
    A = limited_integrate(fa, fd, [], DE)
    if A is None:
        return None
    (va, vd), _ = A
    return (va, vd)


def parametric_log_deriv_heu(fa, fd, wa, wd, DE, c1=None):
    """
    Parametric logarithmic derivative heuristic.

    Explanation
    ===========

    Given a derivation D on k[t], f in k(t), and a hyperexponential monomial
    theta over k(t), raises either NotImplementedError, in which case the
    heuristic failed, or returns None, in which case it has proven that no
    solution exists, or returns a solution (n, m, v) of the equation
    n*f == Dv/v + m*Dtheta/theta, with v in k(t)* and n, m in ZZ with n != 0.

    If this heuristic fails, the structure theorem approach will need to be
    used.

    The argument w == Dtheta/theta

    This is ``ParametricLogarithmicDerivative`` from Section 7.3 of
    Bronstein's book.
    """
    # Note: ideally the code here should let n = 1 whenever possible, as the
    # check in special_denom() only succeeds in that case.

    # Special case when f and w are rational numbers (not in the book, but
    # this fails this heuristic and comes up often enough to add here). In
    # this case, we can set z = 0.
    f = fa.as_expr()/fd.as_expr()
    w = wa.as_expr()/wd.as_expr()
    if f.is_Rational and w.is_Rational:
        if f == 0 and w == 0:
            # Any n works: n*0 == Dv/v + m*0 with v == 1, m == 0.
            return (1, 0, Poly(1, DE.t))
        # solve n*x = m*y in integers, i.e., n=y, m=x (after dividing through
        # by gcd(x, y))
        x = f.p*w.q
        y = w.p*f.q
        g = igcd(x, y)
        x //=g
        y //=g
        n = y
        m = x
        if n == 0:
            return None
        if n < 0:
            n = -n
            m = -m
        return (n, m, Poly(1, DE.t))

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

        M, N = s[c1].as_numer_denom()
        M_poly = M.as_poly(q.gens)
        N_poly = N.as_poly(q.gens)

        nfmwa = N_poly*fa*wd - M_poly*wa*fd
        nfmwd = fd*wd
        Qv = is_log_deriv_k_t_radical_in_field(nfmwa, nfmwd, DE, 'auto')
        if Qv is None:
            # (N*f - M*w) is not the logarithmic derivative of a k(t)-radical.
            return None

        Q, v = Qv

        if Q.is_zero or v.is_zero:
            return None

        return (Q*N, Q*M, v)

    if p.degree(DE.t) > B:
        return None

    c = lcm(fd.as_poly(DE.t).LC(), wd.as_poly(DE.t).LC())
    l = fd.monic().lcm(wd.monic())*Poly(c, DE.t)
    ln, ls = splitfactor(l, DE)
    z = ls*ln.gcd(ln.diff(DE.t))

    if z.degree(DE.t) < 1:
        # The residue equations below would be vacuous, so the heuristic
        # cannot determine the candidate for m/n from them ("failed" in
        # Bronstein's ParametricLogarithmicDerivative).  Solutions may well
        # exist, e.g. f == -1 + 1/(x + 1), w == 1 has the solution
        # (1, -1, x + 1), so returning None here (proven no solution) would
        # be wrong.
        if DE.case == 'base':
            # At the base level, Dv/v is a proper fraction for any
            # v in k(x)*, so n*f == Dv/v + m*w requires the polynomial
            # parts to cancel exactly: n*p == m*q.  This determines
            # c1 == m/n == p/q (a constant), and the rest is decided by
            # is_log_deriv_k_t_radical_in_field(), so this case is
            # complete.
            if q.is_zero:
                if not p.is_zero:
                    # n*p == 0 forces n == 0, which is not allowed.
                    return None
                # p == q == 0: m/n is not determined by the polynomial
                # parts; this needs the structure theorems.
                raise NotImplementedError("parametric_log_deriv_heu() "
                    "heuristic failed: the full parametric logarithmic "
                    "derivative problem (using the structure theorems) is "
                    "not yet implemented.")
            if p.is_zero:
                cc = S.Zero
            else:
                if p*Poly(q.LC(), DE.t) != q*Poly(p.LC(), DE.t):
                    # p is not a constant multiple of q, so n*p == m*q has
                    # no solution with n != 0.
                    return None
                cc = cancel(p.LC()/q.LC())
                if not cc.is_Rational:
                    # If cc contains symbolic constants, it may be
                    # rational only for special values of them (a
                    # condition like "y in QQ" that is not expressible
                    # as a Piecewise relational); None here is correct
                    # generically and at worst yields a false
                    # nonelementary proof at those values, never a
                    # wrong antiderivative, so no Piecewise handling is
                    # needed.
                    return None
            M, N = cc.as_numer_denom()

            nfmwa = Poly(N, DE.t)*fa*wd - Poly(M, DE.t)*wa*fd
            nfmwd = fd*wd
            Qv = is_log_deriv_k_t_radical_in_field(nfmwa, nfmwd, DE)
            if Qv is None:
                # (N*f - M*w) is not the logarithmic derivative of a
                # k(t)-radical.
                return None

            Q, v = Qv

            if Q.is_zero or v.is_zero:
                return None

            return (Q*N, Q*M, v)

        raise NotImplementedError("parametric_log_deriv_heu() heuristic "
            "failed: the full parametric logarithmic derivative problem "
            "(using the structure theorems) is not yet implemented.")

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


def parametric_log_deriv_structure(fa, fd, wa, wd, DE):
    """
    Parametric logarithmic derivative problem via the structure theorems.

    Explanation
    ===========

    Tries to solve n*f == Dv/v + m*w for n, m in ZZ, n != 0, and v in
    k(t)*, by solving the structure equation

        f == (m/n)*w + Sum(ri*Dti/ti, i in E) + Sum(ri*Dti, i in L)

    for rational (m/n, r1, ..., rk), where E and L are the
    hyperexponential and logarithmic monomials of the current tower (up
    to the current level), as in equation (7.44) of Section 7.3 of
    Bronstein's book, with theta adjoined as an extra hyperexponential
    generator.  A solution yields (n, m, v) with
    v == Product(termi**(n*ri)) directly.

    Either returns (n, m, v), or None, which means that f - (m/n)*w is
    not a QQ-linear combination of the logarithmic derivatives available
    from the current tower for any m/n in QQ.  Note that None does NOT
    prove that no solution exists: a solution whose v requires monomials
    not in the tower (e.g. new logarithms) is not found by this method.
    Callers must treat None as inconclusive.
    """
    # The same assumptions as in is_log_deriv_k_t_radical()
    if len(DE.exts) != len(DE.D) - 1:
        if [i for i in DE.cases if i == 'tan'] or \
                ({i for i, case in enumerate(DE.cases) if case == 'primitive'} -
                        set(DE.indices('log'))):
            raise NotImplementedError("Real version of the structure "
                "theorems with hypertangent support is not yet implemented.")
        raise NotImplementedError("Nonelementary extensions not supported "
            "in the structure theorems.")

    # Use only the part of the tower visible at the current level: the
    # callers pose this problem inside DecrementLevel, and v must lie in
    # the current field (w itself is usually the log derivative of the
    # monomial one level up, which must not appear in v).
    top = len(DE.T) + DE.level
    E_indices = [i for i in DE.indices('exp') if i <= top]
    L_indices = [i for i in DE.indices('log') if i <= top]
    E_part = [DE.D[i].quo(Poly(DE.T[i], DE.T[i])).as_expr() for i in E_indices]
    L_part = [DE.D[i].as_expr() for i in L_indices]

    w = wa.as_expr()/wd.as_expr()
    f = fa.as_expr()/fd.as_expr()

    dum = Dummy()
    lhs = Matrix([[w] + E_part + L_part], dum)
    rhs = Matrix([f], dum)

    A, u = constant_system(lhs, rhs, DE)

    # A*x == u is now a linear system with constant coefficients for the
    # rational unknowns x == (m/n, r1, ..., rk).  Note that u is the
    # reduced right hand side, not a solution: the system may be
    # underdetermined (e.g. when w is a combination of the generators),
    # so solve it explicitly, setting any free parameters to zero.
    Am = A.to_Matrix()
    um = u.to_Matrix()
    if not (all(i.is_Rational for i in Am) and
            all(i.is_Rational for i in um)):
        if not all(derivation(i, DE, basic=True).is_zero for i in Am.vec()) \
                or not all(derivation(i, DE, basic=True).is_zero for i in um):
            # The system could not be reduced to one over the constants
            return None
        raise NotImplementedError("Cannot work with non-rational "
            "coefficients in this case.")
    from sympy.matrices import Matrix as _Matrix
    try:
        xs, params = _Matrix(Am).gauss_jordan_solve(_Matrix(um))
    except ValueError:
        # Inconsistent system: f - (m/n)*w is not a combination of the
        # available logarithmic derivatives for any m/n.
        # With symbolic constants, inconsistency is decided only
        # generically; at special values of the constants a solution may
        # exist.  The wrapper then raises NotImplementedError, which is
        # sound (an honest error, never a wrong answer), so no Piecewise
        # handling is needed here.
        return None
    xs = xs.subs(dict.fromkeys(params, S.Zero))

    n = S.One*reduce(ilcm, [i.as_numer_denom()[1] for i in xs], S.One)
    xs *= n  # now integers: [m, n*r1, ..., n*rk]
    m = xs[0]
    terms = ([DE.T[i] for i in E_indices] + [DE.extargs[i - 1] for i in L_indices])
    v = Mul(*[Pow(i, j) for i, j in zip(terms, xs[1:])])
    if n == 0:
        return None
    # Defensive check of n*f == Dv/v + m*w (Dv/v is the same QQ-linear
    # combination of the generator logarithmic derivatives by construction)
    lhs_check = Add(*[Mul(j, i) for i, j in zip(E_part + L_part, xs[1:])])
    if cancel(n*f - m*w - lhs_check) != 0:
        return None
    if n < 0:
        n, m, v = -n, -m, 1/v
    return (n, m, v)


def parametric_log_deriv(fa, fd, wa, wd, DE):
    """
    Solves the parametric logarithmic derivative problem.

    This is the parametric logarithmic derivative problem from Section 7.3
    of Bronstein's book: n*f == Dv/v + m*w for n, m in ZZ, n != 0, and
    v in k(t)*.  The heuristic (``ParametricLogarithmicDerivative``) is
    tried first; if it cannot decide, the structure theorem method of the
    same section is tried.
    """
    try:
        A = parametric_log_deriv_heu(fa, fd, wa, wd, DE)
    except NotImplementedError:
        A = parametric_log_deriv_structure(fa, fd, wa, wd, DE)
        if A is None:
            # The structure method proves nonexistence only when the
            # elementary extension containing Integral(f) needs no
            # monomials outside the current tower, which we cannot check
            # here, so an inconclusive result must not be reported as
            # "no solution".
            raise NotImplementedError("parametric_log_deriv() could not "
                "decide: the heuristic failed and f - (m/n)*w is not a "
                "combination of logarithmic derivatives from the current "
                "tower for any m/n in QQ.")
    return A


def _structure_system_solve(lhs, rhs, DE):
    """
    Find a constant solution x of the structure system lhs*x == rhs.

    Explanation
    ===========

    lhs and rhs are matrices over K.  Returns a list of constant
    solutions (free parameters set to zero), or None if there is no
    constant solution.  constant_system() returns a reduced *system*,
    not a solution, and treating its right hand side as a solution (as
    this module once did) is only accidentally correct when the
    reduction is identity-like, so the reduced system is solved
    explicitly here and the solution is verified against the original
    equation.

    Raises NotImplementedError if the verification fails, which
    indicates constants not computable by cancel() (see the comment in
    constant_system()).
    """
    A, u = constant_system(lhs, rhs, DE)
    if not A:
        return None
    um = u.to_Matrix()
    if not all(derivation(i, DE, basic=True).is_zero for i in um):
        # No constant solution
        return None
    from sympy.matrices import Matrix as EMatrix
    try:
        xs, params = EMatrix(A.to_Matrix()).gauss_jordan_solve(EMatrix(um))
    except ValueError:
        # With symbolic constants in the system, inconsistency (and
        # the pivoting above) is decided only generically; at special
        # values of the constants a solution may exist.  Returning None
        # here errs on the incomplete-but-sound side (a spurious
        # "not a derivative"/"not a radical" answer, at worst a false
        # nonelementary proof, never a wrong antiderivative), so no
        # Piecewise handling is needed here.
        return None
    xs = xs.subs(dict.fromkeys(params, S.Zero))
    # Defensive verification against the original system
    lhs_m = lhs.to_Matrix()
    rhs_m = rhs.to_Matrix()
    for i in range(lhs_m.rows):
        if cancel(Add(*[lhs_m[i, j]*xs[j] for j in range(lhs_m.cols)])
                - rhs_m[i]) != 0:
            raise NotImplementedError("The candidate solution of the "
                "structure system failed verification; the constant field "
                "may not be computable by cancel().")
    return list(xs)


def is_deriv_k(fa, fd, DE):
    r"""
    Checks if Df/f is the derivative of an element of k(t).

    Explanation
    ===========

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

    This is an application of the Risch structure theorems from Section 9.3
    of Bronstein's book (neither edition gives it as named pseudocode).

    See also
    ========
    is_log_deriv_k_t_radical_in_field, is_log_deriv_k_t_radical
    """
    # Compute Df/f
    dfa, dfd = (fd*derivation(fa, DE) - fa*derivation(fd, DE)), fd*fa
    dfa, dfd = dfa.cancel(dfd, include=True)

    # Our assumption here is that each monomial is recursively transcendental
    if len(DE.exts) != len(DE.D) - 1:
        if [i for i in DE.cases if i == 'tan'] or \
                ({i for i, case in enumerate(DE.cases) if case == 'primitive'} -
                        set(DE.indices('log'))):
            raise NotImplementedError("Real version of the structure "
                "theorems with hypertangent support is not yet implemented.")

        # TODO: What should really be done in this case?
        raise NotImplementedError("Nonelementary extensions not supported "
            "in the structure theorems.")

    E_part = [DE.D[i].quo(Poly(DE.T[i], DE.T[i])).as_expr() for i in DE.indices('exp')]
    L_part = [DE.D[i].as_expr() for i in DE.indices('log')]

    # The expression dfa/dfd might not be polynomial in any of its symbols so we
    # use a Dummy as the generator for PolyMatrix.
    dum = Dummy()
    lhs = Matrix([E_part + L_part], dum)
    rhs = Matrix([dfa.as_expr()/dfd.as_expr()], dum)

    u = _structure_system_solve(lhs, rhs, DE)

    if u is None:
        # No constant solution
        return None
    else:
        if not all(i.is_Rational for i in u):
            raise NotImplementedError("Cannot work with non-rational "
                "coefficients in this case.")
        else:
            terms = ([DE.extargs[i - 1] for i in DE.indices('exp')] +
                    [DE.T[i] for i in DE.indices('log')])
            ans = list(zip(terms, u))
            result = Add(*[Mul(i, j) for i, j in ans])
            argterms = ([DE.T[i] for i in DE.indices('exp')] +
                    [DE.extargs[i - 1] for i in DE.indices('log')])
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
    r"""
    Checks if Df is the logarithmic derivative of a k(t)-radical.

    Explanation
    ===========

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

    This is an application of the structure theorems from Sections 9.3 and
    9.4 of Bronstein's book (neither edition gives it as named pseudocode).

    See also
    ========

    is_log_deriv_k_t_radical_in_field, is_deriv_k
    """
    if Df:
        dfa, dfd = (fd*derivation(fa, DE) - fa*derivation(fd, DE)).cancel(fd**2,
            include=True)
    else:
        dfa, dfd = fa, fd

    # Our assumption here is that each monomial is recursively transcendental
    if len(DE.exts) != len(DE.D) - 1:
        if [i for i in DE.cases if i == 'tan'] or \
                ({i for i, case in enumerate(DE.cases) if case == 'primitive'} -
                        set(DE.indices('log'))):
            raise NotImplementedError("Real version of the structure "
                "theorems with hypertangent support is not yet implemented.")

        # TODO: What should really be done in this case?
        raise NotImplementedError("Nonelementary extensions not supported "
            "in the structure theorems.")

    E_part = [DE.D[i].quo(Poly(DE.T[i], DE.T[i])).as_expr() for i in DE.indices('exp')]
    L_part = [DE.D[i].as_expr() for i in DE.indices('log')]

    # The expression dfa/dfd might not be polynomial in any of its symbols so we
    # use a Dummy as the generator for PolyMatrix.
    dum = Dummy()
    lhs = Matrix([E_part + L_part], dum)
    rhs = Matrix([dfa.as_expr()/dfd.as_expr()], dum)

    u = _structure_system_solve(lhs, rhs, DE)

    if u is None:
        # No constant solution
        return None
    else:
        if not all(i.is_Rational for i in u):
            # TODO: But maybe we can tell if they're not rational, like
            # log(2)/log(3). Also, there should be an option to continue
            # anyway, even if the result might potentially be wrong.
            raise NotImplementedError("Cannot work with non-rational "
                "coefficients in this case.")
        else:
            n = S.One*reduce(ilcm, [i.as_numer_denom()[1] for i in u])
            u = [n*i for i in u]
            terms = ([DE.T[i] for i in DE.indices('exp')] +
                    [DE.extargs[i - 1] for i in DE.indices('log')])
            ans = list(zip(terms, u))
            result = Mul(*[Pow(i, j) for i, j in ans])

            # exp(f) will be the same as result up to a multiplicative
            # constant.  We now find the log of that constant.
            argterms = ([DE.extargs[i - 1] for i in DE.indices('exp')] +
                    [DE.T[i] for i in DE.indices('log')])
            const = cancel(fa.as_expr()/fd.as_expr() -
                Add(*[Mul(i, j/n) for i, j in zip(argterms, u)]))

            return (ans, result, n, const)


def is_log_deriv_k_t_radical_in_field(fa, fd, DE, case='auto', z=None):
    """
    Checks if f can be written as the logarithmic derivative of a k(t)-radical.

    Explanation
    ===========

    It differs from is_log_deriv_k_t_radical(fa, fd, DE, Df=False)
    for any given fa, fd, DE in that it finds the solution in the
    given field not in some (possibly unspecified extension) and
    "in_field" with the function name is used to indicate that.

    f in k(t) can be written as the logarithmic derivative of a k(t) radical if
    there exist n in ZZ and u in k(t) with n, u != 0 such that n*f == Du/u.
    Either returns (n, u) or None, which means that f cannot be written as the
    logarithmic derivative of a k(t)-radical.

    case is one of {'primitive', 'exp', 'tan', 'auto'} for the primitive,
    hyperexponential, and hypertangent cases, respectively.  If case is 'auto',
    it will attempt to determine the type of the derivation automatically.

    This is the in-field logarithmic derivative of a k(t)-radical problem
    from Section 5.12 of Bronstein's book (the "Recognizing Logarithmic
    Derivatives of k(t)-radicals" subsection; neither edition gives it as
    named pseudocode).

    See also
    ========
    is_log_deriv_k_t_radical, is_deriv_k
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
        n = S.One*reduce(ilcm, [i.as_numer_denom()[1] for _, i in residueterms], 1)
        u = Mul(*[Pow(i, j*n) for i, j in residueterms])
        return (n, u)

    elif case == 'tan':
        raise NotImplementedError("The hypertangent case is "
        "not yet implemented for is_log_deriv_k_t_radical_in_field()")

    elif case in ('other_linear', 'other_nonlinear'):
        # XXX: If these are supported by the structure theorems, change to NotImplementedError.
        raise ValueError("The %s case is not supported in this function." % case)

    else:
        raise ValueError("case must be one of {'primitive', 'exp', 'tan', "
        "'base', 'auto'}, not %s" % case)

    common_denom = S.One*reduce(ilcm, [i.as_numer_denom()[1] for i in [j for _, j in
        residueterms]] + [n], 1)
    residueterms = [(i, j*common_denom) for i, j in residueterms]
    m = common_denom//n
    if common_denom != n*m:  # Verify exact division
        raise ValueError("Inexact division")
    u = cancel(u**m*Mul(*[Pow(i, j) for i, j in residueterms]))

    return (common_denom, u)

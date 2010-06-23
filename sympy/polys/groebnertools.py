"""Sparse distributed multivariate polynomials and Groebner bases. """

from sympy.polys.monomialtools import (
    monomial_mul,
    monomial_div,
    monomial_lcm,
    monomial_lex_key as O_lex,
    monomial_grlex_key as O_grlex,
    monomial_grevlex_key as O_grevlex,
)

from sympy.polys.polyerrors import (
    ExactQuotientFailed, DomainError,
)

from sympy.utilities import any, all
from operator import itemgetter

def sdp_LC(f, K):
    """Returns the leading coeffcient of `f`. """
    if not f:
        return K.zero
    else:
        return f[0][1]

def sdp_LM(f, u):
    """Returns the leading monomial of `f`. """
    if not f:
        return (0,)*(u+1)
    else:
        return f[0][0]

def sdp_LT(f, u, K):
    """Returns the leading term of `f`. """
    if f:
        return f[0]
    else:
        return (0,)*(u+1), K.zero

def sdp_del_LT(f):
    """Removes the leading from `f`. """
    return f[1:]

def sdp_coeffs(f):
    """Returns a list of monomials in `f`. """
    return [ coeff for _, coeff in f ]

def sdp_monoms(f):
    """Returns a list of monomials in `f`. """
    return [ monom for monom, _ in f ]

def sdp_sort(f, O):
    """Sort terms in `f` using the given monomial order `O`. """
    return sorted(f, key=lambda (m, _): O(m), reverse=True)

def sdp_strip(f):
    """Remove terms with zero coefficients from `f` in `K[X]`. """
    return [ (monom, coeff) for monom, coeff in f if coeff ]

def sdp_normal(f, K):
    """Normalize distributed polynomial in the given domain. """
    return [ (monom, K.convert(coeff)) for monom, coeff in f if coeff ]

def sdp_from_dict(f, O):
    """Make a distributed polynomial from a dictionary. """
    return sdp_sort(f.items(), O)

def sdp_to_dict(f):
    """Make a dictionary from a distributed polynomial. """
    return dict(f)

def sdp_indep_p(f, j, u):
    """Returns `True` if a polynomial is independent of `x_j`. """
    if j < 0 or j > u:
        raise IndexError("-%s <= j < %s expected, got %s" % (u, u, j))
    else:
        return all(not monom[j] for monom in sdp_monoms(h))

def sdp_one_p(f, u, K):
    """Returns True if `f` is a multivariate one in `K[X]`. """
    return f == sdp_one(u, K)

def sdp_one(u, K):
    """Returns a multivariate one in `K[X]`. """
    return (((0,)*(u+1), K.one),)

def sdp_term_p(f):
    """Returns True if `f` has a single term or is zero. """
    return len(f) <= 1

def sdp_abs(f, u, O, K):
    """Make all coefficients positive in `K[X]`. """
    return [ (monom, K.abs(coeff)) for monom, coeff in f ]

def sdp_neg(f, u, O, K):
    """Negate a polynomial in `K[X]`. """
    return [ (monom, -coeff) for monom, coeff in f ]

def sdp_add_term(f, (M, c), u, O, K):
    """Add a single term using bisection method. """
    if not c:
        return f
    if not f:
        return [(M, c)]

    monoms = sdp_monoms(f)

    if cmp(O(M), O(monoms[ 0])) > 0:
        return [(M, c)] + f
    if cmp(O(M), O(monoms[-1])) < 0:
        return f + [(M, c)]

    lo, hi = 0, len(monoms)-1

    while lo <= hi:
        i = (lo + hi) // 2
        j = cmp(O(M), O(monoms[i]))

        if not j:
            coeff = f[i][1] + c

            if not coeff:
                return f[:i] + f[i+1:]
            else:
                return f[:i] + [(M, coeff)] + f[i+1:]
        else:
            if j > 0:
                hi = i - 1
            else:
                lo = i + 1
    else:
        return f[:i] + [(M, c)] + f[i+1:]

def sdp_sub_term(f, (M, c), u, O, K):
    """Sub a single term using bisection method. """
    if not c:
        return f
    if not f:
        return [(M, -c)]

    monoms = sdp_monoms(f)

    if cmp(O(M), O(monoms[ 0])) > 0:
        return [(M, -c)] + f
    if cmp(O(M), O(monoms[-1])) < 0:
        return f + [(M, -c)]

    lo, hi = 0, len(monoms)-1

    while lo <= hi:
        i = (lo + hi) // 2
        j = cmp(O(M), O(monoms[i]))

        if not j:
            coeff = f[i][1] - c

            if not coeff:
                return f[:i] + f[i+1:]
            else:
                return f[:i] + [(M, coeff)] + f[i+1:]
        else:
            if j > 0:
                hi = i - 1
            else:
                lo = i + 1
    else:
        return f[:i] + [(M, -c)] + f[i+1:]

def sdp_mul_term(f, (M, c), u, O, K):
    """Multiply a distributed polynomial by a term. """
    if not f or not c:
        return []
    else:
        if K.is_one(c):
            return [ (monomial_mul(f_M, M), f_c) for f_M, f_c in f ]
        else:
            return [ (monomial_mul(f_M, M), f_c*c) for f_M, f_c in f ]

def sdp_add(f, g, u, O, K):
    """Add distributed polynomials in `K[X]`. """
    h = dict(f)

    for monom, c in g:
        if h.has_key(monom):
            coeff = h[monom] + c

            if not coeff:
                del h[monom]
            else:
                h[monom] = coeff
        else:
            h[monom] = c

    return sdp_from_dict(h, O)

def sdp_sub(f, g, u, O, K):
    """Subtract distributed polynomials in `K[X]`. """
    h = dict(f)

    for monom, c in g:
        if h.has_key(monom):
            coeff = h[monom] - c

            if not coeff:
                del h[monom]
            else:
                h[monom] = coeff
        else:
            h[monom] = -c

    return sdp_from_dict(h, O)

def sdp_mul(f, g, u, O, K):
    """Multiply distributed polynomials in `K[X]`. """
    if sdp_term_p(f):
        if not f:
            return f
        else:
            return sdp_mul_term(g, f[0], u, O, K)

    if sdp_term_p(g):
        if not g:
            return g
        else:
            return sdp_mul_term(f, g[0], u, O, K)

    h = {}

    for fm, fc in f:
        for gm, gc in g:
            monom = monomial_mul(fm, gm)
            coeff = fc*gc

            if h.has_key(monom):
                coeff += h[monom]

                if not coeff:
                    del h[monom]
                    continue

            h[monom] = coeff

    return sdp_from_dict(h, O)

def sdp_sqr(f, u, O, K):
    """Square a distributed polynomial in `K[X]`. """
    h = {}

    for fm, fc in f:
        for Fm, Fc in f:
            monom = monomial_mul(fm, Fm)
            coeff = fc*Fc

            if h.has_key(monom):
                coeff += h[monom]

                if not coeff:
                    del h[monom]
                    continue

            h[monom] = coeff

    return sdp_from_dict(h, O)

def sdp_pow(f, n, u, O, K):
    """Raise `f` to the n-th power in `K[X]`. """
    if not n:
        return sdp_one(u, K)
    if n < 0:
        raise ValueError("can't raise a polynomial to negative power")
    if n == 1 or not f or sdp_one_p(f, u, K):
        return f

    g = sdp_one(u, K)

    while True:
        n, m = n//2, n

        if m & 1:
            g = sdp_mul(g, f, u, O, K)

            if not n:
                break

        f = sdp_sqr(f, u, O, K)

    return g

def sdp_monic(f, K):
    """Divides all coefficients by `LC(f)` in `K[X]`. """
    if not f:
        return f

    lc_f = sdp_LC(f, K)

    if K.is_one(lc_f):
        return f
    else:
        return [ (m, K.quo(c, lc_f)) for m, c in f ]

def sdp_content(f, K):
    """Returns GCD of coefficients in `K[X]`. """
    if K.has_Field:
        return K.one
    else:
        cont = K.zero

        for _, c in f:
            cont = K.gcd(cont, c)

            if K.is_one(cont):
                break

        return cont

def sdp_primitive(f, K):
    """Returns content and a primitive polynomial in `K[X]`. """
    if K.has_Field:
        return K.one, f
    else:
        cont = sdp_content(f, K)

        if K.is_one(cont):
            return cont, f
        else:
            return cont, [ (m, K.exquo(c, cont)) for m, c in f ]

def _term_rr_div(a, b, K):
    """Division of two terms in over a ring. """
    a_lm, a_lc = a
    b_lm, b_lc = b

    monom = monomial_div(a_lm, b_lm)

    if not (monom is None or a_lc % b_lc):
        return monom, K.exquo(a_lc, b_lc)
    else:
        return None

def _term_ff_div(a, b, K):
    """Division of two terms in over a field. """
    a_lm, a_lc = a
    b_lm, b_lc = b

    monom = monomial_div(a_lm, b_lm)

    if monom is not None:
        return monom, K.exquo(a_lc, b_lc)
    else:
        return None

def sdp_div(f, G, u, O, K):
    """Generalized polynomial division with remainder in `K[X]`.

       Given polynomial `f` and a set of polynomials `g = (g_1, ..., g_n)`
       compute a set of quotients `q = (q_1, ..., q_n)` and remainder `r`
       such that `f = q_1*f_1 + ... + q_n*f_n + r`, where `r = 0` or `r`
       is a completely reduced polynomial with respect to `g`.

       References
       ==========

       .. [Cox97] D. Cox, J. Little, D. O'Shea, Ideals, Varieties and
           Algorithms, Springer, Second Edition, 1997, pp. 62

       .. [Ajwa95] I.A. Ajwa, Z. Liu, P.S. Wang, Groebner Bases Algorithm,
           http://citeseer.ist.psu.edu/ajwa95grbner.html, 1995

    """
    Q, r = [ [] for _ in xrange(len(G)) ], []

    if K.has_Field:
        term_div = _term_ff_div
    else:
        term_div = _term_rr_div

    while f:
        for i, g in enumerate(G):
            tq = term_div(sdp_LT(f, u, K), sdp_LT(g, u, K), K)

            if tq is not None:
                Q[i] = sdp_add_term(Q[i], tq, u, O, K)
                f = sdp_sub(f, sdp_mul_term(g, tq, u, O, K), u, O, K)

                break
        else:
            r = sdp_add_term(r, sdp_LT(f, u, K), u, O, K)
            f = sdp_del_LT(f)

    return Q, r

def sdp_rem(f, g, u, O, K):
    """Returns polynomial remainder in `K[X]`. """
    return sdp_div(f, g, u, O, K)[1]

def sdp_quo(f, g, u, O, K):
    """Returns polynomial quotient in `K[X]`. """
    q, r = sdp_div(f, g, u, O, K)

    if not r:
        return q
    else:
        raise ExactQuotientFailed('%s does not divide %s in %s' % (g, f, K))

def sdp_exquo(f, g, u, O, K):
    """Returns exact polynomial quotient in `K[x]`. """
    return sdp_div(f, g, u, O, K)[0]

def sdp_lcm(f, g, u, O, K):
    """Computes LCM of two polynomials in `K[X]`.

       The LCM is computed as the unique generater of the intersection
       of the two ideals generated by `f` and `g`. The approach is to
       compute a Groebner basis with respect to lexicographic ordering
       of `t*f` and `(1 - t)*g`, where `t` is an unrealted variable and
       then filtering out the solution that doesn't contain `t`.

       References
       ==========

       .. [Cox97] D. Cox, J. Little, D. O'Shea, Ideals, Varieties and
           Algorithms, Springer, Second Edition, 1997, pp. 187

    """
    if not f or not g:
        return []

    if sdp_term_p(f) and sdp_term_p(g):
        monom = monomial_lcm(sdp_LM(f, u), sdp_LM(g, u))

        fc, gc = sdp_LC(f, K), sdp_LC(g, K)

        if K.has_Field:
            coeff = K.one
        else:
            coeff = K.lcm(fc, gc)

        return [(monom, coeff)]

    if not K.has_Field:
        lcm = K.one
    else:
        fc, f = sdp_primitive(f, K)
        gc, g = sdp_primitive(g, K)

        lcm = K.lcm(fc, gc)

    f_terms = tuple( ((1,) + m,  c) for m, c in f )
    g_terms = tuple( ((0,) + m,  c) for m, c in g ) \
            + tuple( ((1,) + m, -c) for m, c in g )

    F = sdp_sort(f_terms, O_lex)
    G = sdp_sort(g_terms, O_lex)

    basis = sdp_groebner([F, G], u, O_lex, K)

    H = [ h for h in basis if sdp_indep_p(h, 0, u) ]

    if K.is_one(lcm):
        h = [ (m[1:], c)     for m, c in H[0] ]
    else:
        h = [ (m[1:], c*lcm) for m, c in H[0] ]

    return sdp_sort(h, O)

def sdp_gcd(f, g, u, O, K):
    """Compute GCD of two polynomials in `K[X]` via LCM. """
    if not K.has_Field:
        fc, f = sdp_primitive(f, K)
        gc, g = sdp_primitive(g, K)

        gcd = K.gcd(fc, gc)

    h = sdp_quo(sdp_mul(f, g, u, O, K),
               sdp_lcm(f, g, u, O, K), u, O, K)

    if not K.has_Field:
        if K.is_one(gcd):
            return h
        else:
            return [ (m, c*gcd) for m, c in h ]
    else:
        return sdp_monic(h, K)

def sdp_groebner(F, u, O, K, monic=True):
    """Computes Groebner basis for a set of polynomials in `K[X]`.

       Given a set of multivariate polynomials `F`, finds another
       set `G`, such that Ideal `F = Ideal G` and `G` is a reduced
       Groebner basis.

       The resulting basis is unique and has monic generators if the
       ground domains is a field. Otherwise the result is non-unique
       but Groebner bases over e.g. integers can be computed (if the
       input polynomials are monic).

       Groebner bases can be used to choose specific generators for a
       polynomial ideal. Because these bases are unique you can check
       for ideal equality by comparing the Groebner bases.  To see if
       one polynomial lies in an ideal, divide by the elements in the
       base and see if the remainder vanishes.

       They can also be used to  solve systems of polynomial equations
       as,  by choosing lexicographic ordering,  you can eliminate one
       variable at a time, provided that the ideal is zero-dimensional
       (finite number of solutions).

       References
       ==========

       .. [Bose03] N.K. Bose, B. Buchberger, J.P. Guiver, Multidimensional
           Systems Theory and Applications, Springer, 2003, pp. 98+

       .. [Giovini91] A. Giovini, T. Mora, "One sugar cube, please" or
           Selection strategies in Buchberger algorithm, ISSAC '91, ACM

       .. [Ajwa95] I.A. Ajwa, Z. Liu, P.S. Wang, Groebner Bases Algorithm,
           http://citeseer.ist.psu.edu/ajwa95grbner.html, 1995

       .. [Cox97] D. Cox, J. Little, D. O'Shea, Ideals, Varieties and
           Algorithms, Springer, Second Edition, 1997, pp. 62

    """
    if not K.has_Field:
        raise DomainError("can't compute a Groebner basis over %s" % K)

    F = [ f for f in F if f ]

    if not F:
        return [[]]

    R, P, G, B, I = set(), set(), set(), {}, {}

    for i, f in enumerate(F):
        I[tuple(f)] = i
        R.add(i)

    def normal(g, J):
        h = sdp_rem(g, [ F[j] for j in J ], u, O, K)

        if not h:
            return None
        else:
            H = tuple(h)

            if not H in I:
                I[H] = len(F)
                F.append(h)

            return I[H], sdp_LM(h, u)

    def generate(R, P, G, B):
        while R:
            h = normal(F[R.pop()], G | P)

            if h is not None:
                k, LM = h

                G0 = set(g for g in G if monomial_div(sdp_LM(F[g], u), LM))
                P0 = set(p for p in P if monomial_div(sdp_LM(F[p], u), LM))

                G, P, R = G - G0, P - P0 | set([k]), R | G0 | P0

                for i, j in set(B):
                    if i in G0 or j in G0:
                        del B[(i, j)]

        G |= P

        for i in G:
            for j in P:
                if i == j:
                    continue

                if i < j:
                   k = (i, j)
                else:
                   k = (j, i)

                if k not in B:
                    B[k] = monomial_lcm(sdp_LM(F[i], u), sdp_LM(F[j], u))

        G = set([ normal(F[g], G - set([g]))[0] for g in G ])

        return R, P, G, B

    R, P, G, B = generate(R, P, G, B)

    while B:
        k, M = B.items()[0]

        for l, N in B.iteritems():
            if cmp(O(M), O(N)) == 1:
                k, M = l, N

        del B[k]

        i, j = k[0], k[1]
        p, q = F[i], F[j]

        p_LM, q_LM = sdp_LM(p, u), sdp_LM(q, u)

        if M == monomial_mul(p_LM, q_LM):
            continue

        criterion = False

        for g in G:
            if g == i or g == j:
                continue

            if (min(i, g), max(i, g)) not in B:
                continue

            if (min(j, g), max(j, g)) not in B:
                continue

            if not monomial_div(M, sdp_LM(F[g], u)):
                continue

            criterion = True
            break

        if criterion:
            continue

        p = sdp_mul_term(p, (monomial_div(M, p_LM), K.quo(K.one, sdp_LC(p, K))), u, O, K)
        q = sdp_mul_term(q, (monomial_div(M, q_LM), K.quo(K.one, sdp_LC(q, K))), u, O, K)

        h = normal(sdp_sub(p, q, u, O, K), G)

        if h is not None:
            k, LM = h

            G0 = set(g for g in G if monomial_div(sdp_LM(F[g], u), LM))

            R, P, G = G0, set([k]), G - G0

            for i, j in set(B):
                if i in G0 or j in G0:
                    del B[(i, j)]

            R, P, G, B = generate(R, P, G, B)

    if not monic:
        basis = [ F[g] for g in G ]
    else:
        basis = [ sdp_monic(F[g], K) for g in G ]

    return sorted(basis, key=lambda f: O(sdp_LM(f, u)), reverse=True)


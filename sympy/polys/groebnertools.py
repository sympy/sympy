"""Sparse distributed multivariate polynomials and Groebner bases. """

from sympy.core.compatibility import cmp

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
    return sorted(f, key=lambda term: O(term[0]), reverse=True)

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

def sdp_add_term(f, term, u, O, K):
    """Add a single term using bisection method. """
    M, c = term

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

def sdp_sub_term(f, term, u, O, K):
    """Sub a single term using bisection method. """
    M, c = term

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

def sdp_mul_term(f, term, u, O, K):
    """Multiply a distributed polynomial by a term. """
    M, c = term

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
        if monom in h:
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
        if monom in h:
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

            if monom in h:
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

            if monom in h:
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
            return cont, [ (m, K.quo(c, cont)) for m, c in f ]

def _term_rr_div(a, b, K):
    """Division of two terms in over a ring. """
    a_lm, a_lc = a
    b_lm, b_lc = b

    monom = monomial_div(a_lm, b_lm)

    if not (monom is None or a_lc % b_lc):
        return monom, K.quo(a_lc, b_lc)
    else:
        return None

def _term_ff_div(a, b, K):
    """Division of two terms in over a field. """
    a_lm, a_lc = a
    b_lm, b_lc = b

    monom = monomial_div(a_lm, b_lm)

    if monom is not None:
        return monom, K.quo(a_lc, b_lc)
    else:
        return None

def sdp_div(f, G, u, O, K):
    """
    Generalized polynomial division with remainder in `K[X]`.

    Given polynomial `f` and a set of polynomials `g = (g_1, ..., g_n)`
    compute a set of quotients `q = (q_1, ..., q_n)` and remainder `r`
    such that `f = q_1*f_1 + ... + q_n*f_n + r`, where `r = 0` or `r`
    is a completely reduced polynomial with respect to `g`.

    **References**

    1. [Cox97]_
    2. [Ajwa95]_

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
    """Returns polynomial quotient in `K[x]`. """
    return sdp_div(f, g, u, O, K)[0]

def sdp_exquo(f, g, u, O, K):
    """Returns exact polynomial quotient in `K[X]`. """
    q, r = sdp_div(f, g, u, O, K)

    if not r:
        return q
    else:
        raise ExactQuotientFailed(f, g)

def sdp_lcm(f, g, u, O, K):
    """
    Computes LCM of two polynomials in `K[X]`.

    The LCM is computed as the unique generater of the intersection
    of the two ideals generated by `f` and `g`. The approach is to
    compute a Groebner basis with respect to lexicographic ordering
    of `t*f` and `(1 - t)*g`, where `t` is an unrelated variable and
    then filtering out the solution that doesn't contain `t`.

    **References**

    1. [Cox97]_

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

def sdp_groebner(f, u, O, K, gens='', verbose=False):
    """
    Computes Groebner basis for a set of polynomials in `K[X]`.

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

    **References**

    1. [Bose03]_
    2. [Giovini91]_
    3. [Ajwa95]_
    4. [Cox97]_

    Algorithm used: an improved version of Buchberger's algorithm
    as presented in T. Becker, V. Weispfenning, Groebner Bases: A
    Computational Approach to Commutative Algebra, Springer, 1993,
    page 232.

    Added optional ``gens`` argument to apply :func:`sdp_str` for
    the purpose of debugging the algorithm.

    """
    if not K.has_Field:
        raise DomainError("can't compute a Groebner basis over %s" % K)

    def select(P):
        # normal selection strategy
        # select the pair with minimum LCM(LM(f), LM(g))
        pr = min(P, key=lambda pair: O(monomial_lcm(sdp_LM(f[pair[0]], u), sdp_LM(f[pair[1]], u))))
        return pr

    def normal(g, J):
        h = sdp_rem(g, [ f[j] for j in J ], u, O, K)

        if not h:
            return None
        else:
            h = sdp_monic(h, K)
            h = tuple(h)

            if not h in I:
                I[h] = len(f)
                f.append(h)

            return sdp_LM(h, u), I[h]

    def update(G, B, ih):
        # update G using the set of critical pairs B and h
        # [BW] page 230
        h = f[ih]
        mh = sdp_LM(h, u)

        # filter new pairs (h, g), g in G
        C = G.copy()
        D = set()

        while C:
            # select a pair (h, g) by popping an element from C
            ig = C.pop()
            g = f[ig]
            mg = sdp_LM(g, u)
            LCMhg = monomial_lcm(mh, mg)

            def lcm_divides(ip):
                # LCM(LM(h), LM(p)) divides LCM(LM(h), LM(g))
                m = monomial_lcm(mh, sdp_LM(f[ip], u))
                return monomial_div(LCMhg, m)

            # HT(h) and HT(g) disjoint: mh*mg == LCMhg
            if monomial_mul(mh, mg) == LCMhg or (
                not any(lcm_divides(ipx) for ipx in C) and
                not any(lcm_divides(pr[1]) for pr in D)):
                  D.add((ih, ig))

        E = set()

        while D:
            # select h, g from D (h the same as above)
            ih, ig = D.pop()
            mg = sdp_LM(f[ig], u)
            LCMhg = monomial_lcm(mh, mg)

            if not monomial_mul(mh, mg) == LCMhg:
                E.add((ih, ig))

        # filter old pairs
        B_new = set()

        while B:
            # select g1, g2 from B (-> CP)
            ig1, ig2 = B.pop()
            mg1 = sdp_LM(f[ig1], u)
            mg2 = sdp_LM(f[ig2], u)
            LCM12 = monomial_lcm(mg1, mg2)

            # if HT(h) does not divide lcm(HT(g1), HT(g2))
            if not monomial_div(LCM12, mh) or \
                monomial_lcm(mg1, mh) == LCM12 or \
                monomial_lcm(mg2, mh) == LCM12:
              B_new.add((ig1, ig2))

        B_new |= E

        # filter polynomials
        G_new = set()

        while G:
            ig = G.pop()
            mg = sdp_LM(f[ig], u)

            if not monomial_div(mg, mh):
                G_new.add(ig)

        G_new.add(ih)

        return G_new, B_new
      # end of update ################################

    if not f:
        return []

    # replace f with a reduced list of initial polynomials; see [BW] page 203
    f1 = f[:]

    while True:
        f = f1[:]
        f1 = []

        for i in range(len(f)):
            p = f[i]
            r = sdp_rem(p, f[:i], u, O, K)

            if r:
               f1.append(sdp_monic(r, K))

        if f == f1:
            break

    f = [tuple(p) for p in f]
    I = {}            # ip = I[p]; p = f[ip]
    F = set()         # set of indices of polynomials
    G = set()         # set of indices of intermediate would-be Groebner basis
    CP = set()        # set of pairs of indices of critical pairs

    for i, h in enumerate(f):
        I[h] = i
        F.add(i)

    #####################################
    # algorithm GROEBNERNEWS2 in [BW] page 232
    while F:
        # select p with minimum monomial according to the monomial ordering O
        h = min([f[x] for x in F], key=lambda f: O(sdp_LM(f, u)))
        ih = I[h]
        F.remove(ih)
        G, CP = update(G, CP, ih)

    # count the number of critical pairs which reduce to zero
    reductions_to_zero = 0

    while CP:
        ig1, ig2 = select(CP)
        CP.remove((ig1, ig2))

        h = sdp_spoly(f[ig1], f[ig2], u, O, K)
        # ordering divisors is on average more efficient [Cox] page 111
        G1 = sorted(G, key=lambda g: O(sdp_LM(f[g], u)))
        ht = normal(h, G1)

        if ht:
            G, CP = update(G, CP, ht[1])
        else:
            reductions_to_zero += 1

    ######################################
    # now G is a Groebner basis; reduce it
    Gr = set()

    for ig in G:
        ht = normal(f[ig], G - set([ig]))

        if ht:
            Gr.add(ht[1])

    Gr = [list(f[ig]) for ig in Gr]

    # order according to the monomial ordering
    Gr = sorted(Gr, key=lambda f: O(sdp_LM(f, u)), reverse=True)

    if verbose:
        print 'reductions_to_zero = %d' % reductions_to_zero

    return Gr

def sdp_str(f, gens):
    if isinstance(gens, basestring):
        gens = gens.split(',')
    ngens = len(gens)
    z = (0,)*ngens
    s = ''
    for expv, c in f:
        if c > 0:
            s += ' +'
        else:
            s += ' -'
        if c < 0:
            c = -c
        if c != 1: # and expv != z:
            cnt1 = str(c)
        else:
            cnt1 = ''
        sa = []
        for i in range(ngens):
            exp = expv[i]
            if exp > 1:
                sa.append('%s^%d' % (gens[i], exp))
            if exp == 1:
                sa.append('%s' % gens[i])
        if cnt1:
            sa = [cnt1] + sa
        s += '*'.join(sa)
    return s

def sdp_spoly(p1, p2, u, O, K):
    """
    Compute LCM(LM(p1), LM(p2))/LM(p1)*p1 - LCM(LM(p1), LM(p2))/LM(p2)*p2
    This is the S-poly provided p1 and p2 are monic
    """
    LM1 = sdp_LM(p1, u)
    LM2 = sdp_LM(p2, u)
    LCM12 = monomial_lcm(LM1, LM2)
    m1 = monomial_div(LCM12, LM1)
    m2 = monomial_div(LCM12, LM2)
    s1 = sdp_mul_term(p1, (m1, K.one), u, O, K)
    s2 = sdp_mul_term(p2, (m2, K.one), u, O, K)
    s = sdp_sub(s1, s2, u, O, K)
    return s


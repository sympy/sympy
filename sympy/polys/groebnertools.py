"""Groebner bases algorithms. """

from sympy.polys.monomialtools import (
    monomial_mul, monomial_div, monomial_lcm,
)

from sympy.polys.distributedpolys import (
    sdp_LC, sdp_LM, sdp_LT, sdp_mul_term,
    sdp_sub, sdp_mul_term, sdp_monic,
    sdp_rem, sdp_strip, sdp_sort,
    _term_ff_div, _term_rr_div,
)

from sympy.polys.polyerrors import (
    DomainError,
)

from sympy.polys.polyconfig import query

def sdp_groebner(f, u, O, K, gens='', verbose=False, method=None):
    """
    Computes Groebner basis for a set of polynomials in `K[X]`.

    Wrapper around the (default) improved Buchberger and the other algorithms
    for computing Groebner bases. The choice of algorithm can be changed via
    ``method`` argument or :func:`setup` from :mod:`sympy.polys.polyconfig`,
    where ``method`` can be either ``buchberger`` or ``f5b``.

    """
    if method is None:
        method = query('GB_METHOD')

    _groebner_methods = {
        'buchberger': buchberger,
        'f5b': f5b,
    }

    try:
        func = _groebner_methods[method]
    except KeyError:
        raise ValueError("'%s' is not a valid Groebner bases algorithm (valid are 'buchberger' and 'f5b')" % method)
    else:
        return func(f, u, O, K, gens, verbose)

# Buchberger algorithm

def buchberger(f, u, O, K, gens='', verbose=False):
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

    They can also be used to solve systems of polynomial equations
    as,  by choosing lexicographic ordering,  you can eliminate one
    variable at a time, provided that the ideal is zero-dimensional
    (finite number of solutions).

    References
    ==========

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
    s = ''
    for expv, c in f:
        if c > 0:
            s += ' +'
        else:
            s += ' -'
        if c < 0:
            c = -c
        if c != 1:  # and expv != z:
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

# F5B

# convenience functions


def Sign(f):
    return f[0]


def Polyn(f):
    return f[1]


def Num(f):
    return f[2]


def sig(monomial, index):
    return (monomial, index)


def lbp(signature, polynomial, number):
    return (signature, polynomial, number)

# signature functions


def sig_cmp(u, v, O):
    """
    Compare two signatures by extending the term order to K[X]^n.

    u < v iff
        - the index of v is greater than the index of u
    or
        - the index of v is equal to the index of u and u[0] < v[0] w.r.t. O

    u > v otherwise
    """
    if u[1] > v[1]:
        return -1
    if u[1] == v[1]:
        #if u[0] == v[0]:
        #    return 0
        if O(u[0]) < O(v[0]):
            return -1
    return 1


def sig_key(s, O):
    """
    Key for comparing two signatures.

    s = (m, k), t = (n, l)

    s < t iff [k > l] or [k == l and m < n]
    s > t otherwise
    """
    return (-s[1], O(s[0]))


def sig_mult(s, m):
    """
    Multiply a signature by a monomial.

    The product of a signature (m, i) and a monomial n is defined as
    (m * t, i).
    """
    return sig(monomial_mul(s[0], m), s[1])

# labeled polynomial functions


def lbp_sub(f, g, u, O, K):
    """
    Subtract labeled polynomial g from f.

    The signature and number of the difference of f and g are signature
    and number of the maximum of f and g, w.r.t. lbp_cmp.
    """
    if sig_cmp(Sign(f), Sign(g), O) < 0:
        max_poly = g
    else:
        max_poly = f

    ret = sdp_sub(Polyn(f), Polyn(g), u, O, K)

    return lbp(Sign(max_poly), ret, Num(max_poly))


def lbp_mul_term(f, cx, u, O, K):
    """
    Multiply a labeled polynomial with a term.

    The product of a labeled polynomial (s, p, k) by a monomial is
    defined as (m * s, m * p, k).
    """
    return lbp(sig_mult(Sign(f), cx[0]), sdp_mul_term(Polyn(f), cx, u, O, K), Num(f))


def lbp_cmp(f, g, O):
    """
    Compare two labeled polynomials.

    f < g iff
        - Sign(f) < Sign(g)
    or
        - Sign(f) == Sign(g) and Num(f) > Num(g)

    f > g otherwise
    """
    if sig_cmp(Sign(f), Sign(g), O) == -1:
        return -1
    if Sign(f) == Sign(g):
        if Num(f) > Num(g):
            return -1
        #if Num(f) == Num(g):
        #    return 0
    return 1


def lbp_key(f, O):
    """
    Key for comparing two labeled polynomials.
    """
    return (sig_key(Sign(f), O), -Num(f))

# algorithm and helper functions


def critical_pair(f, g, u, O, K):
    """
    Compute the critical pair corresponding to two labeled polynomials.

    A critical pair is a tuple (um, f, vm, g), where um and vm are
    terms such that um * f - vm * g is the S-polynomial of f and g (so,
    wlog assume um * f > vm * g).
    For performance sake, a critical pair is represented as a tuple
    (Sign(um * f), um, f, Sign(vm * g), vm, g), since um * f creates
    a new, relatively expensive object in memory, whereas Sign(um *
    f) and um are lightweight and f (in the tuple) is a reference to
    an already existing object in memory.
    """
    ltf = sdp_LT(Polyn(f), u, K)
    ltg = sdp_LT(Polyn(g), u, K)
    lt = (monomial_lcm(ltf[0], ltg[0]), K.one)

    if K.has_Field:
        term_div = _term_ff_div
    else:
        term_div = _term_rr_div

    um = term_div(lt, ltf, K)
    vm = term_div(lt, ltg, K)

    # The full information is not needed (now), so only the product
    # with the leading term is considered:
    fr = lbp_mul_term(lbp(Sign(f), [sdp_LT(Polyn(f), u, K)], Num(f)), um, u, O, K)
    gr = lbp_mul_term(lbp(Sign(g), [sdp_LT(Polyn(g), u, K)], Num(g)), vm, u, O, K)

    # return in proper order, such that the S-polynomial is just
    # u_first * f_first - u_second * f_second:
    if lbp_cmp(fr, gr, O) == -1:
        return (Sign(gr), vm, g, Sign(fr), um, f)
    else:
        return (Sign(fr), um, f, Sign(gr), vm, g)


def cp_cmp(c, d, O):
    """
    Compare two critical pairs c and d.

    c < d iff
        - lbp(c[0], _, Num(c[2]) < lbp(d[0], _, Num(d[2])) (this
        corresponds to um_c * f_c and um_d * f_d)
    or
        - lbp(c[0], _, Num(c[2]) >< lbp(d[0], _, Num(d[2])) and
        lbp(c[3], _, Num(c[5])) < lbp(d[3], _, Num(d[5])) (this
        corresponds to vm_c * g_c and vm_d * g_d)

    c > d otherwise
    """
    c0 = lbp(c[0], [], Num(c[2]))
    d0 = lbp(d[0], [], Num(d[2]))

    r = lbp_cmp(c0, d0, O)

    if r == -1:
        return -1
    if r == 0:
        c1 = lbp(c[3], [], Num(c[5]))
        d1 = lbp(d[3], [], Num(d[5]))

        r = lbp_cmp(c1, d1, O)

        if r == -1:
            return -1
        #if r == 0:
        #    return 0
    return 1


def cp_key(c, O):
    """
    Key for comparing critical pairs.
    """
    return (lbp_key(lbp(c[0], [], Num(c[2])), O), lbp_key(lbp(c[3], [], Num(c[5])), O))


def s_poly(cp, u, O, K):
    """
    Compute the S-polynomial of a critical pair.

    The S-polynomial of a critical pair cp is cp[1] * cp[2] - cp[4] * cp[5].
    """
    return lbp_sub(lbp_mul_term(cp[2], cp[1], u, O, K), lbp_mul_term(cp[5], cp[4], u, O, K), u, O, K)


def is_rewritable_or_comparable(sign, num, B, u, K):
    """
    Check if a labeled polynomial is redundant by checking if its
    signature and number imply rewritability or comparability.

    (sign, num) is comparable if there exists a labeled polynomial
    h in B, such that sign[1] (the index) is less than Sign(h)[1]
    and sign[0] is divisible by the leading monomial of h.

    (sign, num) is rewritable if there exists a labeled polynomial
    h in B, such thatsign[1] is equal to Sign(h)[1], num < Num(h)
    and sign[0] is divisible by Sign(h)[0].
    """
    for h in B:
        # comparable
        if sign[1] < Sign(h)[1]:
            if monomial_divides(sign[0], sdp_LM(Polyn(h), u)):
                return True

        # rewritable
        if sign[1] == Sign(h)[1]:
            if num < Num(h):
                if monomial_divides(sign[0], Sign(h)[0]):
                    return True
    return False


def f5_reduce(f, B, u, O, K):
    """
    F5-reduce a labeled polynomial f by B.

    Continously searches for non-zero labeled polynomial h in B, such
    that the leading term lt_h of h divides the leading term lt_f of
    f and Sign(lt_h * h) < Sign(f). If such a labeled polynomial h is
    found, f gets replaced by f - lt_f / lt_h * h. If no such h can be
    found or f is 0, f is no further F5-reducible and f gets returned.

    A polynomial that is reducible in the usual sense (sdp_rem)
    need not be F5-reducible, e.g.:

    >>> from sympy.polys.groebnertools import lbp, sig, f5_reduce, Polyn
    >>> from sympy.polys.distributedpolys import sdp_rem
    >>> from sympy.polys.monomialtools import lex
    >>> from sympy import QQ

    >>> f = lbp(sig((1, 1, 1), 4), [((1, 0, 0), QQ(1))], 3)
    >>> g = lbp(sig((0, 0, 0), 2), [((1, 0, 0), QQ(1))], 2)

    >>> sdp_rem(Polyn(f), [Polyn(g)], 2, lex, QQ)
    []
    >>> f5_reduce(f, [g], 2, lex, QQ)
    (((1, 1, 1), 4), [((1, 0, 0), 1/1)], 3)

    """
    if Polyn(f) == []:
        return f

    if K.has_Field:
        term_div = _term_ff_div
    else:
        term_div = _term_rr_div

    while True:
        g = f

        for h in B:
            if Polyn(h) != []:
                if monomial_divides(sdp_LM(Polyn(f), u), sdp_LM(Polyn(h), u)):
                    t = term_div(sdp_LT(Polyn(f), u, K), sdp_LT(Polyn(h), u, K), K)
                    if sig_cmp(sig_mult(Sign(h), t[0]), Sign(f), O) < 0:
                        # The following check need not be done and is in general slower than without.
                        #if not is_rewritable_or_comparable(Sign(gp), Num(gp), B, u, K):
                        hp = lbp_mul_term(h, t, u, O, K)
                        f = lbp_sub(f, hp, u, O, K)
                        break

        if g == f or Polyn(f) == []:
            return f


def f5b(F, u, O, K, gens='', verbose=False):
    """
    Computes a reduced Groebner basis for the ideal generated by F.

    f5b is an implementation of the F5B algorithm by Yao Sun and
    Dingkang Wang. Similarly to Buchberger's algorithm, the algorithm
    proceeds by computing critical pairs, computing the S-polynomial,
    reducing it and adjoining the reduced S-polynomial if it is not 0.

    Unlike Buchberger's algorithm, each polynomial contains additional
    information, namely a signature and a number. The signature
    specifies the path of computation (i.e. from which polynomial in
    the original basis was it derived and how), the number says when
    the polynomial was added to the basis.  With this information it
    is (often) possible to decide if an S-polynomial will reduce to
    0 and can be discarded.

    Optimizations include: Reducing the generators before computing
    a Groebner basis, removing redundant critical pairs when a new
    polynomial enters the basis and sorting the critical pairs and
    the current basis.

    Once a Groebner basis has been found, it gets reduced.

    ** References **
    Yao Sun, Dingkang Wang: "A New Proof for the Correctness of F5
    (F5-Like) Algorithm", http://arxiv.org/abs/1004.0084 (specifically
    v4)

    Thomas Becker, Volker Weispfenning, Groebner bases: A computational
    approach to commutative algebra, 1993, p. 203, 216
    """
    if not K.has_Field:
        raise DomainError("can't compute a Groebner basis over %s" % K)

    # reduce polynomials (like in Mario Pernici's implementation) (Becker, Weispfenning, p. 203)
    B = F
    while True:
        F = B
        B = []

        for i in xrange(len(F)):
            p = F[i]
            r = sdp_rem(p, F[:i], u, O, K)

            if r != []:
                B.append(r)

        if F == B:
            break

    # basis
    B = [lbp(sig((0,) * (u + 1), i + 1), F[i], i + 1) for i in xrange(len(F))]
    B.sort(key=lambda f: O(sdp_LM(Polyn(f), u)), reverse=True)

    # critical pairs
    CP = [critical_pair(B[i], B[j], u, O, K) for i in xrange(len(B)) for j in xrange(i + 1, len(B))]
    CP.sort(key=lambda cp: cp_key(cp, O), reverse=True)

    k = len(B)

    reductions_to_zero = 0

    while len(CP):
        cp = CP.pop()

        # discard redundant critical pairs:
        if is_rewritable_or_comparable(cp[0], Num(cp[2]), B, u, K):
            continue
        if is_rewritable_or_comparable(cp[3], Num(cp[5]), B, u, K):
            continue

        s = s_poly(cp, u, O, K)

        p = f5_reduce(s, B, u, O, K)

        p = lbp(Sign(p), sdp_monic(Polyn(p), K), k + 1)

        if Polyn(p) != []:
            # remove old critical pairs, that become redundant when adding p:
            indices = []
            for i, cp in enumerate(CP):
                if is_rewritable_or_comparable(cp[0], Num(cp[2]), [p], u, K):
                    indices.append(i)
                elif is_rewritable_or_comparable(cp[3], Num(cp[5]), [p], u, K):
                    indices.append(i)

            for i in reversed(indices):
                del CP[i]

            # only add new critical pairs that are not made redundant by p:
            for g in B:
                if Polyn(g) != []:
                    cp = critical_pair(p, g, u, O, K)
                    if is_rewritable_or_comparable(cp[0], Num(cp[2]), [p], u, K):
                        continue
                    elif is_rewritable_or_comparable(cp[3], Num(cp[5]), [p], u, K):
                        continue

                    CP.append(cp)

            # sort (other sorting methods/selection strategies were not as successful)
            CP.sort(key=lambda cp: cp_key(cp, O), reverse=True)

            # insert p into B:
            m = sdp_LM(Polyn(p), u)
            if O(m) <= O(sdp_LM(Polyn(B[-1]), u)):
                B.append(p)
            else:
                for i, q in enumerate(B):
                    if O(m) > O(sdp_LM(Polyn(q), u)):
                        B.insert(i, p)
                        break

            k += 1

            #print(len(B), len(CP), "%d critical pairs removed" % len(indices))
        else:
            reductions_to_zero += 1

    if verbose:
        print("%d reductions to zero" % reductions_to_zero)

    # reduce Groebner basis:
    H = [sdp_monic(Polyn(g), K) for g in B]
    H = red_groebner(H, u, O, K)

    return sorted(H, key=lambda f: O(sdp_LM(f, u)), reverse=True)


def red_groebner(G, u, O, K):
    """
    Compute reduced Groebner basis, from BeckerWeispfenning93, p. 216

    Selects a subset of generators, that already generate the ideal
    and computes a reduced Groebner basis for them.
    """
    def reduction(P, u, O, K):
        """
        The actual reduction algorithm.
        """
        Q = []
        for i, p in enumerate(P):
            h = sdp_rem(p, P[:i] + P[i + 1:], u, O, K)
            if h != []:
                Q.append(h)

        return [sdp_monic(p, K) for p in Q]

    F = G
    H = []

    while F:
        f0 = F.pop()

        if not any(monomial_divides(sdp_LM(f0, u), sdp_LM(f, u)) for f in F + H):
            H.append(f0)

    # Becker, Weispfenning, p. 217: H is Groebner basis of the ideal generated by G.
    return reduction(H, u, O, K)


def is_groebner(G, u, O, K):
    """
    Check if G is a Groebner basis.
    """
    for i in xrange(len(G)):
        for j in xrange(i + 1, len(G)):
            s = sdp_spoly(G[i], G[j], u, O, K)
            s = sdp_rem(s, G, u, O, K)
            if s != []:
                return False

    return True


def is_minimal(G, u, O, K):
    """
    Checks if G is a minimal Groebner basis.
    """
    G.sort(key=lambda g: O(sdp_LM(g, u)))
    for i, g in enumerate(G):
        if sdp_LC(g, K) != K.one:
            return False

        for h in G[:i] + G[i + 1:]:
            if monomial_divides(sdp_LM(g, u), sdp_LM(h, u)):
                return False

    return True


def is_reduced(G, u, O, K):
    """
    Checks if G is a reduced Groebner basis.
    """
    G.sort(key=lambda g: O(sdp_LM(g, u)))
    for i, g in enumerate(G):
        if sdp_LC(g, K) != K.one:
            return False

        for term in g:
            for h in G[:i] + G[i + 1:]:
                if monomial_divides(term[0], sdp_LM(h, u)):
                    return False

    return True


def monomial_divides(m1, m2):
    """
    Returns True if m2 divides m1, False otherwise. Does not create
    the quotient. Does not check if both are have the same length.
    """
    for i in xrange(len(m1)):
        if m1[i] < m2[i]:
            return False

    return True

# FGLM

def matrix_fglm(F, u, O_from, O_to, K):
    """
    Converts the reduced Groebner basis ``F`` of a zero-dimensional
    ideal w.r.t. ``O_from`` to a reduced Groebner basis
    w.r.t. ``O_to``.

    References
    ==========

    J.C. Faugere, P. Gianni, D. Lazard, T. Mora (1994). Efficient
    Computation of Zero-dimensional Groebner Bases by Change of
    Ordering

    J.C. Faugere's lecture notes:
    http://www-salsa.lip6.fr/~jcf/Papers/2010_MPRI5e.pdf
    """
    old_basis = _basis(F, u, O_from, K)
    M = _representing_matrices(old_basis, F, u, O_from, K)

    # V contains the normalforms (wrt O_from) of S
    S = [(0,) * (u + 1)]
    V = [[K.one] + [K.zero] * (len(old_basis) - 1)]
    G = []

    L = [(i, 0) for i in xrange(u + 1)]  # (i, j) corresponds to x_i * S[j]
    L.sort(key=lambda (k, l): O_to(_incr_k(S[l], k)), reverse=True)
    t = L.pop()

    P = _identity_matrix(len(old_basis), K)

    while True:
        s = len(S)
        v = _matrix_mul(M[t[0]], V[t[1]], K)
        _lambda = _matrix_mul(P, v, K)

        if all(_lambda[i] == K.zero for i in xrange(s, len(old_basis))):
            # there is a linear combination of v by V

            lt = [(_incr_k(S[t[1]], t[0]), K.one)]
            rest = sdp_strip(sdp_sort([(S[i], _lambda[i]) for i in xrange(s)], O_to))
            g = sdp_sub(lt, rest, u, O_to, K)

            if g != []:
                G.append(g)

        else:
            # v is linearly independant from V
            P = _update(s, _lambda, P, K)
            S.append(_incr_k(S[t[1]], t[0]))
            V.append(v)

            L.extend([(i, s) for i in xrange(u + 1)])
            L = list(set(L))
            L.sort(key=lambda (k, l): O_to(_incr_k(S[l], k)), reverse=True)

        L = [(k, l) for (k, l) in L if \
            all(monomial_div(_incr_k(S[l], k), sdp_LM(g, u)) is None for g in G)]

        if not L:
            G = [ sdp_monic(g, K) for g in G ]
            return sorted(G, key=lambda g: O_to(sdp_LM(g, u)), reverse=True)

        t = L.pop()


def _incr_k(m, k):
    return tuple(list(m[:k]) + [m[k] + 1] + list(m[k + 1:]))


def _identity_matrix(n, K):
    M = [[K.zero] * n for _ in xrange(n)]

    for i in xrange(n):
        M[i][i] = K.one

    return M


def _matrix_mul(M, v, K):
    return [sum([row[i] * v[i] for i in xrange(len(v))]) for row in M]


def _update(s, _lambda, P, K):
    """
    Update ``P`` such that for the updated `P'` `P' v = e_{s}`.
    """
    k = min([j for j in xrange(s, len(_lambda)) if _lambda[j] != 0])

    for r in xrange(len(_lambda)):
        if r != k:
            P[r] = [P[r][j] - (P[k][j] * _lambda[r]) / _lambda[k] for j in xrange(len(P[r]))]

    P[k] = [P[k][j] / _lambda[k] for j in xrange(len(P[k]))]

    P[k], P[s] = P[s], P[k]

    return P


def _representing_matrices(basis, G, u, O, K):
    """
    Compute the matrices corresponding to the linear maps `m \mapsto
    x_i m` for all variables `x_i`.
    """
    def var(i):
        return tuple([0] * i + [1] + [0] * (u - i))

    def representing_matrix(m):
        M = [[K.zero] * len(basis) for _ in xrange(len(basis))]

        for i, v in enumerate(basis):
            r = sdp_rem([(monomial_mul(m, v), K.one)], G, u, O, K)

            for term in r:
                j = basis.index(term[0])
                M[j][i] = term[1]

        return M

    return [representing_matrix(var(i)) for i in xrange(u + 1)]


def _basis(G, u, O, K):
    """
    Computes a list of monomials which are not divisible by the leading
    monomials wrt to ``O`` of ``G``. These monomials are a basis of
    `K[X_1, \ldots, X_n]/(G)`.
    """
    leading_monomials = [sdp_LM(g, u) for g in G]
    candidates = [(0,) * (u + 1)]
    basis = []

    while candidates:
        t = candidates.pop()
        basis.append(t)

        new_candidates = [_incr_k(t, k) for k in xrange(u + 1) \
            if all(monomial_div(_incr_k(t, k), lmg) is None \
            for lmg in leading_monomials)]
        candidates.extend(new_candidates)
        candidates.sort(key=lambda m: O(m), reverse=True)

    basis = list(set(basis))

    return sorted(basis, key=lambda m: O(m))

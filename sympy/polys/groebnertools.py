"""Groebner bases algorithms. """

from sympy.core.compatibility import minkey, cmp

from sympy.polys.monomialtools import (
    monomial_mul,
    monomial_div,
    monomial_lcm,
    monomial_lex_key as O_lex,
    monomial_grlex_key as O_grlex,
    monomial_grevlex_key as O_grevlex,
)

from sympy.polys.distributedpolys import (
    sdp_LM,
    sdp_LT,
    sdp_mul_term,    
    sdp_sub,
    sdp_mul_term,
    sdp_monic,
    sdp_rem,    
)

from sympy.polys.polyerrors import (
    ExactQuotientFailed, DomainError,
)

from sympy.utilities import any, all
from operator import itemgetter

from sympy.polys.polyconfig import query

def sdp_groebner(f, u, O, K, gens='', verbose=False):
    """
    Wrapper around the (default) Buchberger and other algorithms
    for Groebner bases. The choice of algorithm can be changed via

    >>>> from sympy.polys.polyconfig import setup
    >>>> setup('GB_METHOD', 'method')

    where 'method' can be 'buchberger' or 'f5b'. If an unknown
    method is provided, the default Buchberger algorithm will be
    used.

    """
    if query('GB_METHOD') == 'buchberger':
        return buchberger(f, u, O, K, gens, verbose)
    #elif query('GB_METHOD') == 'f5b':
    #    return f5b(f, u, O, K, gens, verbose)
    else:
        return buchberger(f, u, O, K, gens, verbose)

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
        pr = minkey(P, key=lambda pair: O(monomial_lcm(sdp_LM(f[pair[0]], u), sdp_LM(f[pair[1]], u))))
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
        h = minkey([f[x] for x in F], key=lambda f: O(sdp_LM(f, u)))
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

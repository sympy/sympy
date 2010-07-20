"""Polynomial factorization routines in characteristic zero. """

from sympy.polys.galoistools import (
    gf_from_int_poly, gf_to_int_poly,
    gf_degree, gf_from_dict,
    gf_lshift, gf_add_mul, gf_mul,
    gf_div, gf_rem,
    gf_gcd, gf_gcdex,
    gf_sqf_p,
    gf_factor_sqf, gf_factor
)

from sympy.polys.densebasic import (
    dup_LC, dmp_LC, dmp_ground_LC,
    dup_TC, dmp_TC, dmp_ground_TC,
    dup_convert, dmp_convert,
    dup_degree, dmp_degree,
    dmp_degree_in, dmp_degree_list,
    dup_from_dict, dmp_from_dict,
    dmp_zero, dmp_zero_p,
    dmp_one, dmp_one_p,
    dmp_nest, dmp_raise,
    dup_strip, dmp_strip,
    dmp_ground,
    dup_inflate,
    dmp_exclude, dmp_include,
    dmp_inject, dmp_eject,
    dmp_terms_gcd
)

from sympy.polys.densearith import (
    dup_neg, dmp_neg,
    dup_add, dmp_add,
    dup_sub, dmp_sub,
    dup_mul, dmp_mul,
    dup_pow, dmp_pow,
    dup_div, dmp_div,
    dup_rem, dmp_rem,
    dup_exquo, dmp_exquo,
    dup_expand, dmp_expand,
    dup_add_mul, dmp_add_mul,
    dup_sub_mul, dmp_sub_mul,
    dup_max_norm, dmp_max_norm,
    dup_l1_norm, dmp_l1_norm,
    dup_mul_ground, dmp_mul_ground,
    dup_exquo_ground, dmp_exquo_ground
)

from sympy.polys.densetools import (
    dup_gcd, dmp_gcd,
    dup_sqf_p, dmp_sqf_p,
    dup_sqf_part, dmp_sqf_part,
    dup_trunc, dmp_ground_trunc,
    dup_content, dmp_ground_content,
    dup_monic, dmp_ground_monic,
    dup_primitive, dmp_primitive, dmp_ground_primitive,
    dup_clear_denoms, dmp_clear_denoms,
    dup_eval, dmp_eval_tail,
    dmp_eval_in, dmp_diff_eval_in,
    dup_inner_gcd, dmp_inner_gcd,
    dup_sqf_norm, dmp_sqf_norm,
    dup_compose, dmp_compose,
    dup_shift
)

from sympy.polys.polyutils import _sort_factors
from sympy.polys.polyconfig import query

from sympy.polys.polyerrors import (
    ExtraneousFactors, DomainError, EvaluationFailed
)

from sympy.ntheory import nextprime, isprime, factorint
from sympy.utilities import any, all, subsets, cythonized

from math import ceil, log
from random import randint

@cythonized("k")
def dup_trial_division(f, factors, K):
    """Determine multiplicities of factors using trial division. """
    result = []

    for factor in factors:
        k = 0

        while True:
            q, r = dup_div(f, factor, K)

            if not r:
                f, k = q, k+1
            else:
                break

        result.append((factor, k))

    return _sort_factors(result)

@cythonized("u,k")
def dmp_trial_division(f, factors, u, K):
    """Determine multiplicities of factors using trial division. """
    result = []

    for factor in factors:
        k = 0

        while True:
            q, r = dmp_div(f, factor, u, K)

            if dmp_zero_p(r, u):
                f, k = q, k+1
            else:
                break

        result.append((factor, k))

    return _sort_factors(result)

def dup_zz_mignotte_bound(f, K):
    """Mignotte bound for univariate polynomials in `K[x]`. """
    a = dup_max_norm(f, K)
    b = abs(dup_LC(f, K))
    n = dup_degree(f)

    return K.sqrt(n+1)*2**n*a*b

def dmp_zz_mignotte_bound(f, u, K):
    """Mignotte bound for multivariate polynomials in `K[X]`. """
    a = dmp_max_norm(f, u, K)
    b = abs(dmp_ground_LC(f, u, K))
    n = sum(dmp_degree_list(f, u))

    return K.sqrt(n+1)*2**n*a*b

def dup_zz_hensel_step(m, f, g, h, s, t, K):
    """One step in Hensel lifting in `Z[x]`.

       Given positive integer `m` and `Z[x]` polynomials `f`, `g`, `h`, `s`
       and `t` such that::

           f == g*h (mod m)
           s*g + t*h == 1 (mod m)

           lc(f) is not a zero divisor (mod m)
           lc(h) == 1

           deg(f) == deg(g) + deg(h)
           deg(s) < deg(h)
           deg(t) < deg(g)

       returns polynomials `G`, `H`, `S` and `T`, such that::

           f == G*H (mod m**2)
           S*G + T**H == 1 (mod m**2)

       References
       ==========

       .. [Gathen99] J. von zur Gathen, J. Gerhard, Modern Computer Algebra,
           First Edition, Cambridge University Press, 1999, pp. 418

    """
    M = m**2

    e = dup_sub_mul(f, g, h, K)
    e = dup_trunc(e, M, K)

    q, r = dup_div(dup_mul(s, e, K), h, K)

    q = dup_trunc(q, M, K)
    r = dup_trunc(r, M, K)

    u = dup_add(dup_mul(t, e, K), dup_mul(q, g, K), K)
    G = dup_trunc(dup_add(g, u, K), M, K)
    H = dup_trunc(dup_add(h, r, K), M, K)

    u = dup_add(dup_mul(s, G, K), dup_mul(t, H, K), K)
    b = dup_trunc(dup_sub(u, [K.one], K), M, K)

    c, d = dup_div(dup_mul(s, b, K), H, K)

    c = dup_trunc(c, M, K)
    d = dup_trunc(d, M, K)

    u = dup_add(dup_mul(t, b, K), dup_mul(c, G, K), K)
    S = dup_trunc(dup_sub(s, d, K), M, K)
    T = dup_trunc(dup_sub(t, u, K), M, K)

    return G, H, S, T

@cythonized("l,r,k,d")
def dup_zz_hensel_lift(p, f, f_list, l, K):
    """Multifactor Hensel lifting in `Z[x]`.

       Given a prime `p`, polynomial `f` over `Z[x]` such that `lc(f)`
       is a unit modulo `p`, monic pair-wise coprime polynomials `f_i`
       over `Z[x]` satisfying::

           f = lc(f) f_1 ... f_r (mod p)

       and a positive integer `l`, returns a list of monic polynomials
       `F_1`, `F_2`, ..., `F_r` satisfying::

          f = lc(f) F_1 ... F_r (mod p**l)

          F_i = f_i (mod p), i = 1..r

       References
       ==========

       .. [Gathen99] J. von zur Gathen, J. Gerhard, Modern Computer Algebra,
           First Edition, Cambridge University Press, 1999, pp. 424

    """
    r = len(f_list)
    lc = dup_LC(f, K)

    if r == 1:
        F = dup_mul_ground(f, K.gcdex(lc, p**l)[0], K)
        return [ dup_trunc(F, p**l, K) ]

    m = p
    k = r // 2
    d = int(ceil(log(l, 2)))

    g = gf_from_int_poly([lc], p)

    for f_i in f_list[:k]:
        g = gf_mul(g, gf_from_int_poly(f_i, p), p, K)

    h = gf_from_int_poly(f_list[k], p)

    for f_i in f_list[k+1:]:
        h = gf_mul(h, gf_from_int_poly(f_i, p), p, K)

    s, t, _ = gf_gcdex(g, h, p, K)

    g = gf_to_int_poly(g, p)
    h = gf_to_int_poly(h, p)
    s = gf_to_int_poly(s, p)
    t = gf_to_int_poly(t, p)

    for _ in range(1, d+1):
        (g, h, s, t), m = dup_zz_hensel_step(m, f, g, h, s, t, K), m**2

    return dup_zz_hensel_lift(p, g, f_list[:k], l, K) \
         + dup_zz_hensel_lift(p, h, f_list[k:], l, K)

@cythonized("l,s")
def dup_zz_zassenhaus(f, K):
    """Factor primitive square-free polynomials in `Z[x]`. """
    n = dup_degree(f)

    if n == 1:
        return [f]

    A = dup_max_norm(f, K)
    b = dup_LC(f, K)
    B = int(abs(K.sqrt(n+1)*2**n*A*b))
    C = int((n+1)**(2*n)*A**(2*n-1))
    gamma = int(ceil(2*log(C, 2)))
    bound = int(2*gamma*log(gamma))

    for p in xrange(3, bound+1):
        if not isprime(p) or b % p == 0:
            continue

        p = K.convert(p)

        F = gf_from_int_poly(f, p)

        if gf_sqf_p(F, p, K):
            break

    l = int(ceil(log(2*B + 1, p)))

    modular = []

    for ff in gf_factor_sqf(F, p, K)[1]:
        modular.append(gf_to_int_poly(ff, p))

    g = dup_zz_hensel_lift(p, f, modular, l, K)

    T = set(range(len(g)))
    factors, s = [], 1

    while 2*s <= len(T):
        for S in subsets(T, s):
            G, H = [b], [b]

            S = set(S)

            for i in S:
                G = dup_mul(G, g[i], K)
            for i in T-S:
                H = dup_mul(H, g[i], K)

            G = dup_trunc(G, p**l, K)
            H = dup_trunc(H, p**l, K)

            G_norm = dup_l1_norm(G, K)
            H_norm = dup_l1_norm(H, K)

            if G_norm*H_norm <= B:
                T = T - S

                G = dup_primitive(G, K)[1]
                f = dup_primitive(H, K)[1]

                factors.append(G)
                b = dup_LC(f, K)

                break
        else:
            s += 1

    return factors + [f]

def dup_zz_irreducible_p(f, K):
    """Test irreducibility using Eisenstein's criterion. """
    lc = dup_LC(f, K)
    tc = dup_TC(f, K)

    e_fc = dup_content(f[1:], K)

    if e_fc:
        e_ff = factorint(int(e_fc))

        for p in e_ff.iterkeys():
            if (lc % p) and (tc % p**2):
                return True

@cythonized("n,p,k")
def dup_zz_cyclotomic_poly(n, K):
    """Efficiently generate n-th cyclotomic polnomial. """
    h = [K.one,-K.one]

    for p, k in factorint(n).iteritems():
        h = dup_exquo(dup_inflate(h, p, K), h, K)
        h = dup_inflate(h, p**(k-1), K)

    return h

@cythonized("n,p,k,i")
def _dup_cyclotomic_decompose(n, K):
    H = [[K.one,-K.one]]

    for p, k in factorint(n).iteritems():
        Q = [ dup_exquo(dup_inflate(h, p, K), h, K) for h in H ]
        H.extend(Q)

        for i in xrange(1, k):
            Q = [ dup_inflate(q, p, K) for q in Q ]
            H.extend(Q)

    return H

@cythonized("n")
def dup_zz_cyclotomic_factor(f, K):
    """Efficiently factor polynomials `x**n - 1` and `x**n + 1` in `Z[x]`.

       Given a univariate polynomial `f` in `Z[x]` returns a list of factors
       of `f`, provided that `f` is in the form `x**n - 1` or `x**n + 1` for
       `n >= 1`. Otherwise returns None.

       Factorization is performed using using cyclotomic decomposition of `f`,
       which makes this method much faster that any other direct factorization
       approach (e.g. Zassenhaus's).

       References
       ==========

       .. [Weisstein09] Eric W. Weisstein, Cyclotomic Polynomial, From MathWorld - A
           Wolfram Web Resource, http://mathworld.wolfram.com/CyclotomicPolynomial.html

    """
    lc_f, tc_f = dup_LC(f, K), dup_TC(f, K)

    if dup_degree(f) <= 0:
        return None

    if lc_f != 1 or tc_f not in [-1, 1]:
        return None

    if any([ bool(cf) for cf in f[1:-1] ]):
        return None

    n = dup_degree(f)
    F = _dup_cyclotomic_decompose(n, K)

    if not K.is_one(tc_f):
        return F
    else:
        H = []

        for h in _dup_cyclotomic_decompose(2*n, K):
            if h not in F:
                H.append(h)

        return H

@cythonized("n")
def dup_zz_factor_sqf(f, K):
    """Factor square-free (non-primitive) polyomials in `Z[x]`. """
    cont, g = dup_primitive(f, K)

    n = dup_degree(g)

    if dup_LC(g, K) < 0:
        cont, g = -cont, dup_neg(g, K)

    if n <= 0:
        return cont, []

    if n == 1 or dup_zz_irreducible_p(g, K):
        return cont, [(g, 1)]

    factors = None

    if query('USE_CYCLOTOMIC_FACTOR'):
        factors = dup_zz_cyclotomic_factor(g, K)

    if factors is None:
        factors = dup_zz_zassenhaus(g, K)

    return cont, _sort_factors(factors, multiple=False)

@cythonized("n,k")
def dup_zz_factor(f, K):
    """Factor (non square-free) polynomials in `Z[x]`.

       Given a univariate polynomial `f` in `Z[x]` computes its complete
       factorization `f_1, ..., f_n` into irreducibles over integers::

                   f = content(f) f_1**k_1 ... f_n**k_n

       The factorization is computed by reducing the input polynomial
       into a primitive square-free polynomial and factoring it using
       Zassenhaus algorithm. Trial division is used to recover the
       multiplicities of factors.

       The result is returned as a tuple consisting of::

                 (content(f), [(f_1, k_1), ..., (f_n, k_n))

       Consider polynomial `f = 2*x**4 - 2`::

           >>> from sympy.polys.factortools import dup_zz_factor
           >>> from sympy.polys.domains import ZZ

           >>> dup_zz_factor([2, 0, 0, 0, -2], ZZ)
           (2, [([1, -1], 1), ([1, 1], 1), ([1, 0, 1], 1)])

       In result we got the following factorization::

                    f = 2 (x - 1) (x + 1) (x**2 + 1)

       Note that this is a complete factorization over integers,
       however over Gaussian integers we can factor the last term.

       By default, polynomials `x**n - 1` and `x**n + 1` are factored
       using cyclotomic decomposition to speedup computations. To
       disable this behaviour set cyclotomic=False.

       References
       ==========

       .. [Gathen99] J. von zur Gathen, J. Gerhard, Modern Computer Algebra,
           First Edition, Cambridge University Press, 1999, pp. 427
    """
    cont, g = dup_primitive(f, K)

    n = dup_degree(g)

    if dup_LC(g, K) < 0:
        cont, g = -cont, dup_neg(g, K)

    if n <= 0:
        return cont, []

    if n == 1 or dup_zz_irreducible_p(g, K):
        return cont, [(g, 1)]

    g = dup_sqf_part(g, K)
    H, factors = None, []

    if query('USE_CYCLOTOMIC_FACTOR'):
        H = dup_zz_cyclotomic_factor(g, K)

    if H is None:
        H = dup_zz_zassenhaus(g, K)

    for h in H:
        k = 0

        while True:
            q, r = dup_div(f, h, K)

            if not r:
                f, k = q, k+1
            else:
                break

        factors.append((h, k))

    return cont, _sort_factors(factors)

def dmp_zz_wang_non_divisors(E, cs, ct, K):
    """Wang/EEZ: Compute a set of valid divisors.  """
    result = [ cs*ct ]

    for q in E:
        q = abs(q)

        for r in reversed(result):
            while r != 1:
                r = K.gcd(r, q)
                q = q // r

            if K.is_one(q):
                return None

        result.append(q)

    return result[1:]

@cythonized("u,v")
def dmp_zz_wang_test_points(f, T, ct, A, u, K):
    """Wang/EEZ: Test evaluation points for suitability. """
    if not dmp_eval_tail(dmp_LC(f, K), A, u-1, K):
        raise EvaluationFailed('no luck')

    g = dmp_eval_tail(f, A, u, K)

    if not dup_sqf_p(g, K):
        raise EvaluationFailed('no luck')

    c, h = dup_primitive(g, K)

    if K.is_negative(dup_LC(h, K)):
        c, h = -c, dup_neg(h, K)

    v = u-1

    E = [ dmp_eval_tail(t, A, v, K) for t, _ in T ]
    D = dmp_zz_wang_non_divisors(E, c, ct, K)

    if D is not None:
        return c, h, E
    else:
        raise EvaluationFailed('no luck')

@cythonized("u,v,i,j,k")
def dmp_zz_wang_lead_coeffs(f, T, cs, E, H, A, u, K):
    """Wang/EEZ: Compute correct leading coefficients. """
    C, J, v = [], [0]*len(E), u-1

    for h in H:
        c = dmp_one(v, K)
        d = dup_LC(h, K)*cs

        for i in reversed(xrange(len(E))):
            k, e, (t, _) = 0, E[i], T[i]

            while not (d % e):
                d, k = d//e, k+1

            if k != 0:
                c, J[i] = dmp_mul(c, dmp_pow(t, k, v, K), v, K), 1

        C.append(c)

    if any([ not j for j in J ]):
        raise ExtraneousFactors # pragma: no cover

    CC, HH = [], []

    for c, h in zip(C, H):
        d = dmp_eval_tail(c, A, v, K)
        lc = dup_LC(h, K)

        if K.is_one(cs):
            cc = lc//d
        else:
            g = K.gcd(lc, d)
            d, cc = d//g, lc//g
            h, cs = dup_mul_ground(h, d, K), cs//d

        c = dmp_mul_ground(c, cc, v, K)

        CC.append(c)
        HH.append(h)

    if K.is_one(cs):
        return f, HH, CC

    CCC, HHH = [], []

    for c, h in zip(CC, HH):
        CCC.append(dmp_mul_ground(c, cs, v, K))
        HHH.append(dmp_mul_ground(h, cs, 0, K))

    f = dmp_mul_ground(f, cs**(len(H)-1), u, K)

    return f, HHH, CCC

@cythonized("m")
def dup_zz_diophantine(F, m, p, K):
    """Wang/EEZ: Solve univariate Diophantine equations. """
    if len(F) == 2:
        a, b = F

        f = gf_from_int_poly(a, p)
        g = gf_from_int_poly(b, p)

        s, t, G = gf_gcdex(g, f, p, K)

        s = gf_lshift(s, m, K)
        t = gf_lshift(t, m, K)

        q, s = gf_div(s, f, p, K)

        t = gf_add_mul(t, q, g, p, K)

        s = gf_to_int_poly(s, p)
        t = gf_to_int_poly(t, p)

        result = [s, t]
    else:
        G = [F[-1]]

        for f in reversed(F[1:-1]):
            G.insert(0, dup_mul(f, G[0], K))

        S, T = [], [[1]]

        for f, g in zip(F, G):
            t, s = dmp_zz_diophantine([g, f], T[-1], [], 0, p, 1, K)
            T.append(t)
            S.append(s)

        result, S = [], S + [T[-1]]

        for s, f in zip(S, F):
            s = gf_from_int_poly(s, p)
            f = gf_from_int_poly(f, p)

            r = gf_rem(gf_lshift(s, m, K), f, p, K)
            s = gf_to_int_poly(r, p)

            result.append(s)

    return result

@cythonized("u,v,d,n,i,j,k")
def dmp_zz_diophantine(F, c, A, d, p, u, K):
    """Wang/EEZ: Solve multivariate Diophantine equations. """
    if not A:
        S = [ [] for _ in F ]
        n = dup_degree(c)

        for i, coeff in enumerate(c):
            if not coeff:
                continue

            T = dup_zz_diophantine(F, n-i, p, K)

            for j, (s, t) in enumerate(zip(S, T)):
                t = dup_mul_ground(t, coeff, K)
                S[j] = dup_trunc(dup_add(s, t, K), p, K)
    else:
        n = len(A)
        e = dmp_expand(F, u, K)

        a, A = A[-1], A[:-1]
        B, G = [], []

        for f in F:
            B.append(dmp_exquo(e, f, u, K))
            G.append(dmp_eval_in(f, a, n, u, K))

        C = dmp_eval_in(c, a, n, u, K)

        v = u - 1

        S = dmp_zz_diophantine(G, C, A, d, p, v, K)
        S = [ dmp_raise(s, 1, v, K) for s in S ]

        for s, b in zip(S, B):
            c = dmp_sub_mul(c, s, b, u, K)

        c = dmp_ground_trunc(c, p, u, K)

        m = dmp_nest([K.one, -a], n, K)
        M = dmp_one(n, K)

        for k in xrange(0, d):
            if dmp_zero_p(c, u):
                break

            M = dmp_mul(M, m, u, K)
            C = dmp_diff_eval_in(c, k+1, a, n, u, K)

            if not dmp_zero_p(C, v):
                C = dmp_exquo_ground(C, K.factorial(k+1), v, K)
                T = dmp_zz_diophantine(G, C, A, d, p, v, K)

                for i, t in enumerate(T):
                    T[i] = dmp_mul(dmp_raise(t, 1, v, K), M, u, K)

                for i, (s, t) in enumerate(zip(S, T)):
                    S[i] = dmp_add(s, t, u, K)

                for t, b in zip(T, B):
                    c = dmp_sub_mul(c, t, b, u, K)

                c = dmp_ground_trunc(c, p, u, K)

        S = [ dmp_ground_trunc(s, p, u, K) for s in S ]

    return S

@cythonized("u,v,d,dj,n,i,j,k,w")
def dmp_zz_wang_hensel_lifting(f, H, LC, A, p, u, K):
    """Wang/EEZ: Parallel Hensel lifting algorithm. """
    S, n, v = [f], len(A), u-1

    H = list(H)

    for i, a in enumerate(reversed(A[1:])):
        s = dmp_eval_in(S[0], a, n-i, u-i, K)
        S.insert(0, dmp_ground_trunc(s, p, v-i, K))

    d = max(dmp_degree_list(f, u)[1:])

    for j, s, a in zip(xrange(2, n+2), S, A):
        G, w = list(H), j-1

        I, J = A[:j-2], A[j-1:]

        for i, (h, lc) in enumerate(zip(H, LC)):
            lc = dmp_ground_trunc(dmp_eval_tail(lc, J, v, K), p, w-1, K)
            H[i] = [lc] + dmp_raise(h[1:], 1, w-1, K)

        m = dmp_nest([K.one, -a], w, K)
        M = dmp_one(w, K)

        c = dmp_sub(s, dmp_expand(H, w, K), w, K)

        dj = dmp_degree_in(s, w, w)

        for k in xrange(0, dj):
            if dmp_zero_p(c, w):
                break

            M = dmp_mul(M, m, w, K)
            C = dmp_diff_eval_in(c, k+1, a, w, w, K)

            if not dmp_zero_p(C, w-1):
                C = dmp_exquo_ground(C, K.factorial(k+1), w-1, K)
                T = dmp_zz_diophantine(G, C, I, d, p, w-1, K)

                for i, (h, t) in enumerate(zip(H, T)):
                    h = dmp_add_mul(h, dmp_raise(t, 1, w-1, K), M, w, K)
                    H[i] = dmp_ground_trunc(h, p, w, K)

                h = dmp_sub(s, dmp_expand(H, w, K), w, K)
                c = dmp_ground_trunc(h, p, w, K)

    if dmp_expand(H, u, K) != f:
        raise ExtraneousFactors # pragma: no cover
    else:
        return H

@cythonized("u,mod,i,j,s_arg,negative")
def dmp_zz_wang(f, u, K, mod=None):
    """Factor primitive square-free polynomials in `Z[X]`.

       Given a multivariate polynomial `f` in `Z[x_1,...,x_n]`, which
       is primitive and square-free in `x_1`, computes factorization
       of `f` into irreducibles over integers.

       The procedure is based on Wang's Enhanced Extended Zassenhaus
       algorithm. The algorithm works by viewing `f` as a univariate
       polynomial in `Z[x_2,...,x_n][x_1]`, for which an evaluation
       mapping is computed::

                         x_2 -> a_2, ..., x_n -> a_n

       where `a_i`, for `i = 2, ..., n`, are carefully chosen integers.
       The mapping is used to transform `f` into a univariate polynomial
       in `Z[x_1]`, which can be factored efficiently using Zassenhaus
       algorithm. The last step is to lift univariate factors to obtain
       true multivariate factors. For this purpose a parallel Hensel
       lifting procedure is used.

       References
       ==========

       .. [Wang78] P. S. Wang, An Improved Multivariate Polynomial Factoring
           Algorithm, Math. of Computation 32, 1978, pp. 1215--1231

       .. [Geddes92] K. Geddes, S. R. Czapor, G. Labahn, Algorithms for
           Computer Algebra, Springer, 1992, pp. 264--272
    """
    ct, T = dmp_zz_factor(dmp_LC(f, K), u-1, K)

    b = dmp_zz_mignotte_bound(f, u, K)
    p = K(nextprime(b))

    if mod is None:
        if u == 1:
            mod = 2
        else:
            mod = 1

    history, configs, A, r = set([]), [], [K.zero]*u, None

    try:
        cs, s, E = dmp_zz_wang_test_points(f, T, ct, A, u, K)

        _, H = dup_zz_factor_sqf(s, K)

        r = len(H)

        if r == 1:
            return [f]

        bad_points = set([tuple(A)])
        configs = [(s, cs, E, H, A)]
    except EvaluationFailed:
        pass

    eez_num_configs = query('EEZ_NUMBER_OF_CONFIGS')
    eez_num_tries = query('EEZ_NUMBER_OF_TRIES')
    eez_mod_step = query('EEZ_MODULUS_STEP')

    while len(configs) < eez_num_configs:
        for _ in xrange(eez_num_tries):
            A = [ K(randint(-mod, mod)) for _ in xrange(u) ]

            if tuple(A) not in history:
                history.add(tuple(A))
            else:
                continue

            try:
                cs, s, E = dmp_zz_wang_test_points(f, T, ct, A, u, K)
            except EvaluationFailed:
                continue

            _, H = dup_zz_factor_sqf(s, K)

            rr = len(H)

            if r is not None:
                if rr != r: # pragma: no cover
                    if rr < r:
                        configs, r = [], rr
                    else:
                        continue
            else:
                r = rr

            if r == 1:
                return [f]

            configs.append((s, cs, E, H, A))

            if len(configs) == eez_num_configs:
                break
        else:
            mod += eez_mod_step

    s_norm, s_arg, i = None, 0, 0

    for s, _, _, _, _ in configs:
        _s_norm = dup_max_norm(s, K)

        if s_norm is not None:
            if _s_norm < s_norm:
                s_norm = _s_norm
                s_arg = i
        else:
            s_norm = _s_norm

        i += 1

    _, cs, E, H, A = configs[s_arg]

    try:
        f, H, LC = dmp_zz_wang_lead_coeffs(f, T, cs, E, H, A, u, K)
        factors = dmp_zz_wang_hensel_lifting(f, H, LC, A, p, u, K)
    except ExtraneousFactors: # pragma: no cover
        if query('EEZ_RESTART_IF_NEEDED'):
            return dmp_zz_wang(f, u, K, mod+1)
        else:
            raise ExtraneousFactors("we need to restart algorithm with better parameters")

    negative, result = 0, []

    for f in factors:
        _, f = dmp_ground_primitive(f, u, K)

        if K.is_negative(dmp_ground_LC(f, u, K)):
            f = dmp_neg(f, u, K)

        result.append(f)

    return result

@cythonized("u,d,k")
def dmp_zz_factor(f, u, K):
    """Factor (non square-free) polynomials in `Z[X]`.

       Given a multivariate polynomial `f` in `Z[x]` computes its complete
       factorization `f_1, ..., f_n` into irreducibles over integers::

                    f = content(f) f_1**k_1 ... f_n**k_n

       The factorization is computed by reducing the input polynomial
       into a primitive square-free polynomial and factoring it using
       Enhanced Extended Zassenhaus (EEZ) algorithm. Trial division
       is used to recover the multiplicities of factors.

       The result is returned as a tuple consisting of::

                (content(f), [(f_1, k_1), ..., (f_n, k_n))

       Consider polynomial `f = 2*(x**2 - y**2)`::

           >>> from sympy.polys.factortools import dmp_zz_factor
           >>> from sympy.polys.domains import ZZ

           >>> dmp_zz_factor([[2], [], [-2, 0, 0]], 1, ZZ)
           (2, [([[1], [-1, 0]], 1), ([[1], [1, 0]], 1)])

       In result we got the following factorization::

                       f = 2 (x - y) (x + y)

       References
       ==========

       .. [Gathen99] J. von zur Gathen, J. Gerhard, Modern Computer Algebra,
           First Edition, Cambridge University Press, 1999, pp. 427

    """
    if not u:
        return dup_zz_factor(f, K)

    if dmp_zero_p(f, u):
        return K.zero, []

    cont, g = dmp_ground_primitive(f, u, K)

    if dmp_ground_LC(g, u, K) < 0:
        cont, g = -cont, dmp_neg(g, u, K)

    if all([ d <= 0 for d in dmp_degree_list(g, u) ]):
        return cont, []

    G, g = dmp_primitive(g, u, K)

    factors = []

    if dmp_degree(g, u) > 0:
        g = dmp_sqf_part(g, u, K)
        H = dmp_zz_wang(g, u, K)

        for h in H:
            k = 0

            while True:
                q, r = dmp_div(f, h, u, K)

                if dmp_zero_p(r, u):
                    f, k = q, k+1
                else:
                    break

            factors.append((h, k))

    for g, k in dmp_zz_factor(G, u-1, K)[1]:
        factors.insert(0, ([g], k))

    return cont, _sort_factors(factors)

def dup_ext_factor(f, K):
    """Factor univariate polynomials over algebraic number fields. """
    n, lc = dup_degree(f), dup_LC(f, K)

    f = dup_monic(f, K)

    if n <= 0:
        return lc, []
    if n == 1:
        return lc, [(f, 1)]

    f, F = dup_sqf_part(f, K), f
    s, g, r = dup_sqf_norm(f, K)

    factors = dup_factor_list_include(r, K.dom)

    if len(factors) == 1:
        return lc, [(f, n//dup_degree(f))]

    H = s*K.unit

    for i, (factor, _) in enumerate(factors):
        h = dup_convert(factor, K.dom, K)
        h, _, g = dup_inner_gcd(h, g, K)
        h = dup_shift(h, H, K)
        factors[i] = h

    factors = dup_trial_division(F, factors, K)

    return lc, factors

@cythonized("u")
def dmp_ext_factor(f, u, K):
    """Factor multivariate polynomials over algebraic number fields. """
    if not u:
        return dup_ext_factor(f, K)

    lc = dmp_ground_LC(f, u, K)
    f = dmp_ground_monic(f, u, K)

    if all([ d <= 0 for d in dmp_degree_list(f, u) ]):
        return lc, []

    f, F = dmp_sqf_part(f, u, K), f
    s, g, r = dmp_sqf_norm(f, u, K)

    factors = dmp_factor_list_include(r, u, K.dom)

    if len(factors) == 1:
        coeff, factors = lc, [f]
    else:
        H = dmp_raise([K.one, s*K.unit], u, 0, K)

        for i, (factor, _) in enumerate(factors):
            h = dmp_convert(factor, u, K.dom, K)
            h, _, g = dmp_inner_gcd(h, g, u, K)
            h = dmp_compose(h, H, u, K)
            factors[i] = h

    return lc, dmp_trial_division(F, factors, u, K)

@cythonized("i")
def dup_gf_factor(f, K):
    """Factor univariate polynomials over finite fields. """
    f = dup_convert(f, K, K.dom)

    coeff, factors = gf_factor(f, K.mod, K.dom)

    for i, (f, k) in enumerate(factors):
        factors[i] = (dup_convert(f, K.dom, K), k)

    return K.convert(coeff, K.dom), factors

def dmp_gf_factor(f, u, K):
    """Factor multivariate polynomials over finite fields. """
    raise DomainError('multivariate polynomials over %s' % K)

@cythonized("i,k,u")
def dup_factor_list(f, K0):
    """Factor polynomials into irreducibles in `K[x]`. """
    if not K0.has_CharacteristicZero:
        coeff, factors = dup_gf_factor(f, K0)
    elif K0.is_Algebraic:
        coeff, factors = dup_ext_factor(f, K0)
    else:
        if not K0.is_Exact:
            K0_inexact, K0 = K0, K0.get_exact()
            f = dup_convert(f, K0_inexact, K0)
        else:
            K0_inexact = None

        if K0.has_Field:
            K = K0.get_ring()

            denom, f = dup_clear_denoms(f, K0, K)
            f = dup_convert(f, K0, K)
        else:
            K = K0

        if K.is_ZZ:
            coeff, factors = dup_zz_factor(f, K)
        elif K.is_Poly:
            f, u = dmp_inject(f, 0, K)

            coeff, factors = dmp_factor_list(f, u, K.dom)

            for i, (f, k) in enumerate(factors):
                factors[i] = (dmp_eject(f, u, K), k)

            coeff = K.convert(coeff, K.dom)
        else: # pragma: no cover
            raise DomainError('factorization not supported over %s' % K0)

        if K0.has_Field:
            for i, (f, k) in enumerate(factors):
                factors[i] = (dup_convert(f, K, K0), k)

            coeff = K0.convert(coeff, K)
            denom = K0.convert(denom, K)

            coeff = K0.quo(coeff, denom)

        if K0_inexact is not None:
            for i, (f, k) in enumerate(factors):
                factors[i] = (dup_convert(f, K0, K0_inexact), k)

            coeff = K0_inexact.convert(coeff, K0)

    return coeff, factors

def dup_factor_list_include(f, K):
    """Factor polynomials into irreducibles in `K[x]`. """
    coeff, factors = dup_factor_list(f, K)

    if not factors:
        return [(dup_strip([coeff]), 1)]
    else:
        g = dup_mul_ground(factors[0][0], coeff, K)
        return [(g, factors[0][1])] + factors[1:]

@cythonized("u,v,i,k")
def _dmp_inner_factor(f, u, K):
    """Simplify factorization in `Z[X]` as much as possible. """
    gcd, f = dmp_terms_gcd(f, u, K)
    J, f, v = dmp_exclude(f, u, K)

    coeff, factors = dmp_zz_factor(f, v, K)

    for i, (f, k) in enumerate(factors):
        factors[i] = (dmp_include(f, J, v, K), k)

    for i, g in enumerate(reversed(gcd)):
        if not g:
            continue

        term = {(0,)*(u-i) + (1,) + (0,)*i: K.one}
        factors.insert(0, (dmp_from_dict(term, u, K), g))

    return coeff, factors

@cythonized("u,v,i,k")
def dmp_factor_list(f, u, K0):
    """Factor polynomials into irreducibles in `K[X]`. """
    if not u:
        return dup_factor_list(f, K0)

    if not K0.has_CharacteristicZero: # pragma: no cover
        coeff, factors = dmp_gf_factor(f, u, K0)
    elif K0.is_Algebraic:
        coeff, factors = dmp_ext_factor(f, u, K0)
    else:
        if not K0.is_Exact:
            K0_inexact, K0 = K0, K0.get_exact()
            f = dmp_convert(f, u, K0_inexact, K0)
        else:
            K0_inexact = None

        if K0.has_Field:
            K = K0.get_ring()

            denom, f = dmp_clear_denoms(f, u, K0, K)
            f = dmp_convert(f, u, K0, K)
        else:
            K = K0

        if K.is_ZZ:
            coeff, factors = _dmp_inner_factor(f, u, K)
        elif K.is_Poly:
            f, v = dmp_inject(f, u, K)

            coeff, factors = dmp_factor_list(f, v, K.dom)

            for i, (f, k) in enumerate(factors):
                factors[i] = (dmp_eject(f, v, K), k)

            coeff = K.convert(coeff, K.dom)
        else: # pragma: no cover
            raise DomainError('factorization not supported over %s' % K0)

        if K0.has_Field:
            for i, (f, k) in enumerate(factors):
                factors[i] = (dmp_convert(f, u, K, K0), k)

            coeff = K0.convert(coeff, K)
            denom = K0.convert(denom, K)

            coeff = K0.quo(coeff, denom)

        if K0_inexact is not None:
            for i, (f, k) in enumerate(factors):
                factors[i] = (dmp_convert(f, u, K0, K0_inexact), k)

            coeff = K0_inexact.convert(coeff, K0)

    return coeff, factors

@cythonized("u")
def dmp_factor_list_include(f, u, K):
    """Factor polynomials into irreducibles in `K[X]`. """
    if not u:
        return dup_factor_list_include(f, K)

    coeff, factors = dmp_factor_list(f, u, K)

    if not factors:
        return [(dmp_ground(coeff, u), 1)]
    else:
        g = dmp_mul_ground(factors[0][0], coeff, u, K)
        return [(g, factors[0][1])] + factors[1:]


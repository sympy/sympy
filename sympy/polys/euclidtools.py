"""Euclidean algorithms, GCDs, LCMs and polynomial remainder sequences. """


from collections import defaultdict
from itertools import chain, compress

from sympy.ntheory import nextprime
from sympy.polys.densearith import (
    dup_sub_mul,
    dup_neg, dmp_neg,
    dmp_add,
    dmp_sub,
    dup_mul, dmp_mul,
    dmp_pow,
    dup_div, dmp_div,
    dup_rem,
    dup_quo, dmp_quo,
    dup_prem, dmp_prem,
    dup_mul_ground, dmp_mul_ground,
    dmp_mul_term,
    dup_quo_ground, dmp_quo_ground,
    dup_max_norm, dmp_max_norm)
from sympy.polys.densebasic import (
    dup_strip, dmp_raise,
    dmp_zero, dmp_one, dmp_ground,
    dmp_one_p, dmp_zero_p,
    dmp_zeros,
    dup_degree, dmp_degree, dmp_degree_in,
    dup_LC, dmp_LC, dmp_ground_LC,
    dmp_multi_deflate, dmp_inflate,
    dup_convert, dmp_convert,
    dmp_apply_pairs)
from sympy.polys.densetools import (
    dup_clear_denoms, dmp_clear_denoms,
    dup_diff, dmp_diff,
    dup_eval, dmp_eval, dmp_eval_in,
    dup_trunc, dmp_ground_trunc,
    dup_monic, dmp_ground_monic,
    dup_primitive, dmp_ground_primitive,
    dup_extract, dmp_ground_extract)
from sympy.polys.galoistools import (
    gf_int, gf_crt)
from sympy.polys.polyconfig import query
from sympy.polys.polyerrors import (
    MultivariatePolynomialError,
    HeuristicGCDFailed,
    HomomorphismFailed,
    NotInvertible,
    DomainError)




def dup_half_gcdex(f, g, K):
    """
    Half extended Euclidean algorithm in `F[x]`.

    Returns ``(s, h)`` such that ``h = gcd(f, g)`` and ``s*f = h (mod g)``.

    Examples
    ========

    >>> from sympy.polys import ring, QQ
    >>> R, x = ring("x", QQ)

    >>> f = x**4 - 2*x**3 - 6*x**2 + 12*x + 15
    >>> g = x**3 + x**2 - 4*x - 4

    >>> R.dup_half_gcdex(f, g)
    (-1/5*x + 3/5, x + 1)

    """
    if not K.is_Field:
        raise DomainError("can't compute half extended GCD over %s" % K)

    a, b = [K.one], []

    while g:
        q, r = dup_div(f, g, K)
        f, g = g, r
        a, b = b, dup_sub_mul(a, q, b, K)

    a = dup_quo_ground(a, dup_LC(f, K), K)
    f = dup_monic(f, K)

    return a, f


def dmp_half_gcdex(f, g, u, K):
    """
    Half extended Euclidean algorithm in `F[X]`.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    """
    if not u:
        return dup_half_gcdex(f, g, K)
    else:
        raise MultivariatePolynomialError(f, g)


def dup_gcdex(f, g, K):
    """
    Extended Euclidean algorithm in `F[x]`.

    Returns ``(s, t, h)`` such that ``h = gcd(f, g)`` and ``s*f + t*g = h``.

    Examples
    ========

    >>> from sympy.polys import ring, QQ
    >>> R, x = ring("x", QQ)

    >>> f = x**4 - 2*x**3 - 6*x**2 + 12*x + 15
    >>> g = x**3 + x**2 - 4*x - 4

    >>> R.dup_gcdex(f, g)
    (-1/5*x + 3/5, 1/5*x**2 - 6/5*x + 2, x + 1)

    """
    s, h = dup_half_gcdex(f, g, K)

    F = dup_sub_mul(h, s, f, K)
    t = dup_quo(F, g, K)

    return s, t, h


def dmp_gcdex(f, g, u, K):
    """
    Extended Euclidean algorithm in `F[X]`.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    """
    if not u:
        return dup_gcdex(f, g, K)
    else:
        raise MultivariatePolynomialError(f, g)


def dup_invert(f, g, K):
    """
    Compute multiplicative inverse of `f` modulo `g` in `F[x]`.

    Examples
    ========

    >>> from sympy.polys import ring, QQ
    >>> R, x = ring("x", QQ)

    >>> f = x**2 - 1
    >>> g = 2*x - 1
    >>> h = x - 1

    >>> R.dup_invert(f, g)
    -4/3

    >>> R.dup_invert(f, h)
    Traceback (most recent call last):
    ...
    NotInvertible: zero divisor

    """
    s, h = dup_half_gcdex(f, g, K)

    if h == [K.one]:
        return dup_rem(s, g, K)
    else:
        raise NotInvertible("zero divisor")


def dmp_invert(f, g, u, K):
    """
    Compute multiplicative inverse of `f` modulo `g` in `F[X]`.

    Examples
    ========

    >>> from sympy.polys import ring, QQ
    >>> R, x = ring("x", QQ)

    """
    if not u:
        return dup_invert(f, g, K)
    else:
        raise MultivariatePolynomialError(f, g)


def dup_euclidean_prs(f, g, K):
    """
    Euclidean polynomial remainder sequence (PRS) in `K[x]`.

    Examples
    ========

    >>> from sympy.polys import ring, QQ
    >>> R, x = ring("x", QQ)

    >>> f = x**8 + x**6 - 3*x**4 - 3*x**3 + 8*x**2 + 2*x - 5
    >>> g = 3*x**6 + 5*x**4 - 4*x**2 - 9*x + 21

    >>> prs = R.dup_euclidean_prs(f, g)

    >>> prs[0]
    x**8 + x**6 - 3*x**4 - 3*x**3 + 8*x**2 + 2*x - 5
    >>> prs[1]
    3*x**6 + 5*x**4 - 4*x**2 - 9*x + 21
    >>> prs[2]
    -5/9*x**4 + 1/9*x**2 - 1/3
    >>> prs[3]
    -117/25*x**2 - 9*x + 441/25
    >>> prs[4]
    233150/19773*x - 102500/6591
    >>> prs[5]
    -1288744821/543589225

    """
    prs = [f, g]
    h = dup_rem(f, g, K)

    while h:
        prs.append(h)
        f, g = g, h
        h = dup_rem(f, g, K)

    return prs


def dmp_euclidean_prs(f, g, u, K):
    """
    Euclidean polynomial remainder sequence (PRS) in `K[X]`.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    """
    if not u:
        return dup_euclidean_prs(f, g, K)
    else:
        raise MultivariatePolynomialError(f, g)


def dup_primitive_prs(f, g, K):
    """
    Primitive polynomial remainder sequence (PRS) in `K[x]`.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> f = x**8 + x**6 - 3*x**4 - 3*x**3 + 8*x**2 + 2*x - 5
    >>> g = 3*x**6 + 5*x**4 - 4*x**2 - 9*x + 21

    >>> prs = R.dup_primitive_prs(f, g)

    >>> prs[0]
    x**8 + x**6 - 3*x**4 - 3*x**3 + 8*x**2 + 2*x - 5
    >>> prs[1]
    3*x**6 + 5*x**4 - 4*x**2 - 9*x + 21
    >>> prs[2]
    -5*x**4 + x**2 - 3
    >>> prs[3]
    13*x**2 + 25*x - 49
    >>> prs[4]
    4663*x - 6150
    >>> prs[5]
    1

    """
    prs = [f, g]
    _, h = dup_primitive(dup_prem(f, g, K), K)

    while h:
        prs.append(h)
        f, g = g, h
        _, h = dup_primitive(dup_prem(f, g, K), K)

    return prs


def dmp_primitive_prs(f, g, u, K):
    """
    Primitive polynomial remainder sequence (PRS) in `K[X]`.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    """
    if not u:
        return dup_primitive_prs(f, g, K)
    else:
        raise MultivariatePolynomialError(f, g)


def dup_inner_subresultants(f, g, K):
    """
    Subresultant PRS algorithm in `K[x]`.

    Computes the subresultant polynomial remainder sequence (PRS)
    and the non-zero scalar subresultants of `f` and `g`.
    By [1] Thm. 3, these are the constants '-c' (- to optimize
    computation of sign).
    The first subdeterminant is set to 1 by convention to match
    the polynomial and the scalar subdeterminants.
    If 'deg(f) < deg(g)', the subresultants of '(g,f)' are computed.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_inner_subresultants(x**2 + 1, x**2 - 1)
    ([x**2 + 1, x**2 - 1, -2], [1, 1, 4])

    References
    ==========

    .. [1] W.S. Brown, The Subresultant PRS Algorithm.
           ACM Transaction of Mathematical Software 4 (1978) 237-249

    """
    n = dup_degree(f)
    m = dup_degree(g)

    if n < m:
        f, g = g, f
        n, m = m, n

    if not f:
        return [], []

    if not g:
        return [f], [K.one]

    R = [f, g]
    d = n - m

    b = (-K.one)**(d + 1)

    h = dup_prem(f, g, K)
    h = dup_mul_ground(h, b, K)

    lc = dup_LC(g, K)
    c = lc**d

    # Conventional first scalar subdeterminant is 1
    S = [K.one, c]
    c = -c

    while h:
        k = dup_degree(h)
        R.append(h)

        f, g, m, d = g, h, k, m - k

        b = -lc * c**d

        h = dup_prem(f, g, K)
        h = dup_quo_ground(h, b, K)

        lc = dup_LC(g, K)

        if d > 1:        # abnormal case
            q = c**(d - 1)
            c = K.quo((-lc)**d, q)
        else:
            c = -lc

        S.append(-c)

    return R, S


def dup_subresultants(f, g, K):
    """
    Computes subresultant PRS of two polynomials in `K[x]`.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_subresultants(x**2 + 1, x**2 - 1)
    [x**2 + 1, x**2 - 1, -2]

    """
    return dup_inner_subresultants(f, g, K)[0]


def dup_prs_resultant(f, g, K):
    """
    Resultant algorithm in `K[x]` using subresultant PRS.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_prs_resultant(x**2 + 1, x**2 - 1)
    (4, [x**2 + 1, x**2 - 1, -2])

    """
    if not f or not g:
        return (K.zero, [])

    R, S = dup_inner_subresultants(f, g, K)

    if dup_degree(R[-1]) > 0:
        return (K.zero, R)

    return S[-1], R


def dup_resultant(f, g, K, includePRS=False):
    """
    Computes resultant of two polynomials in `K[x]`.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_resultant(x**2 + 1, x**2 - 1)
    4

    """
    if includePRS:
        return dup_prs_resultant(f, g, K)
    return dup_prs_resultant(f, g, K)[0]


def dmp_inner_subresultants(f, g, u, K):
    """
    Subresultant PRS algorithm in `K[X]`.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> f = 3*x**2*y - y**3 - 4
    >>> g = x**2 + x*y**3 - 9

    >>> a = 3*x*y**4 + y**3 - 27*y + 4
    >>> b = -3*y**10 - 12*y**7 + y**6 - 54*y**4 + 8*y**3 + 729*y**2 - 216*y + 16

    >>> prs = [f, g, a, b]
    >>> sres = [[1], [1], [3, 0, 0, 0, 0], [-3, 0, 0, -12, 1, 0, -54, 8, 729, -216, 16]]

    >>> R.dmp_inner_subresultants(f, g) == (prs, sres)
    True

    """
    if not u:
        return dup_inner_subresultants(f, g, K)

    n = dmp_degree(f, u)
    m = dmp_degree(g, u)

    if n < m:
        f, g = g, f
        n, m = m, n

    if dmp_zero_p(f, u):
        return [], []

    v = u - 1
    if dmp_zero_p(g, u):
        return [f], [dmp_ground(K.one, v)]

    R = [f, g]
    d = n - m

    b = dmp_pow(dmp_ground(-K.one, v), d + 1, v, K)

    h = dmp_prem(f, g, u, K)
    h = dmp_mul_term(h, b, 0, u, K)

    lc = dmp_LC(g, K)
    c = dmp_pow(lc, d, v, K)

    S = [dmp_ground(K.one, v), c]
    c = dmp_neg(c, v, K)

    while not dmp_zero_p(h, u):
        k = dmp_degree(h, u)
        R.append(h)

        f, g, m, d = g, h, k, m - k

        b = dmp_mul(dmp_neg(lc, v, K),
                    dmp_pow(c, d, v, K), v, K)

        h = dmp_prem(f, g, u, K)
        h = [ dmp_quo(ch, b, v, K) for ch in h ]

        lc = dmp_LC(g, K)

        if d > 1:
            p = dmp_pow(dmp_neg(lc, v, K), d, v, K)
            q = dmp_pow(c, d - 1, v, K)
            c = dmp_quo(p, q, v, K)
        else:
            c = dmp_neg(lc, v, K)

        S.append(dmp_neg(c, v, K))

    return R, S


def dmp_subresultants(f, g, u, K):
    """
    Computes subresultant PRS of two polynomials in `K[X]`.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> f = 3*x**2*y - y**3 - 4
    >>> g = x**2 + x*y**3 - 9

    >>> a = 3*x*y**4 + y**3 - 27*y + 4
    >>> b = -3*y**10 - 12*y**7 + y**6 - 54*y**4 + 8*y**3 + 729*y**2 - 216*y + 16

    >>> R.dmp_subresultants(f, g) == [f, g, a, b]
    True

    """
    return dmp_inner_subresultants(f, g, u, K)[0]


def dmp_prs_resultant(f, g, u, K):
    """
    Resultant algorithm in `K[X]` using subresultant PRS.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> f = 3*x**2*y - y**3 - 4
    >>> g = x**2 + x*y**3 - 9

    >>> a = 3*x*y**4 + y**3 - 27*y + 4
    >>> b = -3*y**10 - 12*y**7 + y**6 - 54*y**4 + 8*y**3 + 729*y**2 - 216*y + 16

    >>> res, prs = R.dmp_prs_resultant(f, g)

    >>> res == b             # resultant has n-1 variables
    False
    >>> res == b.drop(x)
    True
    >>> prs == [f, g, a, b]
    True

    """
    if not u:
        return dup_prs_resultant(f, g, K)

    if dmp_zero_p(f, u) or dmp_zero_p(g, u):
        return (dmp_zero(u - 1), [])

    R, S = dmp_inner_subresultants(f, g, u, K)

    if dmp_degree(R[-1], u) > 0:
        return (dmp_zero(u - 1), R)

    return S[-1], R


def dmp_zz_modular_resultant(f, g, p, u, K):
    """
    Compute resultant of `f` and `g` modulo a prime `p`.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> f = x + y + 2
    >>> g = 2*x*y + x + 3

    >>> R.dmp_zz_modular_resultant(f, g, 5)
    -2*y**2 + 1

    """
    if not u:
        return gf_int(dup_prs_resultant(f, g, K)[0] % p, p)

    v = u - 1

    n = dmp_degree(f, u)
    m = dmp_degree(g, u)

    N = dmp_degree_in(f, 1, u)
    M = dmp_degree_in(g, 1, u)

    B = n*M + m*N

    D, a = [K.one], -K.one
    r = dmp_zero(v)

    while dup_degree(D) <= B:
        while True:
            a += K.one

            if a == p:
                raise HomomorphismFailed('no luck')

            F = dmp_eval_in(f, gf_int(a, p), 1, u, K)

            if dmp_degree(F, v) == n:
                G = dmp_eval_in(g, gf_int(a, p), 1, u, K)

                if dmp_degree(G, v) == m:
                    break

        R = dmp_zz_modular_resultant(F, G, p, v, K)
        e = dmp_eval(r, a, v, K)

        if not v:
            R = dup_strip([R])
            e = dup_strip([e])
        else:
            R = [R]
            e = [e]

        d = K.invert(dup_eval(D, a, K), p)
        d = dup_mul_ground(D, d, K)
        d = dmp_raise(d, v, 0, K)

        c = dmp_mul(d, dmp_sub(R, e, v, K), v, K)
        r = dmp_add(r, c, v, K)

        r = dmp_ground_trunc(r, p, v, K)

        D = dup_mul(D, [K.one, -a], K)
        D = dup_trunc(D, p, K)

    return r


def _collins_crt(r, R, P, p, K):
    """Wrapper of CRT for Collins's resultant algorithm. """
    return gf_int(gf_crt([r, R], [P, p], K), P*p)


def dmp_zz_collins_resultant(f, g, u, K):
    """
    Collins's modular resultant algorithm in `Z[X]`.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> f = x + y + 2
    >>> g = 2*x*y + x + 3

    >>> R.dmp_zz_collins_resultant(f, g)
    -2*y**2 - 5*y + 1

    """

    n = dmp_degree(f, u)
    m = dmp_degree(g, u)

    if n < 0 or m < 0:
        return dmp_zero(u - 1)

    A = dmp_max_norm(f, u, K)
    B = dmp_max_norm(g, u, K)

    a = dmp_ground_LC(f, u, K)
    b = dmp_ground_LC(g, u, K)

    v = u - 1

    B = K(2)*K.factorial(K(n + m))*A**m*B**n
    r, p, P = dmp_zero(v), K.one, K.one

    while P <= B:
        p = K(nextprime(p))

        while not (a % p) or not (b % p):
            p = K(nextprime(p))

        F = dmp_ground_trunc(f, p, u, K)
        G = dmp_ground_trunc(g, p, u, K)

        try:
            R = dmp_zz_modular_resultant(F, G, p, u, K)
        except HomomorphismFailed:
            continue

        if K.is_one(P):
            r = R
        else:
            r = dmp_apply_pairs(r, R, _collins_crt, (P, p, K), v, K)

        P *= p

    return r


def dmp_qq_collins_resultant(f, g, u, K0):
    """
    Collins's modular resultant algorithm in `Q[X]`.

    Examples
    ========

    >>> from sympy.polys import ring, QQ
    >>> R, x,y = ring("x,y", QQ)

    >>> f = QQ(1,2)*x + y + QQ(2,3)
    >>> g = 2*x*y + x + 3

    >>> R.dmp_qq_collins_resultant(f, g)
    -2*y**2 - 7/3*y + 5/6

    """
    n = dmp_degree(f, u)
    m = dmp_degree(g, u)

    if n < 0 or m < 0:
        return dmp_zero(u - 1)

    K1 = K0.get_ring()

    cf, f = dmp_clear_denoms(f, u, K0, K1)
    cg, g = dmp_clear_denoms(g, u, K0, K1)

    f = dmp_convert(f, u, K0, K1)
    g = dmp_convert(g, u, K0, K1)

    r = dmp_zz_collins_resultant(f, g, u, K1)
    r = dmp_convert(r, u - 1, K1, K0)

    c = K0.convert(cf**m * cg**n, K1)

    return dmp_quo_ground(r, c, u - 1, K0)


def dmp_resultant(f, g, u, K, includePRS=False):
    """
    Computes resultant of two polynomials in `K[X]`.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> f = 3*x**2*y - y**3 - 4
    >>> g = x**2 + x*y**3 - 9

    >>> R.dmp_resultant(f, g)
    -3*y**10 - 12*y**7 + y**6 - 54*y**4 + 8*y**3 + 729*y**2 - 216*y + 16

    """
    if not u:
        return dup_resultant(f, g, K, includePRS=includePRS)

    if includePRS:
        return dmp_prs_resultant(f, g, u, K)

    if K.is_Field:
        if K.is_QQ and query('USE_COLLINS_RESULTANT'):
            return dmp_qq_collins_resultant(f, g, u, K)
    else:
        if K.is_ZZ and query('USE_COLLINS_RESULTANT'):
            return dmp_zz_collins_resultant(f, g, u, K)

    return dmp_prs_resultant(f, g, u, K)[0]


def dup_discriminant(f, K):
    """
    Computes discriminant of a polynomial in `K[x]`.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_discriminant(x**2 + 2*x + 3)
    -8

    """
    d = dup_degree(f)

    if d <= 0:
        return K.zero
    else:
        s = (-1)**((d*(d - 1)) // 2)
        c = dup_LC(f, K)

        r = dup_resultant(f, dup_diff(f, 1, K), K)

        return K.quo(r, c*K(s))


def dmp_discriminant(f, u, K):
    """
    Computes discriminant of a polynomial in `K[X]`.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y,z,t = ring("x,y,z,t", ZZ)

    >>> R.dmp_discriminant(x**2*y + x*z + t)
    -4*y*t + z**2

    """
    if not u:
        return dup_discriminant(f, K)

    d, v = dmp_degree(f, u), u - 1

    if d <= 0:
        return dmp_zero(v)
    else:
        s = (-1)**((d*(d - 1)) // 2)
        c = dmp_LC(f, K)

        r = dmp_resultant(f, dmp_diff(f, 1, u, K), u, K)
        c = dmp_mul_ground(c, K(s), v, K)

        return dmp_quo(r, c, v, K)


def _dup_rr_trivial_gcd(f, g, K):
    """Handle trivial cases in GCD algorithm over a ring. """
    if not (f or g):
        return [], [], []
    elif not f:
        if K.is_nonnegative(dup_LC(g, K)):
            return g, [], [K.one]
        else:
            return dup_neg(g, K), [], [-K.one]
    elif not g:
        if K.is_nonnegative(dup_LC(f, K)):
            return f, [K.one], []
        else:
            return dup_neg(f, K), [-K.one], []

    return None


def _dup_ff_trivial_gcd(f, g, K):
    """Handle trivial cases in GCD algorithm over a field. """
    if not (f or g):
        return [], [], []
    elif not f:
        return dup_monic(g, K), [], [dup_LC(g, K)]
    elif not g:
        return dup_monic(f, K), [dup_LC(f, K)], []
    else:
        return None


def _dmp_rr_trivial_gcd(f, g, u, K):
    """Handle trivial cases in GCD algorithm over a ring. """
    zero_f = dmp_zero_p(f, u)
    zero_g = dmp_zero_p(g, u)
    if_contain_one = dmp_one_p(f, u, K) or dmp_one_p(g, u, K)

    if zero_f and zero_g:
        return tuple(dmp_zeros(3, u, K))
    elif zero_f:
        if K.is_nonnegative(dmp_ground_LC(g, u, K)):
            return g, dmp_zero(u), dmp_one(u, K)
        else:
            return dmp_neg(g, u, K), dmp_zero(u), dmp_ground(-K.one, u)
    elif zero_g:
        if K.is_nonnegative(dmp_ground_LC(f, u, K)):
            return f, dmp_one(u, K), dmp_zero(u)
        else:
            return dmp_neg(f, u, K), dmp_ground(-K.one, u), dmp_zero(u)
    elif if_contain_one:
        return dmp_one(u, K), f, g
    elif query('USE_SIMPLIFY_GCD'):
        return _dmp_simplify_gcd(f, g, u, K)
    else:
        return None


def _dmp_ff_trivial_gcd(f, g, u, K):
    """Handle trivial cases in GCD algorithm over a field. """
    zero_f = dmp_zero_p(f, u)
    zero_g = dmp_zero_p(g, u)

    if zero_f and zero_g:
        return tuple(dmp_zeros(3, u, K))
    elif zero_f:
        return (dmp_ground_monic(g, u, K),
                dmp_zero(u),
                dmp_ground(dmp_ground_LC(g, u, K), u))
    elif zero_g:
        return (dmp_ground_monic(f, u, K),
                dmp_ground(dmp_ground_LC(f, u, K), u),
                dmp_zero(u))
    elif query('USE_SIMPLIFY_GCD'):
        return _dmp_simplify_gcd(f, g, u, K)
    else:
        return None


def _dmp_simplify_gcd(f, g, u, K):
    """Try to eliminate `x_0` from GCD computation in `K[X]`. """
    df = dmp_degree(f, u)
    dg = dmp_degree(g, u)

    if df > 0 and dg > 0:
        return None

    if not (df or dg):
        F = dmp_LC(f, K)
        G = dmp_LC(g, K)
    else:
        if not df:
            F = dmp_LC(f, K)
            G = dmp_content(g, u, K)
        else:
            F = dmp_content(f, u, K)
            G = dmp_LC(g, K)

    v = u - 1
    h = dmp_gcd(F, G, v, K)

    cff = [ dmp_quo(cf, h, v, K) for cf in f ]
    cfg = [ dmp_quo(cg, h, v, K) for cg in g ]

    return [h], cff, cfg


def dup_rr_prs_gcd(f, g, K):
    """
    Computes polynomial GCD using subresultants over a ring.

    Returns ``(h, cff, cfg)`` such that ``a = gcd(f, g)``, ``cff = quo(f, h)``,
    and ``cfg = quo(g, h)``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_rr_prs_gcd(x**2 - 1, x**2 - 3*x + 2)
    (x - 1, x + 1, x - 2)

    """
    result = _dup_rr_trivial_gcd(f, g, K)

    if result is not None:
        return result

    fc, F = dup_primitive(f, K)
    gc, G = dup_primitive(g, K)

    c = K.gcd(fc, gc)

    h = dup_subresultants(F, G, K)[-1]
    _, h = dup_primitive(h, K)

    c *= K.canonical_unit(dup_LC(h, K))

    h = dup_mul_ground(h, c, K)

    cff = dup_quo(f, h, K)
    cfg = dup_quo(g, h, K)

    return h, cff, cfg


def dup_ff_prs_gcd(f, g, K):
    """
    Computes polynomial GCD using subresultants over a field.

    Returns ``(h, cff, cfg)`` such that ``a = gcd(f, g)``, ``cff = quo(f, h)``,
    and ``cfg = quo(g, h)``.

    Examples
    ========

    >>> from sympy.polys import ring, QQ
    >>> R, x = ring("x", QQ)

    >>> R.dup_ff_prs_gcd(x**2 - 1, x**2 - 3*x + 2)
    (x - 1, x + 1, x - 2)

    """
    result = _dup_ff_trivial_gcd(f, g, K)

    if result is not None:
        return result

    h = dup_subresultants(f, g, K)[-1]
    h = dup_monic(h, K)

    cff = dup_quo(f, h, K)
    cfg = dup_quo(g, h, K)

    return h, cff, cfg


def dmp_rr_prs_gcd(f, g, u, K):
    """
    Computes polynomial GCD using subresultants over a ring.

    Returns ``(h, cff, cfg)`` such that ``a = gcd(f, g)``, ``cff = quo(f, h)``,
    and ``cfg = quo(g, h)``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y, = ring("x,y", ZZ)

    >>> f = x**2 + 2*x*y + y**2
    >>> g = x**2 + x*y

    >>> R.dmp_rr_prs_gcd(f, g)
    (x + y, x + y, x)

    """
    if not u:
        return dup_rr_prs_gcd(f, g, K)

    result = _dmp_rr_trivial_gcd(f, g, u, K)

    if result is not None:
        return result

    fc, F = dmp_primitive(f, u, K)
    gc, G = dmp_primitive(g, u, K)

    h = dmp_subresultants(F, G, u, K)[-1]
    c, _, _ = dmp_rr_prs_gcd(fc, gc, u - 1, K)

    if K.is_negative(dmp_ground_LC(h, u, K)):
        h = dmp_neg(h, u, K)

    _, h = dmp_primitive(h, u, K)
    h = dmp_mul_term(h, c, 0, u, K)

    cff = dmp_quo(f, h, u, K)
    cfg = dmp_quo(g, h, u, K)

    return h, cff, cfg


def dmp_ff_prs_gcd(f, g, u, K):
    """
    Computes polynomial GCD using subresultants over a field.

    Returns ``(h, cff, cfg)`` such that ``a = gcd(f, g)``, ``cff = quo(f, h)``,
    and ``cfg = quo(g, h)``.

    Examples
    ========

    >>> from sympy.polys import ring, QQ
    >>> R, x,y, = ring("x,y", QQ)

    >>> f = QQ(1,2)*x**2 + x*y + QQ(1,2)*y**2
    >>> g = x**2 + x*y

    >>> R.dmp_ff_prs_gcd(f, g)
    (x + y, 1/2*x + 1/2*y, x)

    """
    if not u:
        return dup_ff_prs_gcd(f, g, K)

    result = _dmp_ff_trivial_gcd(f, g, u, K)

    if result is not None:
        return result

    fc, F = dmp_primitive(f, u, K)
    gc, G = dmp_primitive(g, u, K)

    h = dmp_subresultants(F, G, u, K)[-1]
    c, _, _ = dmp_ff_prs_gcd(fc, gc, u - 1, K)

    _, h = dmp_primitive(h, u, K)
    h = dmp_mul_term(h, c, 0, u, K)
    h = dmp_ground_monic(h, u, K)

    cff = dmp_quo(f, h, u, K)
    cfg = dmp_quo(g, h, u, K)

    return h, cff, cfg

HEU_GCD_MAX = 6


def _dup_zz_gcd_interpolate(h, x, K):
    """Interpolate polynomial GCD from integer GCD. """
    f = []

    while h:
        g = h % x

        if g > x // 2:
            g -= x

        f.insert(0, g)
        h = (h - g) // x

    return f


def dup_zz_heu_gcd(f, g, K):
    """
    Heuristic polynomial GCD in `Z[x]`.

    Given univariate polynomials `f` and `g` in `Z[x]`, returns
    their GCD and cofactors, i.e. polynomials ``h``, ``cff`` and ``cfg``
    such that::

          h = gcd(f, g), cff = quo(f, h) and cfg = quo(g, h)

    The algorithm is purely heuristic which means it may fail to compute
    the GCD. This will be signaled by raising an exception. In this case
    you will need to switch to another GCD method.

    The algorithm computes the polynomial GCD by evaluating polynomials
    f and g at certain points and computing (fast) integer GCD of those
    evaluations. The polynomial GCD is recovered from the integer image
    by interpolation.  The final step is to verify if the result is the
    correct GCD. This gives cofactors as a side effect.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_zz_heu_gcd(x**2 - 1, x**2 - 3*x + 2)
    (x - 1, x + 1, x - 2)

    References
    ==========

    .. [1] [Liao95]_

    """
    result = _dup_rr_trivial_gcd(f, g, K)

    if result is not None:
        return result

    df = dup_degree(f)
    dg = dup_degree(g)

    gcd, f, g = dup_extract(f, g, K)

    if df == 0 or dg == 0:
        return [gcd], f, g

    f_norm = dup_max_norm(f, K)
    g_norm = dup_max_norm(g, K)

    B = K(2*min(f_norm, g_norm) + 29)

    x = max(min(B, 99*K.sqrt(B)),
            2*min(f_norm // abs(dup_LC(f, K)),
                  g_norm // abs(dup_LC(g, K))) + 2)

    for i in range(0, HEU_GCD_MAX):
        ff = dup_eval(f, x, K)
        gg = dup_eval(g, x, K)

        if ff and gg:
            h = K.gcd(ff, gg)

            cff = ff // h
            cfg = gg // h

            h = _dup_zz_gcd_interpolate(h, x, K)
            h = dup_primitive(h, K)[1]

            cff_, r = dup_div(f, h, K)

            if not r:
                cfg_, r = dup_div(g, h, K)

                if not r:
                    h = dup_mul_ground(h, gcd, K)
                    return h, cff_, cfg_

            cff = _dup_zz_gcd_interpolate(cff, x, K)

            h, r = dup_div(f, cff, K)

            if not r:
                cfg_, r = dup_div(g, h, K)

                if not r:
                    h = dup_mul_ground(h, gcd, K)
                    return h, cff, cfg_

            cfg = _dup_zz_gcd_interpolate(cfg, x, K)

            h, r = dup_div(g, cfg, K)

            if not r:
                cff_, r = dup_div(f, h, K)

                if not r:
                    h = dup_mul_ground(h, gcd, K)
                    return h, cff_, cfg

        x = 73794*x * K.sqrt(K.sqrt(x)) // 27011

    raise HeuristicGCDFailed('no luck')


def _dmp_zz_gcd_interpolate(h, x, v, K):
    """Interpolate polynomial GCD from integer GCD. """
    f = []

    while not dmp_zero_p(h, v):
        g = dmp_ground_trunc(h, x, v, K)
        f.insert(0, g)

        h = dmp_sub(h, g, v, K)
        h = dmp_quo_ground(h, x, v, K)

    if K.is_negative(dmp_ground_LC(f, v + 1, K)):
        return dmp_neg(f, v + 1, K)
    else:
        return f


def dmp_zz_heu_gcd(f, g, u, K):
    """
    Heuristic polynomial GCD in `Z[X]`.

    Given univariate polynomials `f` and `g` in `Z[X]`, returns
    their GCD and cofactors, i.e. polynomials ``h``, ``cff`` and ``cfg``
    such that::

          h = gcd(f, g), cff = quo(f, h) and cfg = quo(g, h)

    The algorithm is purely heuristic which means it may fail to compute
    the GCD. This will be signaled by raising an exception. In this case
    you will need to switch to another GCD method.

    The algorithm computes the polynomial GCD by evaluating polynomials
    f and g at certain points and computing (fast) integer GCD of those
    evaluations. The polynomial GCD is recovered from the integer image
    by interpolation. The evaluation process reduces f and g variable by
    variable into a large integer.  The final step is to verify if the
    interpolated polynomial is the correct GCD. This gives cofactors of
    the input polynomials as a side effect.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y, = ring("x,y", ZZ)

    >>> f = x**2 + 2*x*y + y**2
    >>> g = x**2 + x*y

    >>> R.dmp_zz_heu_gcd(f, g)
    (x + y, x + y, x)

    References
    ==========

    .. [1] [Liao95]_

    """
    if not u:
        return dup_zz_heu_gcd(f, g, K)

    result = _dmp_rr_trivial_gcd(f, g, u, K)

    if result is not None:
        return result

    gcd, f, g = dmp_ground_extract(f, g, u, K)

    f_norm = dmp_max_norm(f, u, K)
    g_norm = dmp_max_norm(g, u, K)

    B = K(2*min(f_norm, g_norm) + 29)

    x = max(min(B, 99*K.sqrt(B)),
            2*min(f_norm // abs(dmp_ground_LC(f, u, K)),
                  g_norm // abs(dmp_ground_LC(g, u, K))) + 2)

    for i in range(0, HEU_GCD_MAX):
        ff = dmp_eval(f, x, u, K)
        gg = dmp_eval(g, x, u, K)

        v = u - 1

        if not (dmp_zero_p(ff, v) or dmp_zero_p(gg, v)):
            h, cff, cfg = dmp_zz_heu_gcd(ff, gg, v, K)

            h = _dmp_zz_gcd_interpolate(h, x, v, K)
            h = dmp_ground_primitive(h, u, K)[1]

            cff_, r = dmp_div(f, h, u, K)

            if dmp_zero_p(r, u):
                cfg_, r = dmp_div(g, h, u, K)

                if dmp_zero_p(r, u):
                    h = dmp_mul_ground(h, gcd, u, K)
                    return h, cff_, cfg_

            cff = _dmp_zz_gcd_interpolate(cff, x, v, K)

            h, r = dmp_div(f, cff, u, K)

            if dmp_zero_p(r, u):
                cfg_, r = dmp_div(g, h, u, K)

                if dmp_zero_p(r, u):
                    h = dmp_mul_ground(h, gcd, u, K)
                    return h, cff, cfg_

            cfg = _dmp_zz_gcd_interpolate(cfg, x, v, K)

            h, r = dmp_div(g, cfg, u, K)

            if dmp_zero_p(r, u):
                cff_, r = dmp_div(f, h, u, K)

                if dmp_zero_p(r, u):
                    h = dmp_mul_ground(h, gcd, u, K)
                    return h, cff_, cfg

        x = 73794*x * K.sqrt(K.sqrt(x)) // 27011

    raise HeuristicGCDFailed('no luck')


def dup_qq_heu_gcd(f, g, K0):
    """
    Heuristic polynomial GCD in `Q[x]`.

    Returns ``(h, cff, cfg)`` such that ``a = gcd(f, g)``,
    ``cff = quo(f, h)``, and ``cfg = quo(g, h)``.

    Examples
    ========

    >>> from sympy.polys import ring, QQ
    >>> R, x = ring("x", QQ)

    >>> f = QQ(1,2)*x**2 + QQ(7,4)*x + QQ(3,2)
    >>> g = QQ(1,2)*x**2 + x

    >>> R.dup_qq_heu_gcd(f, g)
    (x + 2, 1/2*x + 3/4, 1/2*x)

    """
    result = _dup_ff_trivial_gcd(f, g, K0)

    if result is not None:
        return result

    K1 = K0.get_ring()

    cf, f = dup_clear_denoms(f, K0, K1)
    cg, g = dup_clear_denoms(g, K0, K1)

    f = dup_convert(f, K0, K1)
    g = dup_convert(g, K0, K1)

    h, cff, cfg = dup_zz_heu_gcd(f, g, K1)

    h = dup_convert(h, K1, K0)

    c = dup_LC(h, K0)
    h = dup_monic(h, K0)

    cff = dup_convert(cff, K1, K0)
    cfg = dup_convert(cfg, K1, K0)

    cff = dup_mul_ground(cff, K0.quo(c, cf), K0)
    cfg = dup_mul_ground(cfg, K0.quo(c, cg), K0)

    return h, cff, cfg


def dmp_qq_heu_gcd(f, g, u, K0):
    """
    Heuristic polynomial GCD in `Q[X]`.

    Returns ``(h, cff, cfg)`` such that ``a = gcd(f, g)``,
    ``cff = quo(f, h)``, and ``cfg = quo(g, h)``.

    Examples
    ========

    >>> from sympy.polys import ring, QQ
    >>> R, x,y, = ring("x,y", QQ)

    >>> f = QQ(1,4)*x**2 + x*y + y**2
    >>> g = QQ(1,2)*x**2 + x*y

    >>> R.dmp_qq_heu_gcd(f, g)
    (x + 2*y, 1/4*x + 1/2*y, 1/2*x)

    """
    result = _dmp_ff_trivial_gcd(f, g, u, K0)

    if result is not None:
        return result

    K1 = K0.get_ring()

    cf, f = dmp_clear_denoms(f, u, K0, K1)
    cg, g = dmp_clear_denoms(g, u, K0, K1)

    f = dmp_convert(f, u, K0, K1)
    g = dmp_convert(g, u, K0, K1)

    h, cff, cfg = dmp_zz_heu_gcd(f, g, u, K1)

    h = dmp_convert(h, u, K1, K0)

    c = dmp_ground_LC(h, u, K0)
    h = dmp_ground_monic(h, u, K0)

    cff = dmp_convert(cff, u, K1, K0)
    cfg = dmp_convert(cfg, u, K1, K0)

    cff = dmp_mul_ground(cff, K0.quo(c, cf), u, K0)
    cfg = dmp_mul_ground(cfg, K0.quo(c, cg), u, K0)

    return h, cff, cfg


def dup_inner_gcd(f, g, K):
    """
    Computes polynomial GCD and cofactors of `f` and `g` in `K[x]`.

    Returns ``(h, cff, cfg)`` such that ``a = gcd(f, g)``,
    ``cff = quo(f, h)``, and ``cfg = quo(g, h)``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_inner_gcd(x**2 - 1, x**2 - 3*x + 2)
    (x - 1, x + 1, x - 2)

    """
    if not K.is_Exact:
        try:
            exact = K.get_exact()
        except DomainError:
            return [K.one], f, g

        f = dup_convert(f, K, exact)
        g = dup_convert(g, K, exact)

        h, cff, cfg = dup_inner_gcd(f, g, exact)

        h = dup_convert(h, exact, K)
        cff = dup_convert(cff, exact, K)
        cfg = dup_convert(cfg, exact, K)

        return h, cff, cfg
    elif K.is_Field:
        if K.is_QQ and query('USE_HEU_GCD'):
            try:
                return dup_qq_heu_gcd(f, g, K)
            except HeuristicGCDFailed:
                pass

        return dup_ff_prs_gcd(f, g, K)
    else:
        if K.is_ZZ and query('USE_HEU_GCD'):
            try:
                return dup_zz_heu_gcd(f, g, K)
            except HeuristicGCDFailed:
                pass

        return dup_rr_prs_gcd(f, g, K)


def _dmp_inner_gcd(f, g, u, K):
    """Helper function for `dmp_inner_gcd()`. """
    if not K.is_Exact:
        try:
            exact = K.get_exact()
        except DomainError:
            return dmp_one(u, K), f, g

        f = dmp_convert(f, u, K, exact)
        g = dmp_convert(g, u, K, exact)

        h, cff, cfg = _dmp_inner_gcd(f, g, u, exact)

        h = dmp_convert(h, u, exact, K)
        cff = dmp_convert(cff, u, exact, K)
        cfg = dmp_convert(cfg, u, exact, K)

        return h, cff, cfg
    elif K.is_Field:
        if K.is_QQ and query('USE_HEU_GCD'):
            try:
                return dmp_qq_heu_gcd(f, g, u, K)
            except HeuristicGCDFailed:
                pass

        return dmp_ff_prs_gcd(f, g, u, K)
    else:
        if K.is_ZZ and query('USE_HEU_GCD'):
            try:
                return dmp_zz_heu_gcd(f, g, u, K)
            except HeuristicGCDFailed:
                pass

        return dmp_rr_prs_gcd(f, g, u, K)


def dmp_inner_gcd(f, g, u, K):
    """
    Computes polynomial GCD and cofactors of `f` and `g` in `K[X]`.

    Returns ``(h, cff, cfg)`` such that ``a = gcd(f, g)``,
    ``cff = quo(f, h)``, and ``cfg = quo(g, h)``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y, = ring("x,y", ZZ)

    >>> f = x**2 + 2*x*y + y**2
    >>> g = x**2 + x*y

    >>> R.dmp_inner_gcd(f, g)
    (x + y, x + y, x)

    """
    if not u:
        return dup_inner_gcd(f, g, K)

    J, (f, g) = dmp_multi_deflate((f, g), u, K)
    h, cff, cfg = _dmp_inner_gcd(f, g, u, K)

    return (dmp_inflate(h, J, u, K),
            dmp_inflate(cff, J, u, K),
            dmp_inflate(cfg, J, u, K))


def dup_gcd(f, g, K):
    """
    Computes polynomial GCD of `f` and `g` in `K[x]`.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_gcd(x**2 - 1, x**2 - 3*x + 2)
    x - 1

    """
    return dup_inner_gcd(f, g, K)[0]


def dmp_gcd(f, g, u, K):
    """
    Computes polynomial GCD of `f` and `g` in `K[X]`.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y, = ring("x,y", ZZ)

    >>> f = x**2 + 2*x*y + y**2
    >>> g = x**2 + x*y

    >>> R.dmp_gcd(f, g)
    x + y

    """
    return dmp_inner_gcd(f, g, u, K)[0]


def dup_rr_lcm(f, g, K):
    """
    Computes polynomial LCM over a ring in `K[x]`.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_rr_lcm(x**2 - 1, x**2 - 3*x + 2)
    x**3 - 2*x**2 - x + 2

    """
    fc, f = dup_primitive(f, K)
    gc, g = dup_primitive(g, K)

    c = K.lcm(fc, gc)

    h = dup_quo(dup_mul(f, g, K),
                dup_gcd(f, g, K), K)

    return dup_mul_ground(h, c, K)


def dup_ff_lcm(f, g, K):
    """
    Computes polynomial LCM over a field in `K[x]`.

    Examples
    ========

    >>> from sympy.polys import ring, QQ
    >>> R, x = ring("x", QQ)

    >>> f = QQ(1,2)*x**2 + QQ(7,4)*x + QQ(3,2)
    >>> g = QQ(1,2)*x**2 + x

    >>> R.dup_ff_lcm(f, g)
    x**3 + 7/2*x**2 + 3*x

    """
    h = dup_quo(dup_mul(f, g, K),
                dup_gcd(f, g, K), K)

    return dup_monic(h, K)


def dup_lcm(f, g, K):
    """
    Computes polynomial LCM of `f` and `g` in `K[x]`.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_lcm(x**2 - 1, x**2 - 3*x + 2)
    x**3 - 2*x**2 - x + 2

    """
    if K.is_Field:
        return dup_ff_lcm(f, g, K)
    else:
        return dup_rr_lcm(f, g, K)


def dmp_rr_lcm(f, g, u, K):
    """
    Computes polynomial LCM over a ring in `K[X]`.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y, = ring("x,y", ZZ)

    >>> f = x**2 + 2*x*y + y**2
    >>> g = x**2 + x*y

    >>> R.dmp_rr_lcm(f, g)
    x**3 + 2*x**2*y + x*y**2

    """
    fc, f = dmp_ground_primitive(f, u, K)
    gc, g = dmp_ground_primitive(g, u, K)

    c = K.lcm(fc, gc)

    h = dmp_quo(dmp_mul(f, g, u, K),
                dmp_gcd(f, g, u, K), u, K)

    return dmp_mul_ground(h, c, u, K)


def dmp_ff_lcm(f, g, u, K):
    """
    Computes polynomial LCM over a field in `K[X]`.

    Examples
    ========

    >>> from sympy.polys import ring, QQ
    >>> R, x,y, = ring("x,y", QQ)

    >>> f = QQ(1,4)*x**2 + x*y + y**2
    >>> g = QQ(1,2)*x**2 + x*y

    >>> R.dmp_ff_lcm(f, g)
    x**3 + 4*x**2*y + 4*x*y**2

    """
    h = dmp_quo(dmp_mul(f, g, u, K),
                dmp_gcd(f, g, u, K), u, K)

    return dmp_ground_monic(h, u, K)


def dmp_lcm(f, g, u, K):
    """
    Computes polynomial LCM of `f` and `g` in `K[X]`.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y, = ring("x,y", ZZ)

    >>> f = x**2 + 2*x*y + y**2
    >>> g = x**2 + x*y

    >>> R.dmp_lcm(f, g)
    x**3 + 2*x**2*y + x*y**2

    """
    if not u:
        return dup_lcm(f, g, K)

    if K.is_Field:
        return dmp_ff_lcm(f, g, u, K)
    else:
        return dmp_rr_lcm(f, g, u, K)


def dmp_content(f, u, K):
    """
    Returns GCD of multivariate coefficients.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y, = ring("x,y", ZZ)

    >>> R.dmp_content(2*x*y + 6*x + 4*y + 12)
    2*y + 6

    """
    cont, v = dmp_LC(f, K), u - 1

    if dmp_zero_p(f, u):
        return cont

    for c in f[1:]:
        cont = dmp_gcd(cont, c, v, K)

        if dmp_one_p(cont, v, K):
            break

    if K.is_negative(dmp_ground_LC(cont, v, K)):
        return dmp_neg(cont, v, K)
    else:
        return cont


def dmp_primitive(f, u, K):
    """
    Returns multivariate content and a primitive polynomial.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y, = ring("x,y", ZZ)

    >>> R.dmp_primitive(2*x*y + 6*x + 4*y + 12)
    (2*y + 6, x + 2)

    """
    cont, v = dmp_content(f, u, K), u - 1

    if dmp_zero_p(f, u) or dmp_one_p(cont, v, K):
        return cont, f
    else:
        return cont, [ dmp_quo(c, cont, v, K) for c in f ]


def dup_cancel(f, g, K, include=True):
    """
    Cancel common factors in a rational function `f/g`.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_cancel(2*x**2 - 2, x**2 - 2*x + 1)
    (2*x + 2, x - 1)

    """
    return dmp_cancel(f, g, 0, K, include=include)


def dmp_cancel(f, g, u, K, include=True):
    """
    Cancel common factors in a rational function `f/g`.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dmp_cancel(2*x**2 - 2, x**2 - 2*x + 1)
    (2*x + 2, x - 1)

    """
    K0 = None

    if K.is_Field and K.has_assoc_Ring:
        K0, K = K, K.get_ring()

        cq, f = dmp_clear_denoms(f, u, K0, K, convert=True)
        cp, g = dmp_clear_denoms(g, u, K0, K, convert=True)
    else:
        cp, cq = K.one, K.one

    _, p, q = dmp_inner_gcd(f, g, u, K)

    if K0 is not None:
        _, cp, cq = K.cofactors(cp, cq)

        p = dmp_convert(p, u, K, K0)
        q = dmp_convert(q, u, K, K0)

        K = K0

    p_neg = K.is_negative(dmp_ground_LC(p, u, K))
    q_neg = K.is_negative(dmp_ground_LC(q, u, K))

    if p_neg and q_neg:
        p, q = dmp_neg(p, u, K), dmp_neg(q, u, K)
    elif p_neg:
        cp, p = -cp, dmp_neg(p, u, K)
    elif q_neg:
        cp, q = -cp, dmp_neg(q, u, K)

    if not include:
        return cp, cq, p, q

    p = dmp_mul_ground(p, cp, u, K)
    q = dmp_mul_ground(q, cq, u, K)

    return p, q


def free_symbols_sparse(p):
    exponents = list(map(sum, zip(*p)))
    return {n for n, e in enumerate(exponents) if e}


def _coeffs_sparse(p1, syms):
    N = len(p1.ring.gens)
    nsyms = set(range(N)) - syms
    p2 = defaultdict(dict)
    r = range(N)
    for m1, c1 in p1.items():
        # m21 = tuple(m1[i] if i in syms else 0 for i in range(N))
        # m22 = tuple(m1[i] if i in nsyms else 0 for i in range(N))
        # p2[m21][m22] = c1
        sym_indices = set(compress(r, m1))
        m21 = [0] * N
        m22 = [0] * N
        for i in sym_indices & syms:
            m21[i] = m1[i]
        for i in sym_indices & nsyms:
            m22[i] = m1[i]
        p2[tuple(m21)][tuple(m22)] = c1
    return p2


def coeffs_sparse(p1, syms):
    p2 = _coeffs_sparse(p1, syms)
    return [p1.ring(pi) for pi in p2.values()]


def coeff_sparse(p1, sym, deg):
    p2 = _coeffs_sparse(p1, {sym})
    m = [0] * len(p1.ring.gens)
    m[sym] = deg
    m = tuple(m)
    return p1.ring(p2[m])


def gcd_sparse(p):
    R = p[0].ring
    K = R.domain

    # 1-term gcd?
    if any(len(pi) == 1 for pi in p):
        return gcd_terms(p, R, K)

    # Extract the monomial gcd
    p, d = gcd_monomial(p)

    # Eliminate nonshared symbols
    p, common = gcd_coeffs(p)

    # Use subresultant PRS
    g = p[0]
    for pi in p[1:]:
        g = gcd_prs_sparse(g, pi)

    if d is not None:
        g = g * d

    return g


def gcd_terms(p, R, K):
    #
    # If any of the polynomials only has a single term then the gcd does as
    # well and must divide all terms. We can bypass all of the more
    # complicated algorithms and find the gcd of the monomials and the gcd of
    # the coefficients.
    #
    terms = chain.from_iterable(pi.terms() for pi in p)
    monoms = set()
    coeffs = set()
    for m, c in terms:
        monoms.add(m)
        coeffs.add(c)
    monom_gcd = gcd_monom(monoms)
    coeff_gcd = gcd_ground(coeffs, K)
    term_gcd = R({monom_gcd: coeff_gcd})
    return term_gcd


def gcd_ground(c, K):
    #
    # gcd in the ground domain
    #
    c = list(c)
    gcd = K.gcd
    d = c[0]
    for ci in c[1:]:
        d = gcd(d, ci)
        if d == K.one:
            break
    return d


def gcd_monom(monoms):
    #
    # gcd of a list of monomials
    #
    monom_gcd = tuple(map(min, zip(*monoms)))
    return monom_gcd


def gcd_monomial(p):
    #
    # Extract any common monomial from the polynomials in p
    #
    R = p[0].ring
    monoms = chain(*p)
    monom_gcd = tuple(map(min, zip(*monoms)))
    if monom_gcd == R.zero_monom:
        return p, None
    else:
        d = R({monom_gcd: R.domain.one})
        p = [pi.exquo(d) for pi in p]
        return p, d


def gcd_coeffs(p):
    #
    # Simplify a list of polys whose gcd is wanted. Returns a possibly longer
    # list of simpler polys having the same gcd as the input. This is done by
    # eliminating symbols that can not be part of the gcd because they do not
    # appear in each item of the input. In the output list all items have
    # exactly the same symbols. The set of those symbols is also returned.
    #
    # Given d = gcd(p(x, y), q(y, z)) then d = r(y) does not depend on x or z.
    # Write p(x, y) as a poly in x with coeffs that are polys in y:
    #
    #   p(x, y) = p0(y) + p1(y)*x + p2(y)*x**2 + ...
    #   q(x, y) = q0(y) + q1(y)*z + q2(y)*z**2 + ...
    #
    # Now r(y) = gcd(p0, p1, ..., q0, q1, ...)
    #
    # We can eliminate all symbols that do not appear in all polys until every
    # poly has exactly the same symbols in it.
    #
    # This can very quickly pick out a term with a coefficient of 1 without
    # the need to do any poly arithmetic. If 1 or -1 is encountered the
    # routine returns [1] immediately.
    #
    # In computing gcd of multivariate polynomials this routine is an
    # important optimisation for the common case where the gcd is just 1. Even
    # if the gcd is not 1 we can still reduce the number of symbols
    # drastically in some cases before proceeding to e.g. subresultant PRS
    # where performance depends strongly on the number of symbols.
    #
    # Example:
    #
    #   >>> x = symbols('x:10')
    #   >>> K = QQ[x]
    #   >>> p1 = K.from_sympy(sum(x[:8]))
    #   >>> p2 = K.from_sympy(sum(x[2:]))
    #   >>> p1
    #   x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7
    #   >>> p2
    #   x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
    #   >>> gcd_coeffs([p1, p2])
    #   ([1], set())
    #
    all_coeffs = p

    while True:
        # Quick exits are most efficient if we start from the simplest polys
        p = sorted(set(all_coeffs), key=len)

        # Find the intersection of symbols for each poly:
        common = free_symbols_sparse(p[0])
        allsame = True
        for pi in p[1:]:

            # Quick exit
            if not common:
                R = p[0].ring
                K = R.domain
                gcd = gcd_terms(p, R, K)
                return [gcd], None

            syms = free_symbols_sparse(pi)
            if allsame and syms != common:
                allsame = False
            common &= syms

        # The loop is complete if they all have the same symbols.
        if allsame:
            return p, common

        # Extract coefficients as polys containing only the common symbols.
        all_coeffs = []
        for i, pi in enumerate(p):
            coeffs_i = coeffs_sparse(pi, free_symbols_sparse(pi) - common)
            all_coeffs.extend(coeffs_i)

            # Quick exit:
            if any(len(c) == 1 for c in coeffs_i):
                R = p[0].ring
                K = R.domain
                gcd = gcd_terms(all_coeffs + p[i+1:], R, K)
                return [gcd], None


def prem_sparse(f, g, x):
    df = f.degree(x)
    dg = g.degree(x)

    if dg < 0:
        raise ZeroDivisionError

    r, dr = f, df

    if df < dg:
        return r

    N = df - dg + 1

    lc_g = coeff_sparse(g, x, dg)

    xp = f.ring.gens[x]

    while True:
        lc_r = coeff_sparse(r, x, dr)
        j, N = dr - dg, N - 1

        R = r * lc_g
        G = g * lc_r * xp**j
        r = R - G

        _dr, dr = dr, r.degree(x)

        if dr < dg:
            break
        elif not (dr < _dr):
            raise ValueError

    c = lc_g ** N

    return r * c


def subresultants_sparse(f, g, x):
    n = f.degree(x)
    m = g.degree(x)

    if n < m:
        f, g = g, f
        n, m = m, n

    if f == 0:
        return 0, 0

    if g == 0:
        return f, 1

    R = [f, g]
    d = n - m

    b = (-1) ** (d + 1)

    h = prem_sparse(f, g, x)
    h = h * b

    lc = coeff_sparse(g, x, m)
    c = lc ** d

    S = [1, c]
    c = -c

    while h:
        k = h.degree(x)
        R.append(h)

        f, g, m, d = g, h, k, m-k

        b = -lc * c**d

        h = prem_sparse(f, g, x)
        h = h.exquo(b)

        lc = coeff_sparse(g, x, k)

        if d > 1:
            p = (-lc) ** d
            q = c ** (d-1)
            c = p.exquo(q)
        else:
            c = -lc

        S.append(-c)

    return R


def primitive_sparse(p, x):
    coeffs = coeffs_sparse(p, {x})
    content = gcd_sparse(coeffs)
    primitive = p.exquo(content)
    return content, primitive


def main_variable_sparse(p):
    syms = free_symbols_sparse(p)
    return min(syms)


def gcd_prs_sparse(p1, p2):
    x = main_variable_sparse(p1)

    c1, pp1 = primitive_sparse(p1, x)
    c2, pp2 = primitive_sparse(p2, x)

    h = subresultants_sparse(pp1, pp2, x)[-1]
    c = gcd_sparse([c1, c2])

    K = p1.ring.to_domain()
    if K.is_negative(coeff_sparse(h, x, h.degree(x))):
        h = -h

    _, h = primitive_sparse(h, x)
    h = h * c

    return h

"""Euclidean algorithms, GCDs, LCMs and polynomial remainder sequences. """

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

from sympy.polys.densearith import (
    dup_sub_mul,
    dup_neg, dmp_neg,
    dup_add, dmp_add,
    dup_sub, dmp_sub,
    dup_mul, dmp_mul,
    dup_pow, dmp_pow,
    dup_div, dmp_div,
    dup_rem, dmp_rem,
    dup_exquo, dmp_exquo,
    dup_prem, dmp_prem,
    dup_mul_ground, dmp_mul_ground,
    dup_mul_term, dmp_mul_term,
    dup_exquo_ground, dmp_exquo_ground,
    dup_max_norm, dmp_max_norm)

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

from sympy.polys.polyerrors import (
    MultivariatePolynomialError,
    HeuristicGCDFailed,
    HomomorphismFailed,
    NotInvertible,
    DomainError)

from sympy.polys.polyconfig import query

from sympy.utilities import cythonized

from sympy.ntheory import nextprime

def dup_half_gcdex(f, g, K):
    """
    Half extended Euclidean algorithm in ``F[x]``.

    Returns ``(s, h)`` such that ``h = gcd(f, g)`` and ``s*f = h (mod g)``.

    Example
    =======

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.euclidtools import dup_half_gcdex

    >>> f = QQ.map([1, -2, -6, 12, 15])
    >>> g = QQ.map([1, 1, -4, -4])

    >>> dup_half_gcdex(f, g, QQ)
    ([-1/5, 3/5], [1/1, 1/1])

    """
    if not (K.has_Field or not K.is_Exact):
        raise DomainError("can't compute half extended GCD over %s" % K)

    a, b = [K.one], []

    while g:
        q, r = dup_div(f, g, K)
        f, g = g, r
        a, b = b, dup_sub_mul(a, q, b, K)

    a = dup_exquo_ground(a, dup_LC(f, K), K)
    f = dup_monic(f, K)

    return a, f

def dmp_half_gcdex(f, g, u, K):
    """
    Half extended Euclidean algorithm in ``F[X]``.

    Example
    =======

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.euclidtools import dmp_half_gcdex

    """
    if not u:
        return dup_half_gcdex(f, g, K)
    else:
        raise MultivariatePolynomialError(f, g)

def dup_gcdex(f, g, K):
    """
    Extended Euclidean algorithm in ``F[x]``.

    Returns ``(s, t, h)`` such that ``h = gcd(f, g)`` and ``s*f + t*g = h``.

    Example
    =======

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.euclidtools import dup_gcdex

    >>> f = QQ.map([1, -2, -6, 12, 15])
    >>> g = QQ.map([1, 1, -4, -4])

    >>> dup_gcdex(f, g, QQ)
    ([-1/5, 3/5], [1/5, -6/5, 2/1], [1/1, 1/1])

    """
    s, h = dup_half_gcdex(f, g, K)

    F = dup_sub_mul(h, s, f, K)
    t = dup_exquo(F, g, K)

    return s, t, h

def dmp_gcdex(f, g, u, K):
    """
    Extended Euclidean algorithm in ``F[X]``.

    Example
    =======

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.euclidtools import dmp_gcdex

    """
    if not u:
        return dup_gcdex(f, g, K)
    else:
        raise MultivariatePolynomialError(f, g)

def dup_invert(f, g, K):
    """
    Compute multiplicative inverse of ``f`` modulo ``g`` in ``F[x]``.

    Example
    =======

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.euclidtools import dup_invert

    >>> f = QQ.map([1, 0, -1])
    >>> g = QQ.map([2, -1])
    >>> h = QQ.map([1, -1])

    >>> dup_invert(f, g, QQ)
    [-4/3]

    >>> dup_invert(f, h, QQ)
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
    Compute multiplicative inverse of ``f`` modulo ``g`` in ``F[X]``.

    Example
    =======

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.euclidtools import dmp_invert

    """
    if not u:
        return dup_invert(f, g, K)
    else:
        raise MultivariatePolynomialError(f, g)

def dup_euclidean_prs(f, g, K):
    """
    Euclidean polynomial remainder sequence (PRS) in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.euclidtools import dup_euclidean_prs

    >>> f = QQ.map([1, 0, 1, 0, -3, -3, 8, 2, -5])
    >>> g = QQ.map([3, 0, 5, 0, -4, -9, 21])

    >>> prs = dup_euclidean_prs(f, g, QQ)

    >>> prs[0]
    [1/1, 0/1, 1/1, 0/1, -3/1, -3/1, 8/1, 2/1, -5/1]
    >>> prs[1]
    [3/1, 0/1, 5/1, 0/1, -4/1, -9/1, 21/1]
    >>> prs[2]
    [-5/9, 0/1, 1/9, 0/1, -1/3]
    >>> prs[3]
    [-117/25, -9/1, 441/25]
    >>> prs[4]
    [233150/19773, -102500/6591]
    >>> prs[5]
    [-1288744821/543589225]

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
    Euclidean polynomial remainder sequence (PRS) in ``K[X]``.

    Example
    =======

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.euclidtools import dmp_euclidean_prs

    """
    if not u:
        return dup_euclidean_prs(f, g, K)
    else:
        raise MultivariatePolynomialError(f, g)

def dup_primitive_prs(f, g, K):
    """
    Primitive polynomial remainder sequence (PRS) in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dup_primitive_prs

    >>> f = ZZ.map([1, 0, 1, 0, -3, -3, 8, 2, -5])
    >>> g = ZZ.map([3, 0, 5, 0, -4, -9, 21])

    >>> prs = dup_primitive_prs(f, g, ZZ)

    >>> prs[0]
    [1, 0, 1, 0, -3, -3, 8, 2, -5]
    >>> prs[1]
    [3, 0, 5, 0, -4, -9, 21]
    >>> prs[2]
    [-5, 0, 1, 0, -3]
    >>> prs[3]
    [13, 25, -49]
    >>> prs[4]
    [4663, -6150]
    >>> prs[5]
    [1]

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
    Primitive polynomial remainder sequence (PRS) in ``K[X]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dmp_primitive_prs

    """
    if not u:
        return dup_primitive_prs(f, g, K)
    else:
        raise MultivariatePolynomialError(f, g)

@cythonized("n,m,d,k")
def dup_inner_subresultants(f, g, K):
    """
    Subresultant PRS algorithm in ``K[x]``.

    Computes the subresultant polynomial remainder sequence (PRS) of ``f``
    and ``g``, and the values for $\\beta_i$ and $\\delta_i$. The last two
    sequences of values are necessary for computing the resultant in
    :func:`dup_prs_resultant`.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dup_inner_subresultants

    >>> f = ZZ.map([1, 0, 1])
    >>> g = ZZ.map([1, 0, -1])

    >>> dup_inner_subresultants(f, g, ZZ)
    ([[1, 0, 1], [1, 0, -1], [-2]], [-1, -1], [0, 2])

    """
    n = dup_degree(f)
    m = dup_degree(g)

    if n < m:
        f, g = g, f
        n, m = m, n

    R = [f, g]
    d = n - m

    b = (-K.one)**(d+1)
    c =  -K.one

    B, D = [b], [d]

    if not f or not g:
        return R, B, D

    h = dup_prem(f, g, K)
    h = dup_mul_ground(h, b, K)

    while h:
        k = dup_degree(h)
        R.append(h)

        lc = dup_LC(g, K)

        if not d:
            q = c
        else:
            q = c**(d-1)

        c = K.exquo((-lc)**d, q)
        b = -lc * c**(m-k)

        f, g, m, d = g, h, k, m-k

        B.append(b)
        D.append(d)

        h = dup_prem(f, g, K)
        h = dup_exquo_ground(h, b, K)

    return R, B, D

def dup_subresultants(f, g, K):
    """
    Computes subresultant PRS of two polynomials in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dup_subresultants

    >>> f = ZZ.map([1, 0, 1])
    >>> g = ZZ.map([1, 0, -1])

    >>> dup_subresultants(f, g, ZZ)
    [[1, 0, 1], [1, 0, -1], [-2]]

    """
    return dup_inner_subresultants(f, g, K)[0]

@cythonized("s,i,du,dv,dw")
def dup_prs_resultant(f, g, K):
    """
    Resultant algorithm in ``K[x]`` using subresultant PRS.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dup_prs_resultant

    >>> f = ZZ.map([1, 0, 1])
    >>> g = ZZ.map([1, 0, -1])

    >>> dup_prs_resultant(f, g, ZZ)
    (4, [[1, 0, 1], [1, 0, -1], [-2]])

    """
    if not f or not g:
        return (K.zero, [])

    R, B, D = dup_inner_subresultants(f, g, K)

    if dup_degree(R[-1]) > 0:
        return (K.zero, R)
    if R[-2] == [K.one]:
        return (dup_LC(R[-1], K), R)

    s, i = 1, 1
    p, q = K.one, K.one

    for b, d in zip(B, D)[:-1]:
        du = dup_degree(R[i-1])
        dv = dup_degree(R[i  ])
        dw = dup_degree(R[i+1])

        if du % 2 and dv % 2:
            s = -s

        lc, i = dup_LC(R[i], K), i+1

        p *= b**dv * lc**(du-dw)
        q *= lc**(dv*(1+d))

    if s < 0:
        p = -p

    i = dup_degree(R[-2])

    res = dup_LC(R[-1], K)**i
    res = K.exquo(res*p, q)

    return res, R

def dup_resultant(f, g, K):
    """
    Computes resultant of two polynomials in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dup_resultant

    >>> f = ZZ.map([1, 0, 1])
    >>> g = ZZ.map([1, 0, -1])

    >>> dup_resultant(f, g, ZZ)
    4

    """
    return dup_prs_resultant(f, g, K)[0]

@cythonized("u,v,n,m,d,k")
def dmp_inner_subresultants(f, g, u, K):
    """
    Subresultant PRS algorithm in ``K[X]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dmp_inner_subresultants

    >>> f = ZZ.map([[3, 0], [], [-1, 0, 0, -4]])
    >>> g = ZZ.map([[1], [1, 0, 0, 0], [-9]])

    >>> a = [[3, 0, 0, 0, 0], [1, 0, -27, 4]]
    >>> b = [[-3, 0, 0, -12, 1, 0, -54, 8, 729, -216, 16]]

    >>> R = ZZ.map([f, g, a, b])
    >>> B = ZZ.map([[-1], [1], [9, 0, 0, 0, 0, 0, 0, 0, 0]])
    >>> D = ZZ.map([0, 1, 1])

    >>> dmp_inner_subresultants(f, g, 1, ZZ) == (R, B, D)
    True

    """
    if not u:
        return dup_inner_subresultants(f, g, K)

    n = dmp_degree(f, u)
    m = dmp_degree(g, u)

    if n < m:
        f, g = g, f
        n, m = m, n

    R = [f, g]
    d = n - m
    v = u - 1

    b = dmp_pow(dmp_ground(-K.one, v), d+1, v, K)
    c = dmp_ground(-K.one, v)

    B, D = [b], [d]

    if dmp_zero_p(f, u) or dmp_zero_p(g, u):
        return R, B, D

    h = dmp_prem(f, g, u, K)
    h = dmp_mul_term(h, b, 0, u, K)

    while not dmp_zero_p(h, u):
        k = dmp_degree(h, u)
        R.append(h)

        lc = dmp_LC(g, K)

        p = dmp_pow(dmp_neg(lc, v, K), d, v, K)

        if not d:
            q = c
        else:
            q = dmp_pow(c, d-1, v, K)

        c = dmp_exquo(p, q, v, K)
        b = dmp_mul(dmp_neg(lc, v, K),
                    dmp_pow(c, m-k, v, K), v, K)

        f, g, m, d = g, h, k, m-k

        B.append(b)
        D.append(d)

        h = dmp_prem(f, g, u, K)
        h = [ dmp_exquo(ch, b, v, K) for ch in h ]

    return R, B, D

@cythonized("u")
def dmp_subresultants(f, g, u, K):
    """
    Computes subresultant PRS of two polynomials in ``K[X]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dmp_subresultants

    >>> f = [[3, 0], [], [-1, 0, 0, -4]]
    >>> g = [[1], [1, 0, 0, 0], [-9]]

    >>> a = [[3, 0, 0, 0, 0], [1, 0, -27, 4]]
    >>> b = [[-3, 0, 0, -12, 1, 0, -54, 8, 729, -216, 16]]

    >>> dmp_subresultants(f, g, 1, ZZ) == [f, g, a, b]
    True

    """
    return dmp_inner_subresultants(f, g, u, K)[0]

@cythonized("u,v,s,i,d,du,dv,dw")
def dmp_prs_resultant(f, g, u, K):
    """
    Resultant algorithm in ``K[X]`` using subresultant PRS.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dmp_prs_resultant

    >>> f = ZZ.map([[3, 0], [], [-1, 0, 0, -4]])
    >>> g = ZZ.map([[1], [1, 0, 0, 0], [-9]])

    >>> a = ZZ.map([[3, 0, 0, 0, 0], [1, 0, -27, 4]])
    >>> b = ZZ.map([[-3, 0, 0, -12, 1, 0, -54, 8, 729, -216, 16]])

    >>> dmp_prs_resultant(f, g, 1, ZZ) == (b[0], [f, g, a, b])
    True

    """
    if not u:
        return dup_prs_resultant(f, g, K)

    if dmp_zero_p(f, u) or dmp_zero_p(g, u):
        return (dmp_zero(u-1), [])

    R, B, D = dmp_inner_subresultants(f, g, u, K)

    if dmp_degree(R[-1], u) > 0:
        return (dmp_zero(u-1), R)
    if dmp_one_p(R[-2], u, K):
        return (dmp_LC(R[-1], K), R)

    s, i, v = 1, 1, u-1

    p = dmp_one(v, K)
    q = dmp_one(v, K)

    for b, d in zip(B, D)[:-1]:
        du = dmp_degree(R[i-1], u)
        dv = dmp_degree(R[i  ], u)
        dw = dmp_degree(R[i+1], u)

        if du % 2 and dv % 2:
            s = -s

        lc, i = dmp_LC(R[i], K), i+1

        p = dmp_mul(dmp_mul(p, dmp_pow(b, dv, v, K), v, K),
                               dmp_pow(lc, du-dw, v, K), v, K)
        q = dmp_mul(q, dmp_pow(lc, dv*(1+d), v, K), v, K)

        _, p, q = dmp_inner_gcd(p, q, v, K)

    if s < 0:
        p = dmp_neg(p, v, K)

    i = dmp_degree(R[-2], u)

    res = dmp_pow(dmp_LC(R[-1], K), i, v, K)
    res = dmp_exquo(dmp_mul(res, p, v, K), q, v, K)

    return res, R

@cythonized("u,v,n,m,N,M,B")
def dmp_zz_modular_resultant(f, g, p, u, K):
    """
    Compute resultant of ``f`` and ``g`` modulo a prime ``p``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dmp_zz_modular_resultant

    >>> f = ZZ.map([[1], [1, 2]])
    >>> g = ZZ.map([[2, 1], [3]])

    >>> dmp_zz_modular_resultant(f, g, ZZ(5), 1, ZZ)
    [-2, 0, 1]

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

@cythonized("u,v,n,m")
def dmp_zz_collins_resultant(f, g, u, K):
    """
    Collins's modular resultant algorithm in ``Z[X]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dmp_zz_collins_resultant

    >>> f = ZZ.map([[1], [1, 2]])
    >>> g = ZZ.map([[2, 1], [3]])

    >>> dmp_zz_collins_resultant(f, g, 1, ZZ)
    [-2, -5, 1]

    """

    n = dmp_degree(f, u)
    m = dmp_degree(g, u)

    if n < 0 or m < 0:
        return dmp_zero(u-1)

    A = dmp_max_norm(f, u, K)
    B = dmp_max_norm(g, u, K)

    a = dmp_ground_LC(f, u, K)
    b = dmp_ground_LC(g, u, K)

    v = u - 1

    B = K(2)*K.factorial(n+m)*A**m*B**n
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

@cythonized("u,n,m")
def dmp_qq_collins_resultant(f, g, u, K0):
    """
    Collins's modular resultant algorithm in ``Q[X]``.

    Example
    =======

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.euclidtools import dmp_qq_collins_resultant

    >>> f = [[QQ(1,2)], [QQ(1), QQ(2,3)]]
    >>> g = [[QQ(2), QQ(1)], [QQ(3)]]

    >>> dmp_qq_collins_resultant(f, g, 1, QQ)
    [-2/1, -7/3, 5/6]

    """
    n = dmp_degree(f, u)
    m = dmp_degree(g, u)

    if n < 0 or m < 0:
        return dmp_zero(u-1)

    K1 = K0.get_ring()

    cf, f = dmp_clear_denoms(f, u, K0, K1)
    cg, g = dmp_clear_denoms(g, u, K0, K1)

    f = dmp_convert(f, u, K0, K1)
    g = dmp_convert(g, u, K0, K1)

    r = dmp_zz_collins_resultant(f, g, u, K1)
    r = dmp_convert(r, u-1, K1, K0)

    c = K0.convert(cf**m * cg**n, K1)

    return dmp_exquo_ground(r, c, u-1, K0)

@cythonized("u")
def dmp_resultant(f, g, u, K):
    """
    Computes resultant of two polynomials in ``K[X]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dmp_resultant

    >>> f = ZZ.map([[3, 0], [], [-1, 0, 0, -4]])
    >>> g = ZZ.map([[1], [1, 0, 0, 0], [-9]])

    >>> dmp_resultant(f, g, 1, ZZ)
    [-3, 0, 0, -12, 1, 0, -54, 8, 729, -216, 16]

    """
    if not u:
        return dup_resultant(f, g, K)

    if K.has_Field:
        if K.is_QQ and query('USE_COLLINS_RESULTANT'):
            return dmp_qq_collins_resultant(f, g, u, K)
    else:
        if K.is_ZZ and query('USE_COLLINS_RESULTANT'):
            return dmp_zz_collins_resultant(f, g, u, K)

    return dmp_prs_resultant(f, g, u, K)[0]

@cythonized("d,s")
def dup_discriminant(f, K):
    """
    Computes discriminant of a polynomial in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dup_discriminant

    >>> dup_discriminant([ZZ(1), ZZ(2), ZZ(3)], ZZ)
    -8

    """
    d = dup_degree(f)

    if d <= 0:
        return K.zero
    else:
        s = (-1)**((d*(d-1)) // 2)
        c = dup_LC(f, K)

        r = dup_resultant(f, dup_diff(f, 1, K), K)

        return K.exquo(r, c*K(s))

@cythonized("u,v,d,s")
def dmp_discriminant(f, u, K):
    """
    Computes discriminant of a polynomial in ``K[X]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dmp_discriminant

    >>> f = ZZ.map([[[[1]], [[]]], [[[1], []]], [[[1, 0]]]])

    >>> dmp_discriminant(f, 3, ZZ)
    [[[-4, 0]], [[1], [], []]]

    """
    if not u:
        return dup_discriminant(f, K)

    d, v = dmp_degree(f, u), u-1

    if d <= 0:
        return dmp_zero(v)
    else:
        s = (-1)**((d*(d-1)) // 2)
        c = dmp_LC(f, K)

        r = dmp_resultant(f, dmp_diff(f, 1, u, K), u, K)
        c = dmp_mul_ground(c, K(s), v, K)

        return dmp_exquo(r, c, v, K)

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

@cythonized("u")
def _dmp_rr_trivial_gcd(f, g, u, K):
    """Handle trivial cases in GCD algorithm over a ring. """
    zero_f = dmp_zero_p(f, u)
    zero_g = dmp_zero_p(g, u)

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
    elif query('USE_SIMPLIFY_GCD'):
        return _dmp_simplify_gcd(f, g, u, K)
    else:
        return None

@cythonized("u")
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

@cythonized("u,v,df,dg")
def _dmp_simplify_gcd(f, g, u, K):
    """Try to eliminate ``x_0`` from GCD computation in ``K[X]``. """
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

    cff = [ dmp_exquo(cf, h, v, K) for cf in f ]
    cfg = [ dmp_exquo(cg, h, v, K) for cg in g ]

    return [h], cff, cfg

def dup_rr_prs_gcd(f, g, K):
    """
    Computes polynomial GCD using subresultants over a ring.

    Returns ``(h, cff, cfg)`` such that ``a = gcd(f, g)``, ``cff = quo(f, h)``,
    and ``cfg = quo(g, h)``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dup_rr_prs_gcd

    >>> f = ZZ.map([1, 0, -1])
    >>> g = ZZ.map([1, -3, 2])

    >>> dup_rr_prs_gcd(f, g, ZZ)
    ([1, -1], [1, 1], [1, -2])

    """
    result = _dup_rr_trivial_gcd(f, g, K)

    if result is not None:
        return result

    fc, F = dup_primitive(f, K)
    gc, G = dup_primitive(g, K)

    c = K.gcd(fc, gc)

    h = dup_subresultants(F, G, K)[-1]
    _, h = dup_primitive(h, K)

    if K.is_negative(dup_LC(h, K)):
        c = -c

    h = dup_mul_ground(h, c, K)

    cff = dup_exquo(f, h, K)
    cfg = dup_exquo(g, h, K)

    return h, cff, cfg

def dup_ff_prs_gcd(f, g, K):
    """
    Computes polynomial GCD using subresultants over a field.

    Returns ``(h, cff, cfg)`` such that ``a = gcd(f, g)``, ``cff = quo(f, h)``,
    and ``cfg = quo(g, h)``.

    Example
    =======

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.euclidtools import dup_ff_prs_gcd

    >>> f = QQ.map([1, 0, -1])
    >>> g = QQ.map([1, -3, 2])

    >>> dup_ff_prs_gcd(f, g, QQ)
    ([1/1, -1/1], [1/1, 1/1], [1/1, -2/1])

    """
    result = _dup_ff_trivial_gcd(f, g, K)

    if result is not None:
        return result

    h = dup_subresultants(f, g, K)[-1]
    h = dup_monic(h, K)

    cff = dup_exquo(f, h, K)
    cfg = dup_exquo(g, h, K)

    return h, cff, cfg

@cythonized("u")
def dmp_rr_prs_gcd(f, g, u, K):
    """
    Computes polynomial GCD using subresultants over a ring.

    Returns ``(h, cff, cfg)`` such that ``a = gcd(f, g)``, ``cff = quo(f, h)``,
    and ``cfg = quo(g, h)``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dmp_rr_prs_gcd

    >>> f = ZZ.map([[1], [2, 0], [1, 0, 0]])
    >>> g = ZZ.map([[1], [1, 0], []])

    >>> dmp_rr_prs_gcd(f, g, 1, ZZ)
    ([[1], [1, 0]], [[1], [1, 0]], [[1], []])

    """
    if not u:
        return dup_rr_prs_gcd(f, g, K)

    result = _dmp_rr_trivial_gcd(f, g, u, K)

    if result is not None:
        return result

    fc, F = dmp_primitive(f, u, K)
    gc, G = dmp_primitive(g, u, K)

    h = dmp_subresultants(F, G, u, K)[-1]
    c, _, _ = dmp_rr_prs_gcd(fc, gc, u-1, K)

    if K.is_negative(dmp_ground_LC(h, u, K)):
        h = dmp_neg(h, u, K)

    _, h = dmp_primitive(h, u, K)
    h = dmp_mul_term(h, c, 0, u, K)

    cff = dmp_exquo(f, h, u, K)
    cfg = dmp_exquo(g, h, u, K)

    return h, cff, cfg

@cythonized("u")
def dmp_ff_prs_gcd(f, g, u, K):
    """
    Computes polynomial GCD using subresultants over a field.

    Returns ``(h, cff, cfg)`` such that ``a = gcd(f, g)``, ``cff = quo(f, h)``,
    and ``cfg = quo(g, h)``.

    Example
    =======

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.euclidtools import dmp_ff_prs_gcd

    >>> f = [[QQ(1,2)], [QQ(1), QQ(0)], [QQ(1,2), QQ(0), QQ(0)]]
    >>> g = [[QQ(1)], [QQ(1), QQ(0)], []]

    >>> dmp_ff_prs_gcd(f, g, 1, QQ)
    ([[1/1], [1/1, 0/1]], [[1/2], [1/2, 0/1]], [[1/1], []])

    """
    if not u:
        return dup_ff_prs_gcd(f, g, K)

    result = _dmp_ff_trivial_gcd(f, g, u, K)

    if result is not None:
        return result

    fc, f = dmp_primitive(f, u, K)
    gc, g = dmp_primitive(g, u, K)

    h = dmp_subresultants(f, g, u, K)[-1]
    c, _, _ = dmp_ff_prs_gcd(fc, gc, u-1, K)

    _, h = dmp_primitive(h, u, K)
    h = dmp_mul_term(h, c, 0, u, K)
    h = dmp_ground_monic(h, u, K)

    cff = dmp_exquo(f, h, u, K)
    cfg = dmp_exquo(g, h, u, K)

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
        h = (h-g) // x

    return f

@cythonized("i,df,dg")
def dup_zz_heu_gcd(f, g, K):
    """
    Heuristic polynomial GCD in ``Z[x]``.

    Given univariate polynomials ``f`` and ``g`` in ``Z[x]``, returns
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

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dup_zz_heu_gcd

    >>> f = ZZ.map([1, 0, -1])
    >>> g = ZZ.map([1, -3, 2])

    >>> dup_zz_heu_gcd(f, g, ZZ)
    ([1, -1], [1, 1], [1, -2])

    References
    ==========

    .. [Liao95] Hsin-Chao Liao,  R. Fateman, Evaluation of the heuristic
    polynomial GCD, International Symposium on Symbolic and Algebraic
    Computation (ISSAC), ACM Press, Montreal, Quebec, Canada, 1995,
    pp. 240--247

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

    B = 2*min(f_norm, g_norm) + 29

    x = max(min(B, 99*K.sqrt(B)),
            2*min(f_norm // abs(dup_LC(f, K)),
                  g_norm // abs(dup_LC(g, K))) + 2)

    for i in xrange(0, HEU_GCD_MAX):
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
                    return h, cff, cfg

        x = 73794*x * K.sqrt(K.sqrt(x)) // 27011

    raise HeuristicGCDFailed('no luck')

@cythonized("v")
def _dmp_zz_gcd_interpolate(h, x, v, K):
    """Interpolate polynomial GCD from integer GCD. """
    f = []

    while not dmp_zero_p(h, v):
        g = dmp_ground_trunc(h, x, v, K)
        f.insert(0, g)

        h = dmp_sub(h, g, v, K)
        h = dmp_exquo_ground(h, x, v, K)

    if K.is_negative(dmp_ground_LC(f, v+1, K)):
        return dmp_neg(f, v+1, K)
    else:
        return f

@cythonized("u,v,i,dg,df")
def dmp_zz_heu_gcd(f, g, u, K):
    """
    Heuristic polynomial GCD in ``Z[X]``.

    Given univariate polynomials ``f`` and ``g`` in ``Z[X]``, returns
    their GCD and cofactors, i.e. polynomials ``h``, ``cff`` and ``cfg``
    such that::

          h = gcd(f, g), cff = quo(f, h) and cfg = quo(g, h)

    The algorithm is purely heuristic which means it may fail to compute
    the GCD. This will be signaled by raising an exception. In this case
    you will need to switch to another GCD method.

    The algorithm computes the polynomial GCD by evaluating polynomials
    f and g at certain points and computing (fast) integer GCD of those
    evaluations. The polynomial GCD is recovered from the integer image
    by interpolation. The evaluation proces reduces f and g variable by
    variable into a large integer.  The final step  is to verify if the
    interpolated polynomial is the correct GCD. This gives cofactors of
    the input polynomials as a side effect.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dmp_zz_heu_gcd

    >>> f = ZZ.map([[1], [2, 0], [1, 0, 0]])
    >>> g = ZZ.map([[1], [1, 0], []])

    >>> dmp_zz_heu_gcd(f, g, 1, ZZ)
    ([[1], [1, 0]], [[1], [1, 0]], [[1], []])

    References
    ==========

    .. [Liao95] Hsin-Chao Liao,  R. Fateman, Evaluation of the heuristic
    polynomial GCD, International Symposium on Symbolic and Algebraic
    Computation (ISSAC), ACM Press, Montreal, Quebec, Canada, 1995,
    pp. 240--247

    """
    if not u:
        return dup_zz_heu_gcd(f, g, K)

    result = _dmp_rr_trivial_gcd(f, g, u, K)

    if result is not None:
        return result

    df = dmp_degree(f, u)
    dg = dmp_degree(g, u)

    gcd, f, g = dmp_ground_extract(f, g, u, K)

    f_norm = dmp_max_norm(f, u, K)
    g_norm = dmp_max_norm(g, u, K)

    B = 2*min(f_norm, g_norm) + 29

    x = max(min(B, 99*K.sqrt(B)),
            2*min(f_norm // abs(dmp_ground_LC(f, u, K)),
                  g_norm // abs(dmp_ground_LC(g, u, K))) + 2)

    for i in xrange(0, HEU_GCD_MAX):
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
    Heuristic polynomial GCD in ``Q[x]``.

    Returns ``(h, cff, cfg)`` such that ``a = gcd(f, g)``,
    ``cff = quo(f, h)``, and ``cfg = quo(g, h)``.

    Example
    =======

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.euclidtools import dup_qq_heu_gcd

    >>> f = [QQ(1,2), QQ(7,4), QQ(3,2)]
    >>> g = [QQ(1,2), QQ(1), QQ(0)]

    >>> dup_qq_heu_gcd(f, g, QQ)
    ([1/1, 2/1], [1/2, 3/4], [1/2, 0/1])

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

@cythonized("u")
def dmp_qq_heu_gcd(f, g, u, K0):
    """
    Heuristic polynomial GCD in ``Q[X]``.

    Returns ``(h, cff, cfg)`` such that ``a = gcd(f, g)``,
    ``cff = quo(f, h)``, and ``cfg = quo(g, h)``.

    Example
    =======

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.euclidtools import dmp_qq_heu_gcd

    >>> f = [[QQ(1,4)], [QQ(1), QQ(0)], [QQ(1), QQ(0), QQ(0)]]
    >>> g = [[QQ(1,2)], [QQ(1), QQ(0)], []]

    >>> dmp_qq_heu_gcd(f, g, 1, QQ)
    ([[1/1], [2/1, 0/1]], [[1/4], [1/2, 0/1]], [[1/2], []])

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
    Computes polynomial GCD and cofactors of ``f`` and ``g`` in ``K[x]``.

    Returns ``(h, cff, cfg)`` such that ``a = gcd(f, g)``,
    ``cff = quo(f, h)``, and ``cfg = quo(g, h)``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dup_inner_gcd

    >>> f = ZZ.map([1, 0, -1])
    >>> g = ZZ.map([1, -3, 2])

    >>> dup_inner_gcd(f, g, ZZ)
    ([1, -1], [1, 1], [1, -2])

    """
    if K.has_Field or not K.is_Exact:
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

@cythonized("u")
def _dmp_inner_gcd(f, g, u, K):
    """Helper function for `dmp_inner_gcd()`. """
    if K.has_Field or not K.is_Exact:
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

@cythonized("u")
def dmp_inner_gcd(f, g, u, K):
    """
    Computes polynomial GCD and cofactors of ``f`` and ``g`` in ``K[X]``.

    Returns ``(h, cff, cfg)`` such that ``a = gcd(f, g)``,
    ``cff = quo(f, h)``, and ``cfg = quo(g, h)``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dmp_inner_gcd

    >>> f = ZZ.map([[1], [2, 0], [1, 0, 0]])
    >>> g = ZZ.map([[1], [1, 0], []])

    >>> dmp_inner_gcd(f, g, 1, ZZ)
    ([[1], [1, 0]], [[1], [1, 0]], [[1], []])

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
    Computes polynomial GCD of ``f`` and ``g`` in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dup_gcd

    >>> f = ZZ.map([1, 0, -1])
    >>> g = ZZ.map([1, -3, 2])

    >>> dup_gcd(f, g, ZZ)
    [1, -1]

    """
    return dup_inner_gcd(f, g, K)[0]

@cythonized("u")
def dmp_gcd(f, g, u, K):
    """
    Computes polynomial GCD of ``f`` and ``g`` in ``K[X]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dmp_gcd

    >>> f = ZZ.map([[1], [2, 0], [1, 0, 0]])
    >>> g = ZZ.map([[1], [1, 0], []])

    >>> dmp_gcd(f, g, 1, ZZ)
    [[1], [1, 0]]

    """
    return dmp_inner_gcd(f, g, u, K)[0]

def dup_rr_lcm(f, g, K):
    """
    Computes polynomial LCM over a ring in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dup_rr_lcm

    >>> f = ZZ.map([1, 0, -1])
    >>> g = ZZ.map([1, -3, 2])

    >>> dup_rr_lcm(f, g, ZZ)
    [1, -2, -1, 2]

    """
    fc, f = dup_primitive(f, K)
    gc, g = dup_primitive(g, K)

    c = K.lcm(fc, gc)

    h = dup_exquo(dup_mul(f, g, K),
                  dup_gcd(f, g, K), K)

    return dup_mul_ground(h, c, K)

def dup_ff_lcm(f, g, K):
    """
    Computes polynomial LCM over a field in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.euclidtools import dup_ff_lcm

    >>> f = [QQ(1,2), QQ(7,4), QQ(3,2)]
    >>> g = [QQ(1,2), QQ(1), QQ(0)]

    >>> dup_ff_lcm(f, g, QQ)
    [1/1, 7/2, 3/1, 0/1]

    """
    h = dup_exquo(dup_mul(f, g, K),
                  dup_gcd(f, g, K), K)

    return dup_monic(h, K)

def dup_lcm(f, g, K):
    """
    Computes polynomial LCM of ``f`` and ``g`` in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dup_lcm

    >>> f = ZZ.map([1, 0, -1])
    >>> g = ZZ.map([1, -3, 2])

    >>> dup_lcm(f, g, ZZ)
    [1, -2, -1, 2]

    """
    if K.has_Field or not K.is_Exact:
        return dup_ff_lcm(f, g, K)
    else:
        return dup_rr_lcm(f, g, K)

@cythonized("u")
def dmp_rr_lcm(f, g, u, K):
    """
    Computes polynomial LCM over a ring in ``K[X]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dmp_rr_lcm

    >>> f = ZZ.map([[1], [2, 0], [1, 0, 0]])
    >>> g = ZZ.map([[1], [1, 0], []])

    >>> dmp_rr_lcm(f, g, 1, ZZ)
    [[1], [2, 0], [1, 0, 0], []]

    """
    fc, f = dmp_ground_primitive(f, u, K)
    gc, g = dmp_ground_primitive(g, u, K)

    c = K.lcm(fc, gc)

    h = dmp_exquo(dmp_mul(f, g, u, K),
                  dmp_gcd(f, g, u, K), u, K)

    return dmp_mul_ground(h, c, u, K)

@cythonized("u")
def dmp_ff_lcm(f, g, u, K):
    """
    Computes polynomial LCM over a field in ``K[X]``.

    Example
    =======

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.euclidtools import dmp_ff_lcm

    >>> f = [[QQ(1,4)], [QQ(1), QQ(0)], [QQ(1), QQ(0), QQ(0)]]
    >>> g = [[QQ(1,2)], [QQ(1), QQ(0)], []]

    >>> dmp_ff_lcm(f, g, 1, QQ)
    [[1/1], [4/1, 0/1], [4/1, 0/1, 0/1], []]

    """
    h = dmp_exquo(dmp_mul(f, g, u, K),
                  dmp_gcd(f, g, u, K), u, K)

    return dmp_ground_monic(h, u, K)

@cythonized("u")
def dmp_lcm(f, g, u, K):
    """
    Computes polynomial LCM of ``f`` and ``g`` in ``K[X]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dmp_lcm

    >>> f = ZZ.map([[1], [2, 0], [1, 0, 0]])
    >>> g = ZZ.map([[1], [1, 0], []])

    >>> dmp_lcm(f, g, 1, ZZ)
    [[1], [2, 0], [1, 0, 0], []]

    """
    if not u:
        return dup_lcm(f, g, K)

    if K.has_Field or not K.is_Exact:
        return dmp_ff_lcm(f, g, u, K)
    else:
        return dmp_rr_lcm(f, g, u, K)

@cythonized("u,v")
def dmp_content(f, u, K):
    """
    Returns GCD of multivariate coefficients.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dmp_content

    >>> f = ZZ.map([[2, 6], [4, 12]])

    >>> dmp_content(f, 1, ZZ)
    [2, 6]

    """
    cont, v = dmp_LC(f, K), u-1

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

@cythonized("u,v")
def dmp_primitive(f, u, K):
    """
    Returns multivariate content and a primitive polynomial.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dmp_primitive

    >>> f = ZZ.map([[2, 6], [4, 12]])

    >>> dmp_primitive(f, 1, ZZ)
    ([2, 6], [[1], [2]])

    """
    cont, v = dmp_content(f, u, K), u-1

    if dmp_zero_p(f, u) or dmp_one_p(cont, v, K):
        return cont, f
    else:
        return cont, [ dmp_exquo(c, cont, v, K) for c in f ]

def dup_cancel(f, g, K, multout=True):
    """
    Cancel common factors in a rational function ``f/g``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dup_cancel

    >>> f = ZZ.map([2, 0, -2])
    >>> g = ZZ.map([1, -2, 1])

    >>> dup_cancel(f, g, ZZ)
    ([2, 2], [1, -1])

    """
    return dmp_cancel(f, g, 0, K, multout=multout)

def dmp_cancel(f, g, u, K, multout=True):
    """
    Cancel common factors in a rational function ``f/g``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.euclidtools import dmp_cancel

    >>> f = ZZ.map([[2], [0], [-2]])
    >>> g = ZZ.map([[1], [-2], [1]])

    >>> dmp_cancel(f, g, 1, ZZ)
    ([[2], [2]], [[1], [-1]])

    """
    if dmp_zero_p(f, u) or dmp_zero_p(g, u):
        if multout:
            return f, g
        else:
            return K.one, K.one, f, g

    K0 = None

    if K.has_Field and K.has_assoc_Ring:
        K0, K = K, K.get_ring()

        cq, f = dmp_clear_denoms(f, u, K0, K, convert=True)
        cp, g = dmp_clear_denoms(g, u, K0, K, convert=True)
    else:
        cp, cq = K.one, K.one

    _, p, q = dmp_inner_gcd(f, g, u, K)

    if K0 is not None:
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

    if not multout:
        return cp, cq, p, q

    p = dmp_mul_ground(p, cp, u, K)
    q = dmp_mul_ground(q, cq, u, K)

    return p, q

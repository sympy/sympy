"""Dense univariate polynomials with coefficients in Galois fields. """

from random import uniform
from math import ceil, sqrt, log

from sympy.polys.polyutils import _sort_factors
from sympy.polys.polyconfig import query

from sympy.polys.polyerrors import ExactQuotientFailed

from sympy.utilities import (
    any, all, cythonized
)

from sympy.ntheory import factorint

def gf_crt(U, M, K):
    """
    Chinese Remainder Theorem.

    Given a set of integer residues `u_0,...,u_n` and a set of
    co-prime integer moduli `m_0,...,m_n`, returns an integer
    `u`, such that `u = u_i mod m_i` for `i = `0,...,n`.

    As an example consider a set of residues `U = [49, 76, 65]`
    and a set of moduli `M = [99, 97, 95]`. Then we have::

       >>> from sympy.polys.galoistools import gf_crt
       >>> from sympy.polys.algebratools import ZZ

       >>> gf_crt([49, 76, 65], [99, 97, 95], ZZ)
       639985

    This is the correct result because::

       >>> 639985 % 99
       49

       >>> 639985 % 97
       76

       >>> 639985 % 95
       65

    """
    p, v = K.one, K.zero

    for m in M:
        p *= m

    for u, m in zip(U, M):
        e = p // m
        s, _, _ = K.gcdex(e, m)
        v += e*(u*s % m)

    return v % p

def gf_crt1(M, K):
    """
    First part of the Chinese Remainder Theorem.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_crt1
    >>> gf_crt1([99, 97, 95], ZZ)
    (912285, [9215, 9405, 9603], [62, 24, 12])
    """
    p, E, S = K.one, [], []

    for m in M:
        p *= m

    for m in M:
        E.append(p // m)
        S.append(K.gcdex(E[-1], m)[0] % m)

    return p, E, S

def gf_crt2(U, M, p, E, S, K):
    """
    Second part of the Chinese Remainder Theorem.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_crt2
    >>> gf_crt2([49, 76, 65], [99, 97, 95], 912285, [9215, 9405, 9603],
    ... [62, 24, 12], ZZ)
    639985
    """
    v = K.zero

    for u, m, e, s in zip(U, M, E, S):
        v += e*(u*s % m)

    return v % p

def gf_int(a, p):
    """
    Coerce `a mod p` to an integer in `[-p/2, p/2]` range.

    Example
    =======
    >>> from sympy.polys.galoistools import gf_int
    >>> gf_int(2, 7)
    2
    >>> gf_int(5, 7)
    -2
    """
    if a <= p // 2:
        return a
    else:
        return a - p

def gf_degree(f):
    """
    Returns leading degree of `f`.

    Example
    =======
    >>> from sympy.polys.galoistools import gf_degree
    >>> gf_degree([1, 1, 2, 0])
    3
    >>> gf_degree([])
    -1
    """
    return len(f)-1

def gf_LC(f, K):
    """
    Returns leading coefficient of `f`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_LC
    >>> gf_LC([3, 0, 1], ZZ)
    3
    """
    if not f:
        return K.zero
    else:
        return f[0]

def gf_TC(f, K):
    """
    Returns trailing coefficient of `f`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_TC
    >>> gf_TC([3, 0, 1], ZZ)
    1
    """
    if not f:
        return K.zero
    else:
        return f[-1]

@cythonized("k")
def gf_strip(f):
    """
    Remove leading zeros from `f`.

    Example
    =======
    >>> from sympy.polys.galoistools import gf_strip
    >>> gf_strip([0, 0, 0, 3, 0, 1])
    [3, 0, 1]
    """
    if not f or f[0]:
        return f

    k = 0

    for coeff in f:
        if coeff:
            break
        else:
            k += 1

    return f[k:]

def gf_trunc(f, p):
    """
    Reduce all coefficients modulo `p`.

    Example
    =======
    >>> from sympy.polys.galoistools import gf_trunc
    >>> gf_trunc([7, -2, 3], 5)
    [2, 3, 3]
    """
    return gf_strip([ a % p for a in f ])

def gf_normal(f, p, K):
    """
    Normalize all coefficients in `K`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_normal
    >>> gf_normal([5, 10, 21, -3], 5, ZZ)
    [1, 2]
    """
    return gf_trunc(map(K, f), p)

def gf_convert(f, p, K0, K1):
    """
    Normalize all coefficients in `K`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ, QQ
    >>> from sympy.polys.galoistools import gf_convert
    >>> gf_convert([QQ(1), QQ(4), QQ(0)], 3, QQ, ZZ)
    [1, 1, 0]
    """
    # XXX: The above doctest probably isn't representative of the behavior
    # this function was written for.  What is it really supposed to do?
    return gf_trunc([ K1.convert(c, K0) for c in f ], p)

@cythonized("k,n")
def gf_from_dict(f, p, K):
    """
    Create `GF(p)[x]` polynomial from a dict.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_from_dict
    >>> gf_from_dict({10:4, 4:33, 0:-1}, 5, ZZ)
    [4, 0, 0, 0, 0, 0, 3, 0, 0, 0, 4]
    """
    n, h = max(f.iterkeys()), []

    if type(n) is int:
        for k in xrange(n, -1, -1):
            h.append(f.get(k, K.zero) % p)
    else:
        (n,) = n

        for k in xrange(n, -1, -1):
            h.append(f.get((k,), K.zero) % p)

    return gf_trunc(h, p)

@cythonized("k,n")
def gf_to_dict(f, p, symmetric=True):
    """
    Convert `GF(p)[x]` polynomial to a dict.

    Example
    =======
    >>> from sympy.polys.galoistools import gf_to_dict
    >>> gf_to_dict([4, 0, 0, 0, 0, 0, 3, 0, 0, 0, 4], 5)
    {0: -1, 4: -2, 10: -1}
    >>> gf_to_dict([4, 0, 0, 0, 0, 0, 3, 0, 0, 0, 4], 5, symmetric=False)
    {0: 4, 4: 3, 10: 4}
    """
    n, result = gf_degree(f), {}

    for k in xrange(0, n+1):
        if symmetric:
            a = gf_int(f[n-k], p)
        else:
            a = f[n-k]

        if a: result[k] = a

    return result

def gf_from_int_poly(f, p):
    """
    Create `GF(p)[x]` polynomial from `Z[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_from_int_poly
    >>> gf_from_int_poly([7, -2, 3], 5)
    [2, 3, 3]
    """
    return gf_trunc(f, p)

def gf_to_int_poly(f, p, symmetric=True):
    """
    Convert `GF(p)[x]` polynomial to `Z[x]`.

    Example
    =======
    >>> from sympy.polys.galoistools import gf_to_int_poly
    >>> gf_to_int_poly([2, 3, 3], 5)
    [2, -2, -2]
    >>> gf_to_int_poly([2, 3, 3], 5, symmetric=False)
    [2, 3, 3]
    """
    if symmetric:
        return [ gf_int(c, p) for c in f ]
    else:
        return f

def gf_neg(f, p, K):
    """
    Negate a polynomial in `GF(p)[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_neg
    >>> gf_neg([3, 2, 1, 0], 5, ZZ)
    [2, 3, 4, 0]
    """
    return [ -coeff % p for coeff in f ]

def gf_add_ground(f, a, p, K):
    """
    Returns `f + a` where `f` in `GF(p)[x]` and `a` in `GF(p)`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_add_ground
    >>> gf_add_ground([3, 2, 4], 2, 5, ZZ)
    [3, 2, 1]
    """
    if not f:
        a = a % p
    else:
        a = (f[-1] + a) % p

        if len(f) > 1:
            return f[:-1] + [a]

    if not a:
        return []
    else:
        return [a]

def gf_sub_ground(f, a, p, K):
    """
    Returns `f - a` where `f` in `GF(p)[x]` and `a` in `GF(p)`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_sub_ground
    >>> gf_sub_ground([3, 2, 4], 2, 5, ZZ)
    [3, 2, 2]
    """
    if not f:
        a = -a % p
    else:
        a = (f[-1] - a) % p

        if len(f) > 1:
            return f[:-1] + [a]

    if not a:
        return []
    else:
        return [a]

def gf_mul_ground(f, a, p, K):
    """
    Returns `f * a` where `f` in `GF(p)[x]` and `a` in `GF(p)`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_mul_ground
    >>> gf_mul_ground([3, 2, 4], 2, 5, ZZ)
    [1, 4, 3]
    """
    if not a:
        return []
    else:
        return [ (a*b) % p for b in f ]

def gf_exquo_ground(f, a, p, K):
    """
    Returns `f / a` where `f` in `GF(p)[x]` and `a` in `GF(p)`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_exquo_ground
    >>> gf_exquo_ground([3, 2, 4], 2, 5, ZZ)
    [4, 1, 2]
    """
    return gf_mul_ground(f, K.invert(a, p), p, K)

@cythonized("df,dg,k")
def gf_add(f, g, p, K):
    """
    Add polynomials in `GF(p)[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_add
    >>> gf_add([3, 2, 4], [2, 2, 2], 5, ZZ)
    [4, 1]
    """
    if not f:
        return g
    if not g:
        return f

    df = gf_degree(f)
    dg = gf_degree(g)

    if df == dg:
        return gf_strip([ (a + b) % p for a, b in zip(f, g) ])
    else:
        k = abs(df - dg)

        if df > dg:
            h, f = f[:k], f[k:]
        else:
            h, g = g[:k], g[k:]

        return h + [ (a + b) % p for a, b in zip(f, g) ]

@cythonized("df,dg,k")
def gf_sub(f, g, p, K):
    """
    Subtract polynomials in `GF(p)[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_sub
    >>> gf_sub([3, 2, 4], [2, 2, 2], 5, ZZ)
    [1, 0, 2]
    """
    if not g:
        return f
    if not f:
        return gf_neg(g, p, K)

    df = gf_degree(f)
    dg = gf_degree(g)

    if df == dg:
        return gf_strip([ (a - b) % p for a, b in zip(f, g) ])
    else:
        k = abs(df - dg)

        if df > dg:
            h, f = f[:k], f[k:]
        else:
            h, g = gf_neg(g[:k], p, K), g[k:]

        return h + [ (a - b) % p for a, b in zip(f, g) ]

@cythonized("df,dg,dh,i,j")
def gf_mul(f, g, p, K):
    """
    Multiply polynomials in `GF(p)[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_mul
    >>> gf_mul([3, 2, 4], [2, 2, 2], 5, ZZ)
    [1, 0, 3, 2, 3]
    """
    df = gf_degree(f)
    dg = gf_degree(g)

    dh = df + dg
    h = [0]*(dh+1)

    for i in xrange(0, dh+1):
        coeff = K.zero

        for j in xrange(max(0, i-dg), min(i, df)+1):
            coeff += f[j]*g[i-j]

        h[i] = coeff % p

    return gf_strip(h)

@cythonized("df,dh,i,j,jmin,jmax,n")
def gf_sqr(f, p, K):
    """
    Square polynomials in `GF(p)[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_sqr
    >>> gf_sqr([3, 2, 4], 5, ZZ)
    [4, 2, 3, 1, 1]
    """
    df = gf_degree(f)

    dh = 2*df
    h = [0]*(dh+1)

    for i in xrange(0, dh+1):
        coeff = K.zero

        jmin = max(0, i-df)
        jmax = min(i, df)

        n = jmax - jmin + 1

        jmax = jmin + n // 2 - 1

        for j in xrange(jmin, jmax+1):
            coeff += f[j]*f[i-j]

        coeff += coeff

        if n & 1:
            elem = f[jmax+1]
            coeff += elem**2

        h[i] = coeff % p

    return gf_strip(h)

def gf_add_mul(f, g, h, p, K):
    """
    Returns `f + g*h` where `f`, `g`, `h` in `GF(p)[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_add_mul
    >>> gf_add_mul([3, 2, 4], [2, 2, 2], [1, 4], 5, ZZ)
    [2, 3, 2, 2]
    """
    return gf_add(f, gf_mul(g, h, p, K), p, K)

def gf_sub_mul(f, g, h, p, K):
    """
    Returns `f - g*h` where `f`, `g`, `h` in `GF(p)[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_sub_mul
    >>> gf_sub_mul([3, 2, 4], [2, 2, 2], [1, 4], 5, ZZ)
    [3, 3, 2, 1]
    """
    return gf_sub(f, gf_mul(g, h, p, K), p, K)

@cythonized("k")
def gf_expand(F, p, K):
    """
    Expand results of `factor()` in `GF(p)[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_expand
    >>> gf_expand([([3, 2, 4], 1), ([2, 2], 2), ([3, 1], 3)], 5, ZZ)
    [4, 3, 0, 3, 0, 1, 4, 1]
    """
    if type(F) is tuple:
        lc, F = F
    else:
        lc = K.one

    g = [lc]

    for f, k in F:
        f = gf_pow(f, k, p, K)
        g = gf_mul(g, f, p, K)

    return g

@cythonized("df,dg,dq,dr,i,j")
def gf_div(f, g, p, K):
    """
    Division with remainder in `GF(p)[x]`.

    Given univariate polynomials `f` and `g` with coefficients in a
    finite field with `p` elements, returns polynomials `q` and `r`
    (quotient and remainder) such that `f = q*g + r`.

    Consider polynomials `x**3 + x + 1` and `x**2 + x` in GF(2)::

       >>> from sympy.polys.galoistools import gf_div, gf_add_mul
       >>> from sympy.polys.algebratools import ZZ

       >>> gf_div([1, 0, 1, 1], [1, 1, 0], 2, ZZ)
       ([1, 1], [1])

    As result we obtained quotient `x + 1` and remainder `1`, thus::

       >>> gf_add_mul([1], [1, 1], [1, 1, 0], 2, ZZ)
       [1, 0, 1, 1]

    References
    ==========

    .. [Monagan93] Michael Monagan, In-place Arithmetic for Polynomials
       over Z_n, Proceedings of DISCO '92, Springer-Verlag LNCS, 721,
       1993, pp. 22-34

    .. [Gathen99] J. von zur Gathen, J. Gerhard, Modern Computer Algebra,
       First Edition, Cambridge University Press, 1999, pp. 247

    """
    df = gf_degree(f)
    dg = gf_degree(g)

    if not g:
        raise ZeroDivisionError("polynomial division")
    elif df < dg:
        return [], f

    inv = K.invert(g[0], p)

    h, dq, dr = list(f), df-dg, dg-1

    for i in xrange(0, df+1):
        coeff = h[i]

        for j in xrange(max(0, dg-i), min(df-i, dr)+1):
            coeff -= h[i+j-dg] * g[dg-j]

        if i <= dq:
            coeff *= inv

        h[i] = coeff % p

    return h[:dq+1], gf_strip(h[dq+1:])

def gf_rem(f, g, p, K):
    """
    Returns polynomial remainder in `GF(p)[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_rem
    >>> gf_rem([1, 0, 1, 1], [1, 1, 0], 2, ZZ)
    [1]
    """
    return gf_div(f, g, p, K)[1]

def gf_quo(f, g, p, K):
    """
    Returns polynomial quotient in `GF(p)[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_quo
    >>> gf_quo([1, 0, 1, 1], [1, 1, 0], 2, ZZ)
    Traceback (most recent call last):
    ...
    ExactQuotientFailed: [1, 1, 0] does not divide [1, 0, 1, 1] in ZZ
    >>> gf_quo([1, 0, 3, 2, 3], [2, 2, 2], 5, ZZ)
    [3, 2, 4]
    """
    q, r = gf_div(f, g, p, K)

    if not r:
        return q
    else:
        raise ExactQuotientFailed('%s does not divide %s in %s' % (g, f, K))

@cythonized("df,dg,dq,dr,i,j")
def gf_exquo(f, g, p, K):
    """
    Computes exact quotient in `GF(p)[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_exquo
    >>> gf_exquo([1, 0, 1, 1], [1, 1, 0], 2, ZZ)
    [1, 1]
    >>> gf_exquo([1, 0, 3, 2, 3], [2, 2, 2], 5, ZZ)
    [3, 2, 4]
    """
    df = gf_degree(f)
    dg = gf_degree(g)

    if not g:
        raise ZeroDivisionError("polynomial division")
    elif df < dg:
        return []

    inv = K.invert(g[0], p)

    h, dq, dr = f[:], df-dg, dg-1

    for i in xrange(0, dq+1):
        coeff = h[i]

        for j in xrange(max(0, dg-i), min(df-i, dr)+1):
            coeff -= h[i+j-dg] * g[dg-j]

        h[i] = (coeff * inv) % p

    return h[:dq+1]

@cythonized("n")
def gf_lshift(f, n, K):
    """
    Efficiently multiply `f` by `x**n`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_lshift
    >>> gf_lshift([3, 2, 4], 4, ZZ)
    [3, 2, 4, 0, 0, 0, 0]
    """
    if not f:
        return f
    else:
        return f + [K.zero]*n

@cythonized("n")
def gf_rshift(f, n, K):
    """
    Efficiently divide `f` by `x**n`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_rshift
    >>> gf_rshift([1, 2, 3, 4, 0], 3, ZZ)
    ([1, 2], [3, 4, 0])
    """
    if not n:
        return f, []
    else:
        return f[:-n], f[-n:]

def gf_pow(f, n, p, K):
    """
    Computes `f**n` in `GF(p)[x]` using repeated squaring.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_pow
    >>> gf_pow([3, 2, 4], 3, 5, ZZ)
    [2, 4, 4, 2, 2, 1, 4]
    """
    if not n:
        return [K.one]
    elif n == 1:
        return f
    elif n == 2:
        return gf_sqr(f, p, K)

    h = [K.one]

    while True:
        if n & 1:
            h = gf_mul(h, f, p, K)
            n -= 1

        n >>= 1

        if not n:
            break

        f = gf_sqr(f, p, K)

    return h

def gf_pow_mod(f, n, g, p, K):
    """
    Computes `f**n` in `GF(p)[x]/(g)` using repeated squaring.

    Given polynomials `f` and `g` in `GF(p)[x]` and a non-negative
    integer `n`, efficiently computes `f**n (mod g)` i.e. remainder
    from division `f**n` by `g` using repeated squaring algorithm.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_pow_mod
    >>> gf_pow_mod([3, 2, 4], 3, [1, 1], 5, ZZ)
    []

    References
    ==========

    .. [Gathen99] J. von zur Gathen, J. Gerhard, Modern Computer Algebra,
       First Edition, Cambridge University Press, 1999, pp. 69

    """
    if not n:
        return [K.one]
    elif n == 1:
        return gf_rem(f, g, p, K)
    elif n == 2:
        return gf_rem(gf_sqr(f, p, K), g, p, K)

    h = [K.one]

    while True:
        if n & 1:
            h = gf_mul(h, f, p, K)
            h = gf_rem(h, g, p, K)
            n -= 1

        n >>= 1

        if not n:
            break

        f = gf_sqr(f, p, K)
        f = gf_rem(f, g, p, K)

    return h

def gf_gcd(f, g, p, K):
    """
    Euclidean Algorithm in `GF(p)[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_gcd
    >>> gf_gcd([3, 2, 4], [2, 2, 3], 5, ZZ)
    [1, 3]
    """
    while g:
        f, g = g, gf_rem(f, g, p, K)

    return gf_monic(f, p, K)[1]

def gf_lcm(f, g, p, K):
    """
    Computes polynomial LCM in `GF(p)[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_lcm
    >>> gf_lcm([3, 2, 4], [2, 2, 3], 5, ZZ)
    [1, 2, 0, 4]
    """
    if not f or not g:
        return []

    h = gf_exquo(gf_mul(f, g, p, K),
                 gf_gcd(f, g, p, K), p, K)

    return gf_monic(h, p, K)[1]

def gf_cofactors(f, g, p, K):
    """
    Returns polynomial GCD and cofactors in `GF(p)[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_cofactors
    >>> gf_cofactors([3, 2, 4], [2, 2, 3], 5, ZZ)
    ([1, 3], [3, 3], [2, 1])
    """
    if not f and not g:
        return ([], [], [])

    h = gf_gcd(f, g, p, K)

    return (h, gf_exquo(f, h, p, K),
               gf_exquo(g, h, p, K))

def gf_gcdex(f, g, p, K):
    """
    Extended Euclidean Algorithm in `GF(p)[x]`.

    Given polynomials `f` and `g` in `GF(p)[x]`, computes polynomials
    `s`, `t` and `h`, such that `h = gcd(f, g)` and `s*f + t*g = h`. The
    typical application of EEA is solving polynomial diophantine equations.

    Consider polynomials `f = (x + 7) (x + 1)`, `g = (x + 7) (x**2 + 1)`
    in `GF(11)[x]`. Application of Extended Euclidean Algorithm gives::

       >>> from sympy.polys.galoistools import gf_gcdex, gf_mul, gf_add
       >>> from sympy.polys.algebratools import ZZ

       >>> s, t, g = gf_gcdex([1,8,7], [1,7,1,7], 11, ZZ)
       >>> s, t, g
       ([5, 6], [6], [1, 7])

    As result we obtained polynomials `s = 5*x + 6` and `t = 6`, and
    additionally `gcd(f, g) = x + 7`. This is correct because::

       >>> S = gf_mul(s,   [1,8,7], 11, ZZ)
       >>> T = gf_mul(t, [1,7,1,7], 11, ZZ)

       >>> gf_add(S, T, 11, ZZ) == [1, 7]
       True

    References
    ==========

    .. [Gathen99] J. von zur Gathen, J. Gerhard, Modern Computer Algebra,
       First Edition, Cambridge University Press, 1999, pp. 46

    """
    if not (f or g):
        return [K.one], [], []

    p0, r0 = gf_monic(f, p, K)
    p1, r1 = gf_monic(g, p, K)

    if not f:
        return [], [K.invert(p1, p)], r1
    if not g:
        return [K.invert(p0, p)], [], r0

    s0, s1 = [K.invert(p0, p)], []
    t0, t1 = [], [K.invert(p1, p)]

    while True:
        Q, R = gf_div(r0, r1, p, K)

        if not R:
            break

        (lc, r1), r0 = gf_monic(R, p, K), r1

        inv = K.invert(lc, p)

        s = gf_sub_mul(s0, s1, Q, p, K)
        t = gf_sub_mul(t0, t1, Q, p, K)

        s1, s0 = gf_mul_ground(s, inv, p, K), s1
        t1, t0 = gf_mul_ground(t, inv, p, K), t1

    return s1, t1, r1

def gf_monic(f, p, K):
    """
    Returns LC and a monic polynomial in `GF(p)[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_monic
    >>> gf_monic([3, 2, 4], 5, ZZ)
    (3, [1, 4, 3])
    """
    if not f:
        return K.zero, []
    else:
        lc = f[0]

        if K.is_one(lc):
            return lc, list(f)
        else:
            return lc, gf_exquo_ground(f, lc, p, K)

@cythonized("df,n")
def gf_diff(f, p, K):
    """
    Differentiate polynomial in `GF(p)[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_diff
    >>> gf_diff([3, 2, 4], 5, ZZ)
    [1, 2]
    """
    df = gf_degree(f)

    h, n = [K.zero]*df, df

    for coeff in f[:-1]:
        coeff *= K(n)
        coeff %= p

        if coeff:
            h[df-n] = coeff

        n -= 1

    return gf_strip(h)

def gf_eval(f, a, p, K):
    """
    Evaluate `f(a)` in `GF(p)` using Horner scheme.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_eval
    >>> gf_eval([3, 2, 4], 2, 5, ZZ)
    0
    """
    result = K.zero

    for c in f:
        result *= a
        result += c
        result %= p

    return result

def gf_multi_eval(f, A, p, K):
    """
    Evaluate `f(a)` for `a` in `[a_1, ..., a_n]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_multi_eval
    >>> gf_multi_eval([3, 2, 4], [0, 1, 2, 3, 4], 5, ZZ)
    [4, 4, 0, 2, 0]
    """
    return [ gf_eval(f, a, p, K) for a in A ]

def gf_compose(f, g, p, K):
    """
    Compute polynomial composition `f(g)` in `GF(p)[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_compose
    >>> gf_compose([3, 2, 4], [2, 2, 2], 5, ZZ)
    [2, 4, 0, 3, 0]
    """
    if len(g) <= 1:
        return gf_strip([gf_eval(f, gf_LC(g, K), p, K)])

    if not f:
        return []

    h = [f[0]]

    for c in f[1:]:
        h = gf_mul(h, g, p, K)
        h = gf_add_ground(h, c, p, K)

    return h

def gf_compose_mod(g, h, f, p, K):
    """
    Compute polynomial composition `g(h)` in `GF(p)[x]/(f)`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_compose_mod
    >>> gf_compose_mod([3, 2, 4], [2, 2, 2], [4, 3], 5, ZZ)
    [4]
    """
    if not g:
        return []

    comp = [g[0]]

    for a in g[1:]:
        comp = gf_mul(comp, h, p, K)
        comp = gf_add_ground(comp, a, p, K)
        comp = gf_rem(comp, f, p, K)

    return comp

@cythonized("n")
def gf_trace_map(a, b, c, n, f, p, K):
    """
    Compute polynomial trace map in `GF(p)[x]/(f)`.

    Given a polynomial `f` in `GF(p)[x]`, polynomials `a`, `b`, `c`
    in the quotient ring `GF(p)[x]/(f)` such that `b = c**t (mod f)`
    for some positive power `t` of `p`, and a positive integer `n`,
    returns a mapping::

       a -> a**t**n, a + a**t + a**t**2 + ... + a**t**n (mod f)

    In factorization context, `b = x**p mod f` and `c = x mod f`.
    This way we can efficiently compute trace polynomials in equal
    degree factorization routine, much faster than with other methods,
    like iterated Frobenius algorithm, for large degrees.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_trace_map
    >>> gf_trace_map([1, 2], [4, 4], [1, 1], 4, [3, 2, 4], 5, ZZ)
    ([1, 3], [1, 3])

    References
    ==========

    .. [Gathen92] J. von zur Gathen, V. Shoup, Computing Frobenius Maps
       and Factoring Polynomials, ACM Symposium on Theory of Computing,
       1992, pp. 187-224

    """
    u = gf_compose_mod(a, b, f, p, K)
    v = b

    if n & 1:
        U = gf_add(a, u, p, K)
        V = b
    else:
        U = a
        V = c

    n >>= 1

    while n:
        u = gf_add(u, gf_compose_mod(u, v, f, p, K), p, K)
        v = gf_compose_mod(v, v, f, p, K)

        if n & 1:
            U = gf_add(U, gf_compose_mod(u, V, f, p, K), p, K)
            V = gf_compose_mod(v, V, f, p, K)

        n >>= 1

    return gf_compose_mod(a, V, f, p, K), U

@cythonized("i,n")
def gf_random(n, p, K):
    """
    Generate a random polynomial in `GF(p)[x]` of degree `n`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_random
    >>> gf_random(10, 5, ZZ) #doctest: +SKIP
    [1, 2, 3, 2, 1, 1, 1, 2, 0, 4, 2]
    """
    return [K.one] + [ K(int(uniform(0, p))) for i in xrange(0, n) ]

@cythonized("i,n")
def gf_irreducible(n, p, K):
    """
    Generate random irreducible polynomial of degree `n` in `GF(p)[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_irreducible
    >>> gf_irreducible(10, 5, ZZ) #doctest: +SKIP
    [1, 4, 2, 2, 3, 2, 4, 1, 4, 0, 4]
    """
    while True:
        f = gf_random(n, p, K)

        if gf_irreducible_p(f, p, K):
            return f

@cythonized("i,n")
def gf_irred_p_ben_or(f, p, K):
    """
    Ben-Or's polynomial irreducibility test over finite fields.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_irred_p_ben_or
    >>> gf_irred_p_ben_or([1, 4, 2, 2, 3, 2, 4, 1, 4, 0, 4], 5, ZZ)
    True
    >>> gf_irred_p_ben_or([3, 2, 4], 5, ZZ)
    False
    """
    n = gf_degree(f)

    if n <= 1:
        return True

    _, f = gf_monic(f, p, K)

    H = h = gf_pow_mod([K.one, K.zero], p, f, p, K)

    for i in xrange(0, n//2):
        g = gf_sub(h, [K.one, K.zero], p, K)

        if gf_gcd(f, g, p, K) == [K.one]:
            h = gf_compose_mod(h, H, f, p, K)
        else:
            return False

    return True

@cythonized("i,n,d")
def gf_irred_p_rabin(f, p, K):
    """
    Rabin's polynomial irreducibility test over finite fields.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_irred_p_rabin
    >>> gf_irred_p_rabin([1, 4, 2, 2, 3, 2, 4, 1, 4, 0, 4], 5, ZZ)
    True
    >>> gf_irred_p_rabin([3, 2, 4], 5, ZZ)
    False
    """
    n = gf_degree(f)

    if n <= 1:
        return True

    _, f = gf_monic(f, p, K)

    x = [K.one, K.zero]

    H = h = gf_pow_mod(x, p, f, p, K)

    indices = set([ n//d for d in factorint(n) ])

    for i in xrange(1, n):
        if i in indices:
            g = gf_sub(h, x, p, K)

            if gf_gcd(f, g, p, K) != [K.one]:
                return False

        h = gf_compose_mod(h, H, f, p, K)

    return h == x

_irred_methods = {
    'ben-or' : gf_irred_p_ben_or,
    'rabin'  : gf_irred_p_rabin,
}

def gf_irreducible_p(f, p, K):
    """
    Test irreducibility of a polynomial `f` in `GF(p)[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_irreducible_p
    >>> gf_irreducible_p([1, 4, 2, 2, 3, 2, 4, 1, 4, 0, 4], 5, ZZ)
    True
    >>> gf_irreducible_p([3, 2, 4], 5, ZZ)
    False
    """
    method = query('GF_IRRED_METHOD')

    if method is not None:
        irred = _irred_methods[method](f, p, K)
    else:
        irred = gf_irred_p_rabin(f, p, K)

    return irred

def gf_sqf_p(f, p, K):
    """
    Returns `True` if `f` is square-free in `GF(p)[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_sqf_p
    >>> gf_sqf_p([3, 2, 4], 5, ZZ)
    True
    >>> gf_sqf_p([2, 4, 4, 2, 2, 1, 4], 5, ZZ)
    False
    """
    _, f = gf_monic(f, p, K)

    if not f:
        return True
    else:
        return gf_gcd(f, gf_diff(f, p, K), p, K) == [K.one]

def gf_sqf_part(f, p, K):
    """
    Returns square-free part of a `GF(p)[x]` polynomial.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_sqf_part
    >>> gf_sqf_part([1, 1, 3, 0, 1, 0, 2, 2, 1], 5, ZZ)
    [1, 4, 3]
    """
    _, sqf = gf_sqf_list(f, p, K)

    g = [K.one]

    for f, _ in sqf:
        g = gf_mul(g, f, p, K)

    return g

@cythonized("i,n,d,r")
def gf_sqf_list(f, p, K):
    """
    Returns the square-free decomposition of a `GF(p)[x]` polynomial.

    Given a polynomial `f` in `GF(p)[x]`, returns the leading coefficient
    of `f` and a square-free decomposition `f_1**e_1 f_2**e_2 ... f_k**e_k`
    such that all `f_i` are monic polynomials and `(f_i, f_j)` for `i != j`
    are co-prime and `e_1 ... e_k` are given in increasing order. All
    trivial terms (i.e. `f_i = 1`) aren't included in the output.

    Consider polynomial `f = x**11 + 1` over `GF(11)[x]`::

       >>> from sympy.polys.galoistools import (
       ...     gf_from_dict, gf_diff, gf_sqf_list, gf_pow,
       ... )
       ... # doctest: +NORMALIZE_WHITESPACE

       >>> from sympy.polys.algebratools import ZZ

       >>> f = gf_from_dict({11: 1, 0: 1}, 11, ZZ)

    Note that `f'(x) = 0`::

       >>> gf_diff(f, 11, ZZ)
       []

    This phenomenon doesn't happen in characteristic zero. However we can
    still compute square-free decomposition of `f` using `gf_sqf()`::

       >>> gf_sqf_list(f, 11, ZZ)
       (1, [([1, 1], 11)])

    We obtained factorization `f = (x + 1)**11`. This is correct because::

       >>> gf_pow([1, 1], 11, 11, ZZ) == f
       True

    References
    ==========

    .. [Geddes92] K. Geddes, S. Czapor, G. Labahn, Algorithms for
       Computer Algebra, First Edition, Springer, 1992, pp. 343-347

    """
    n, sqf, factors, r = 1, False, [], int(p)

    lc, f = gf_monic(f, p, K)

    if gf_degree(f) < 1:
        return lc, []

    while True:
        F = gf_diff(f, p, K)

        if F != []:
            g = gf_gcd(f, F, p, K)
            h = gf_exquo(f, g, p, K)

            i = 1

            while h != [K.one]:
                G = gf_gcd(g, h, p, K)
                H = gf_exquo(h, G, p, K)

                if gf_degree(H) > 0:
                    factors.append((H, i*n))

                g, h, i = gf_exquo(g, G, p, K), G, i+1

            if g == [K.one]:
                sqf = True
            else:
                f = g

        if not sqf:
            d = gf_degree(f) // r

            for i in xrange(0, d+1):
                f[i] = f[i*r]

            f, n = f[:d+1], n*r
        else:
            break

    return lc, factors

@cythonized("n,i,j,r")
def gf_Qmatrix(f, p, K):
    """
    Calculate Berlekamp's `Q` matrix.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_Qmatrix
    >>> gf_Qmatrix([3, 2, 4], 5, ZZ)
    [[1, 0],
     [3, 4]]
    >>> gf_Qmatrix([1, 0, 0, 0, 1], 5, ZZ)
    [[1, 0, 0, 0],
     [0, 4, 0, 0],
     [0, 0, 1, 0],
     [0, 0, 0, 4]]
    """
    n, r = gf_degree(f), int(p)

    q = [K.one] + [K.zero]*(n-1)
    Q = [list(q)] + [[]]*(n-1)

    for i in xrange(1, (n-1)*r + 1):
        qq, c = [(-q[-1]*f[-1]) % p], q[-1]

        for j in xrange(1, n):
            qq.append((q[j-1] - c*f[-j-1]) % p)

        if not (i % r):
            Q[i//r] = list(qq)

        q = qq

    return Q

@cythonized("n,i,j,k")
def gf_Qbasis(Q, p, K):
    """
    Compute a basis of the kernel of `Q`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_Qmatrix, gf_Qbasis
    >>> gf_Qbasis(gf_Qmatrix([1, 0, 0, 0, 1], 5, ZZ), 5, ZZ)
    [[1, 0, 0, 0], [0, 0, 1, 0]]
    >>> gf_Qbasis(gf_Qmatrix([3, 2, 4], 5, ZZ), 5, ZZ)
    [[1, 0]]
    """
    Q, n = [ list(q) for q in Q ], len(Q)

    for k in xrange(0, n):
        Q[k][k] = (Q[k][k] - K.one) % p

    for k in xrange(0, n):
        for i in xrange(k, n):
            if Q[k][i]:
                break
        else:
            continue

        inv = K.invert(Q[k][i], p)

        for j in xrange(0, n):
            Q[j][i] = (Q[j][i]*inv) % p

        for j in xrange(0, n):
            t = Q[j][k]
            Q[j][k] = Q[j][i]
            Q[j][i] = t

        for i in xrange(0, n):
            if i != k:
                q = Q[k][i]

                for j in xrange(0, n):
                    Q[j][i] = (Q[j][i] - Q[j][k]*q) % p

    for i in xrange(0, n):
        for j in xrange(0, n):
            if i == j:
                Q[i][j] = (K.one - Q[i][j]) % p
            else:
                Q[i][j] = (-Q[i][j]) % p

    basis = []

    for q in Q:
        if any(q):
            basis.append(q)

    return basis

@cythonized("i,k")
def gf_berlekamp(f, p, K):
    """
    Factor a square-free `f` in `GF(p)[x]` for small `p`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_berlekamp
    >>> gf_berlekamp([1, 0, 0, 0, 1], 5, ZZ)
    [[1, 0, 2], [1, 0, 3]]
    """
    Q = gf_Qmatrix(f, p, K)
    V = gf_Qbasis(Q, p, K)

    for i, v in enumerate(V):
        V[i] = gf_strip(list(reversed(v)))

    factors = [f]

    for k in xrange(1, len(V)):
        for f in list(factors):
            s = K.zero

            while s < p:
                g = gf_sub_ground(V[k], s, p, K)
                h = gf_gcd(f, g, p, K)

                if h != [K.one] and h != f:
                    factors.remove(f)

                    f = gf_exquo(f, h, p, K)
                    factors.extend([f, h])

                if len(factors) == len(V):
                    return _sort_factors(factors, multiple=False)

                s += K.one

    return _sort_factors(factors, multiple=False)

@cythonized("i")
def gf_ddf_zassenhaus(f, p, K):
    """
    Cantor-Zassenhaus: Deterministic Distinct Degree Factorization

    Given a monic square-free polynomial `f` in `GF(p)[x]`, computes
    partial distinct degree factorization `f_1 ... f_d` of `f` where
    `deg(f_i) != deg(f_j)` for `i != j`. The result is returned as a
    list of pairs `(f_i, e_i)` where `deg(f_i) > 0` and `e_i > 0` is
    an argument to the equal degree factorization routine.

    Consider the polynomial `x**15 - 1` in `GF(11)[x]`::

       >>> from sympy.polys.galoistools import gf_from_dict
       >>> from sympy.polys.algebratools import ZZ

       >>> f = gf_from_dict({15: 1, 0: -1}, 11, ZZ)

    Distinct degree factorization gives::

       >>> from sympy.polys.galoistools import gf_ddf_zassenhaus

       >>> gf_ddf_zassenhaus(f, 11, ZZ)
       [([1, 0, 0, 0, 0, 10], 1), ([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1], 2)]

    which means `x**15 - 1 = (x**5 - 1) (x**10 + x**5 + 1)`. To obtain
    factorization into irreducibles, use equal degree factorization
    procedure (EDF) with each of the factors.

    References
    ==========

    .. [Gathen99] J. von zur Gathen, J. Gerhard, Modern Computer Algebra,
       First Edition, Cambridge University Press, 1999, pp. 356

    .. [Geddes92] K. Geddes, S. Czapor, G. Labahn, Algorithms for Computer
       Algebra, First Edition, Springer, 1992, pp. 368-371

    """
    i, g, factors = 1, [K.one, K.zero], []

    while 2*i <= gf_degree(f):
        g = gf_pow_mod(g, int(p), f, p, K)
        h = gf_gcd(f, gf_sub(g, [K.one, K.zero], p, K), p, K)

        if h != [K.one]:
            factors.append((h, i))

            f = gf_exquo(f, h, p, K)
            g = gf_rem(g, f, p, K)

        i += 1

    if f != [K.one]:
        return factors + [(f, gf_degree(f))]
    else:
        return factors

@cythonized("n,N,i")
def gf_edf_zassenhaus(f, n, p, K):
    """
    Cantor-Zassenhaus: Probabilistic Equal Degree Factorization

    Given a monic square-free polynomial `f` in `GF(p)[x]` and integer
    `n` such that `n` divides `deg(f)`, returns all irreducible factors
    `f_1 ... f_d` of `f`, each of degree `n`. This is a complete
    factorization in Galois fields.

    Consider the square-free polynomial `f = x**3 + x**2 + x + 1` in
    `GF(5)[x]`. Let's compute its irreducible factors of degree one::

       >>> from sympy.polys.galoistools import gf_edf_zassenhaus
       >>> from sympy.polys.algebratools import ZZ

       >>> gf_edf_zassenhaus([1,1,1,1], 1, 5, ZZ)
       [[1, 1], [1, 2], [1, 3]]

    References
    ==========

    .. [Gathen99] J. von zur Gathen, J. Gerhard, Modern Computer Algebra,
       First Edition, Cambridge University Press, 1999, pp. 358

    .. [Geddes92] K. Geddes, S. Czapor, G. Labahn, Algorithms for Computer
       Algebra, First Edition, Springer, 1992, pp. 371-373

    """
    factors, q = [f], int(p)

    if gf_degree(f) <= n:
        return factors

    N = gf_degree(f) // n

    while len(factors) < N:
        r = gf_random(2*n-1, p, K)

        if p == 2:
            h = r

            for i in xrange(0, 2**(n*N-1)):
                r = gf_pow_mod(r, 2, f, p, K)
                h = gf_add(h, r, p, K)

            g = gf_gcd(f, h, p, K)
        else:
            h = gf_pow_mod(r, (q**n-1) // 2, f, p, K)
            g = gf_gcd(f, gf_sub_ground(h, K.one, p, K), p, K)

        if g != [K.one] and g != f:
            factors = gf_edf_zassenhaus(g, n, p, K) \
                    + gf_edf_zassenhaus(gf_exquo(f, g, p, K), n, p, K)

    return _sort_factors(factors, multiple=False)

@cythonized("n,k,i,j")
def gf_ddf_shoup(f, p, K):
    """
    Kaltofen-Shoup: Deterministic Distinct Degree Factorization

    Given a monic square-free polynomial `f` in `GF(p)[x]`, computes
    partial distinct degree factorization `f_1 ... f_d` of `f` where
    `deg(f_i) != deg(f_j)` for `i != j`. The result is returned as a
    list of pairs `(f_i, e_i)` where `deg(f_i) > 0` and `e_i > 0` is
    an argument to the equal degree factorization routine.

    This algorithm is an improved version of Zassenhaus algorithm for
    large `deg(f)` and modulus `p` (especially for `deg(f) ~ lg(p)`).

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_ddf_shoup
    >>> f = [1, 235, 2571, 1793, 2162, 2398, 1171, 2174, 1259, 2057,
    ... 2571, 2117, 154]
    >>> g = [([1, 2837, 2277], 1), ([1, 315, 2158, 2657, 243, 1288,
    ... 118, 1672, 515, 1960, 1832], 10)]

    >>> gf_ddf_shoup(f, 2917, ZZ) == g
    True

    References
    ==========

    .. [Kaltofen98] E. Kaltofen, V. Shoup, Subquadratic-time Factoring
       of Polynomials over Finite Fields, Mathematics of Computation,
       Volume 67, Issue 223, 1998, pp. 1179-1197

    .. [Shoup95] V. Shoup, A New Polynomial Factorization Algorithm and
       its Implementation, Journal of Symbolic Computation, Volume 20,
       Issue 4, 1995, pp. 363-397

    .. [Gathen92] J. von zur Gathen, V. Shoup, Computing Frobenius Maps
       and Factoring Polynomials, ACM Symposium on Theory of Computing,
       1992, pp. 187-224

    """
    n = gf_degree(f)
    k = int(ceil(sqrt(n//2)))

    h = gf_pow_mod([K.one, K.zero], int(p), f, p, K)

    U = [[K.one,K.zero], h] + [K.zero]*(k-1)

    for i in xrange(2, k+1):
        U[i] = gf_compose_mod(U[i-1], h, f, p, K)

    h, U = U[k], U[:k]
    V = [h] + [K.zero]*(k-1)

    for i in xrange(1, k):
        V[i] = gf_compose_mod(V[i-1], h, f, p, K)

    factors = []

    for i, v in enumerate(V):
        h, j = [K.one], k-1

        for u in U:
            g = gf_sub(v, u, p, K)
            h = gf_mul(h, g, p, K)
            h = gf_rem(h, f, p, K)

        g = gf_gcd(f, h, p, K)
        f = gf_exquo(f, g, p, K)

        for u in reversed(U):
            h = gf_sub(v, u, p, K)
            F = gf_gcd(g, h, p, K)

            if F != [K.one]:
                factors.append((F, k*(i+1)-j))

            g, j = gf_exquo(g, F, p, K), j-1

    if f != [K.one]:
        factors.append((f, gf_degree(f)))

    return factors

@cythonized("n,N,q")
def gf_edf_shoup(f, n, p, K):
    """
    Gathen-Shoup: Probabilistic Equal Degree Factorization

    Given a monic square-free polynomial `f` in `GF(p)[x]` and integer
    `n` such that `n` divides `deg(f)`, returns all irreducible factors
    `f_1 ... f_d` of `f`, each of degree `n`. This is a complete
    factorization over Galois fields.

    This algorithm is an improved version of Zassenhaus algorithm for
    large `deg(f)` and modulus `p` (especially for `deg(f) ~ lg(p)`).

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_edf_shoup
    >>> gf_edf_shoup([1, 2837, 2277], 1, 2917, ZZ)
    [[1, 852], [1, 1985]]

    References
    ==========

    .. [Shoup91] V. Shoup, A Fast Deterministic Algorithm for Factoring
       Polynomials over Finite Fields of Small Characteristic, In
       Proceedings of International Symposium on Symbolic and
       Algebraic Computation, 1991, pp. 14-21

    .. [Gathen92] J. von zur Gathen, V. Shoup, Computing Frobenius Maps
       and Factoring Polynomials, ACM Symposium on Theory of Computing,
       1992, pp. 187-224

    """
    N, q = gf_degree(f), int(p)

    if not N:
        return []
    if N <= n:
        return [f]

    factors, x = [f], [K.one, K.zero]

    r = gf_random(N-1, p, K)

    h = gf_pow_mod(x, q, f, p, K)
    H = gf_trace_map(r, h, x, n-1, f, p, K)[1]

    if p == 2:
        h1 = gf_gcd(f, H, p, K)
        h2 = gf_exquo(f, h1, p, K)

        factors = gf_edf_shoup(h1, n, p, K) \
                + gf_edf_shoup(h2, n, p, K)
    else:
        h = gf_pow_mod(H, (q-1)//2, f, p, K)

        h1 = gf_gcd(f, h, p, K)
        h2 = gf_gcd(f, gf_sub_ground(h, K.one, p, K), p, K)
        h3 = gf_exquo(f, gf_mul(h1, h2, p, K), p, K)

        factors = gf_edf_shoup(h1, n, p, K) \
                + gf_edf_shoup(h2, n, p, K) \
                + gf_edf_shoup(h3, n, p, K)

    return _sort_factors(factors, multiple=False)

@cythonized("n")
def gf_zassenhaus(f, p, K):
    """
    Factor a square-free `f` in `GF(p)[x]` for medium `p`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_zassenhaus
    >>> gf_zassenhaus([3, 2, 4], 5, ZZ)
    [[3], [1, 1], [1, 3]]
    """
    factors = []

    for factor, n in gf_ddf_zassenhaus(f, p, K):
        factors += gf_edf_zassenhaus(factor, n, p, K)

    return _sort_factors(factors, multiple=False)

@cythonized("n")
def gf_shoup(f, p, K):
    """
    Factor a square-free `f` in `GF(p)[x]` for large `p`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_shoup
    >>> f = [1, 235, 2571, 1793, 2162, 2398, 1171, 2174, 1259, 2057,
    ... 2571, 2117, 154]
    >>> g = [[1, 852], [1, 1985], [1, 315, 2158, 2657, 243, 1288, 118,
    ... 1672, 515, 1960, 1832]]

    >>> gf_shoup(f, 2917, ZZ) == g
    True
    """
    factors = []

    for factor, n in gf_ddf_shoup(f, p, K):
        factors += gf_edf_shoup(factor, n, p, K)

    return _sort_factors(factors, multiple=False)

_factor_methods = {
    'berlekamp'  : gf_berlekamp,  # `p` : small
    'zassenhaus' : gf_zassenhaus, # `p` : medium
    'shoup'      : gf_shoup,      # `p` : large
}

def gf_factor_sqf(f, p, K):
    """
    Factor a square-free polynomial `f` in `GF(p)[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.galoistools import gf_factor_sqf
    >>> gf_factor_sqf([3, 2, 4], 5, ZZ)
    (3, [[1, 1], [1, 3]])
    """
    lc, f = gf_monic(f, p, K)

    if gf_degree(f) < 1:
        return lc, []

    method = query('GF_FACTOR_METHOD')

    if method is not None:
        factors = _factor_methods[method](f, p, K)
    else:
        factors = gf_zassenhaus(f, p, K)

    return lc, factors

@cythonized("n")
def gf_factor(f, p, K):
    """
    Factor (non square-free) polynomials in `GF(p)[x]`.

    Given a possibly non square-free polynomial `f` in `GF(p)[x]`, returns
    its complete factorization into irreducibles::

                 f_1(x)**e_1 f_2(x)**e_2 ... f_d(x)**e_d

    where each `f_i` is a monic polynomial and `gcd(f_i, f_j) == 1`, for
    `i != j`.  The result is given as a tuple consisting of the leading
    coefficient of `f` and a list of factors with their multiplicities.

    The algorithm proceeds by first computing square-free decomposition
    of `f` and then iteratively factoring each of the square-free factors.

    Consider a non square-free polynomial `f = (7*x + 1) (x + 2)**2` in
    `GF(11)[x]`. We obtain its factorization into irreducibles as follows::

       >>> from sympy.polys.galoistools import gf_factor
       >>> from sympy.polys.algebratools import ZZ

       >>> gf_factor([5, 2, 7, 2], 11, ZZ)
       (5, [([1, 2], 1), ([1, 8], 2)])

    We arrived with factorization `f = 5 (x + 2) (x + 8)**2`. We didn't
    recover the exact form of the input polynomial because we requested to
    get monic factors of `f` and its leading coefficient separately.

    Square-free factors of `f` can be factored into irreducibles over
    `GF(p)` using three very different methods:

    1. Berlekamp - efficient for very small values of `p` (usually `p < 25`)
    2. Cantor-Zassenhaus - efficient on average input and with "typical" `p`
    3. Shoup-Kaltofen-Gathen - efficient with very large inputs and modulus

    If you want to use a specific factorization method, instead of the default
    one, set GF_FACTOR_METHOD with one of `berlekamp`, `zassenhaus` or `shoup`
    values.

    References
    ==========

    .. [Gathen99] J. von zur Gathen, J. Gerhard, Modern Computer Algebra,
       First Edition, Cambridge University Press, 1999, pp. 365

    """
    lc, f = gf_monic(f, p, K)

    if gf_degree(f) < 1:
        return lc, []

    factors = []

    for g, n in gf_sqf_list(f, p, K)[1]:
        for h in gf_factor_sqf(g, p, K)[1]:
            factors.append((h, n))

    return lc, _sort_factors(factors)


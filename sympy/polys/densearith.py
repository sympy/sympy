"""Arithmetics for dense recursive polynomials in ``K[x]`` or ``K[X]``. """

from __future__ import annotations

from sympy.polys.domains.domain import Domain, Er, Ef, Eeuclid, Eabs, Eordered

from sympy.polys.densebasic import (
    dup, dmp, _dup, _dmp, _dmp_ground,
    dup_slice, dup_truncate,
    dup_reverse,
    dup_LC, dmp_LC,
    dup_degree, dmp_degree,
    dup_strip, dmp_strip,
    dmp_zero_p, dmp_zero,
    dmp_one_p, dmp_one,
    dmp_ground, dmp_zeros)
from sympy.polys.polyerrors import (ExactQuotientFailed, PolynomialDivisionFailed)


def dup_add_term(f: dup[Er], c: Er, i: int, K: Domain[Er]) -> dup[Er]:
    """
    Add ``c*x**i`` to ``f`` in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_add_term(x**2 - 1, ZZ(2), 4)
    2*x**4 + x**2 - 1

    """
    if not c:
        return f

    n = len(f)
    m = n - i - 1

    if i == n - 1:
        return dup_strip([f[0] + c] + f[1:], K)
    else:
        if i >= n:
            return [c] + [K.zero]*(i - n) + f
        else:
            return f[:m] + [f[m] + c] + f[m + 1:]


def dmp_add_term(f: dmp[Er], c: dmp[Er], i: int, u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Add ``c(x_2..x_u)*x_0**i`` to ``f`` in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dmp_add_term(x*y + 1, 2, 2)
    2*x**2 + x*y + 1

    """
    if not u:
        return _dmp(dup_add_term(_dup(f), _dmp_ground(c), i, K))

    v = u - 1

    if dmp_zero_p(c, v):
        return f

    n = len(f)
    m = n - i - 1

    if i == n - 1:
        return dmp_strip([dmp_add(f[0], c, v, K)] + f[1:], u, K)
    else:
        if i >= n:
            return [c] + dmp_zeros(i - n, v, K) + f
        else:
            return f[:m] + [dmp_add(f[m], c, v, K)] + f[m + 1:]


def dup_sub_term(f: dup[Er], c: Er, i: int, K: Domain[Er]) -> dup[Er]:
    """
    Subtract ``c*x**i`` from ``f`` in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_sub_term(2*x**4 + x**2 - 1, ZZ(2), 4)
    x**2 - 1

    """
    if not c:
        return f

    n = len(f)
    m = n - i - 1

    if i == n - 1:
        return dup_strip([f[0] - c] + f[1:], K)
    else:
        if i >= n:
            return [-c] + [K.zero]*(i - n) + f
        else:
            return f[:m] + [f[m] - c] + f[m + 1:]


def dmp_sub_term(f: dmp[Er], c: dmp[Er], i: int, u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Subtract ``c(x_2..x_u)*x_0**i`` from ``f`` in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dmp_sub_term(2*x**2 + x*y + 1, 2, 2)
    x*y + 1

    """
    if not u:
        return _dmp(dup_sub_term(_dup(f), _dmp_ground(c), i, K))

    v = u - 1

    if dmp_zero_p(c, v):
        return f

    n = len(f)
    m = n - i - 1

    if i == n - 1:
        return dmp_strip([dmp_sub(f[0], c, v, K)] + f[1:], u, K)
    else:
        if i >= n:
            return [dmp_neg(c, v, K)] + dmp_zeros(i - n, v, K) + f
        else:
            return f[:m] + [dmp_sub(f[m], c, v, K)] + f[m + 1:]


def dup_mul_term(f: dup[Er], c: Er, i: int, K: Domain[Er]) -> dup[Er]:
    """
    Multiply ``f`` by ``c*x**i`` in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_mul_term(x**2 - 1, ZZ(3), 2)
    3*x**4 - 3*x**2

    """
    if not c or not f:
        return []
    else:
        return [ cf * c for cf in f ] + [K.zero]*i


def dmp_mul_term(f: dmp[Er], c: dmp[Er], i: int, u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Multiply ``f`` by ``c(x_2..x_u)*x_0**i`` in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dmp_mul_term(x**2*y + x, 3*y, 2)
    3*x**4*y**2 + 3*x**3*y

    """
    if not u:
        return _dmp(dup_mul_term(_dup(f), _dmp_ground(c), i, K))

    v = u - 1

    if dmp_zero_p(f, u):
        return f
    if dmp_zero_p(c, v):
        return dmp_zero(u, K)
    else:
        return [ dmp_mul(cf, c, v, K) for cf in f ] + dmp_zeros(i, v, K)


def dup_add_ground(f: dup[Er], c: Er, K: Domain[Er]) -> dup[Er]:
    """
    Add an element of the ground domain to ``f``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_add_ground(x**3 + 2*x**2 + 3*x + 4, ZZ(4))
    x**3 + 2*x**2 + 3*x + 8

    """
    return dup_add_term(f, c, 0, K)


def dmp_add_ground(f: dmp[Er], c: Er, u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Add an element of the ground domain to ``f``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dmp_add_ground(x**3 + 2*x**2 + 3*x + 4, ZZ(4))
    x**3 + 2*x**2 + 3*x + 8

    """
    return dmp_add_term(f, dmp_ground(c, u - 1, K), 0, u, K)


def dup_sub_ground(f: dup[Er], c: Er, K: Domain[Er]) -> dup[Er]:
    """
    Subtract an element of the ground domain from ``f``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_sub_ground(x**3 + 2*x**2 + 3*x + 4, ZZ(4))
    x**3 + 2*x**2 + 3*x

    """
    return dup_sub_term(f, c, 0, K)


def dmp_sub_ground(f: dmp[Er], c: Er, u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Subtract an element of the ground domain from ``f``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dmp_sub_ground(x**3 + 2*x**2 + 3*x + 4, ZZ(4))
    x**3 + 2*x**2 + 3*x

    """
    return dmp_sub_term(f, dmp_ground(c, u - 1, K), 0, u, K)


def dup_mul_ground(f: dup[Er], c: Er, K: Domain[Er]) -> dup[Er]:
    """
    Multiply ``f`` by a constant value in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_mul_ground(x**2 + 2*x - 1, ZZ(3))
    3*x**2 + 6*x - 3

    """
    if not c or not f:
        return []
    else:
        return [ cf * c for cf in f ]


def dmp_mul_ground(f: dmp[Er], c: Er, u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Multiply ``f`` by a constant value in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dmp_mul_ground(2*x + 2*y, ZZ(3))
    6*x + 6*y

    """
    if not u:
        return _dmp(dup_mul_ground(_dup(f), c, K))

    v = u - 1

    return [ dmp_mul_ground(cf, c, v, K) for cf in f ]


def dup_quo_ground(f: dup[Er], c: Er, K: Domain[Er]) -> dup[Er]:
    """
    Quotient by a constant in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ, QQ

    >>> R, x = ring("x", ZZ)
    >>> R.dup_quo_ground(3*x**2 + 2, ZZ(2))
    x**2 + 1

    >>> R, x = ring("x", QQ)
    >>> R.dup_quo_ground(3*x**2 + 2, QQ(2))
    3/2*x**2 + 1

    """
    if not c:
        raise ZeroDivisionError('polynomial division')
    if not f:
        return f

    if K.is_Field:
        return [ K.quo(cf, c) for cf in f ]
    else:
        return [ cf // c for cf in f ] # type: ignore


def dmp_quo_ground(
    f: dmp[Er], c: Er, u: int, K: Domain[Er]
) -> dmp[Er]:
    """
    Quotient by a constant in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ, QQ

    >>> R, x,y = ring("x,y", ZZ)
    >>> R.dmp_quo_ground(2*x**2*y + 3*x, ZZ(2))
    x**2*y + x

    >>> R, x,y = ring("x,y", QQ)
    >>> R.dmp_quo_ground(2*x**2*y + 3*x, QQ(2))
    x**2*y + 3/2*x

    """
    if not u:
        return _dmp(dup_quo_ground(_dup(f), c, K))

    v = u - 1

    return [ dmp_quo_ground(cf, c, v, K) for cf in f ]


def dup_exquo_ground(f: dup[Er], c: Er, K: Domain[Er]) -> dup[Er]:
    """
    Exact quotient by a constant in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, QQ
    >>> R, x = ring("x", QQ)

    >>> R.dup_exquo_ground(x**2 + 2, QQ(2))
    1/2*x**2 + 1

    """
    if not c:
        raise ZeroDivisionError('polynomial division')
    if not f:
        return f

    return [ K.exquo(cf, c) for cf in f ]


def dmp_exquo_ground(f: dmp[Er], c: Er, u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Exact quotient by a constant in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, QQ
    >>> R, x,y = ring("x,y", QQ)

    >>> R.dmp_exquo_ground(x**2*y + 2*x, QQ(2))
    1/2*x**2*y + x

    """
    if not u:
        return _dmp(dup_exquo_ground(_dup(f), c, K))

    v = u - 1

    return [ dmp_exquo_ground(cf, c, v, K) for cf in f ]


def dup_lshift(f: dup[Er], n: int, K: Domain[Er]) -> dup[Er]:
    """
    Efficiently multiply ``f`` by ``x**n`` in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_lshift(x**2 + 1, 2)
    x**4 + x**2

    """
    if not f:
        return f
    else:
        return f + [K.zero]*n


def dup_rshift(f: dup[Er], n: int, K: Domain[Er]) -> dup[Er]:
    """
    Efficiently divide ``f`` by ``x**n`` in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_rshift(x**4 + x**2, 2)
    x**2 + 1
    >>> R.dup_rshift(x**4 + x**2 + 2, 2)
    x**2 + 1

    """
    return f[:-n]


def dup_abs(f: dup[Eabs], K: Domain[Eabs]) -> dup[Eabs]:
    """
    Make all coefficients positive in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_abs(x**2 - 1)
    x**2 + 1

    """
    return [ K.abs(coeff) for coeff in f ]


def dmp_abs(f: dmp[Eabs], u: int, K: Domain[Eabs]) -> dmp[Eabs]:
    """
    Make all coefficients positive in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dmp_abs(x**2*y - x)
    x**2*y + x

    """
    if not u:
        return _dmp(dup_abs(_dup(f), K))

    v = u - 1

    return [ dmp_abs(cf, v, K) for cf in f ]


def dup_neg(f: dup[Er], K: Domain[Er]) -> dup[Er]:
    """
    Negate a polynomial in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_neg(x**2 - 1)
    -x**2 + 1

    """
    return [ -coeff for coeff in f ]


def dmp_neg(f: dmp[Er], u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Negate a polynomial in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dmp_neg(x**2*y - x)
    -x**2*y + x

    """
    if not u:
        return _dmp(dup_neg(_dup(f), K))

    v = u - 1

    return [ dmp_neg(cf, v, K) for cf in f ]


def dup_add(f: dup[Er], g: dup[Er], K: Domain[Er]) -> dup[Er]:
    """
    Add dense polynomials in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_add(x**2 - 1, x - 2)
    x**2 + x - 3

    """
    if not f:
        return g
    if not g:
        return f

    df = dup_degree(f)
    dg = dup_degree(g)

    if df == dg:
        return dup_strip([ a + b for a, b in zip(f, g) ], K)
    else:
        k = abs(df - dg)

        if df > dg:
            h, f = f[:k], f[k:]
        else:
            h, g = g[:k], g[k:]

        return h + [ a + b for a, b in zip(f, g) ]


def dmp_add(f: dmp[Er], g: dmp[Er], u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Add dense polynomials in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dmp_add(x**2 + y, x**2*y + x)
    x**2*y + x**2 + x + y

    """
    if not u:
        return _dmp(dup_add(_dup(f), _dup(g), K))

    df = dmp_degree(f, u)

    if df < 0:
        return g

    dg = dmp_degree(g, u)

    if dg < 0:
        return f

    v = u - 1

    if df == dg:
        return dmp_strip([ dmp_add(a, b, v, K) for a, b in zip(f, g) ], u, K)
    else:
        k = abs(df - dg)

        if df > dg:
            h, f = f[:k], f[k:]
        else:
            h, g = g[:k], g[k:]

        return h + [ dmp_add(a, b, v, K) for a, b in zip(f, g) ]


def dup_sub(f: dup[Er], g: dup[Er], K: Domain[Er]) -> dup[Er]:
    """
    Subtract dense polynomials in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_sub(x**2 - 1, x - 2)
    x**2 - x + 1

    """
    if not f:
        return dup_neg(g, K)
    if not g:
        return f

    df = dup_degree(f)
    dg = dup_degree(g)

    if df == dg:
        return dup_strip([ a - b for a, b in zip(f, g) ], K)
    else:
        k = abs(df - dg)

        if df > dg:
            h, f = f[:k], f[k:]
        else:
            h, g = dup_neg(g[:k], K), g[k:]

        return h + [ a - b for a, b in zip(f, g) ]


def dmp_sub(f: dmp[Er], g: dmp[Er], u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Subtract dense polynomials in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dmp_sub(x**2 + y, x**2*y + x)
    -x**2*y + x**2 - x + y

    """
    if not u:
        return _dmp(dup_sub(_dup(f), _dup(g), K))

    df = dmp_degree(f, u)

    if df < 0:
        return dmp_neg(g, u, K)

    dg = dmp_degree(g, u)

    if dg < 0:
        return f

    v = u - 1

    if df == dg:
        return dmp_strip([ dmp_sub(a, b, v, K) for a, b in zip(f, g) ], u, K)
    else:
        k = abs(df - dg)

        if df > dg:
            h, f = f[:k], f[k:]
        else:
            h, g = dmp_neg(g[:k], u, K), g[k:]

        return h + [ dmp_sub(a, b, v, K) for a, b in zip(f, g) ]


def dup_add_mul(f: dup[Er], g: dup[Er], h: dup[Er], K: Domain[Er]) -> dup[Er]:
    """
    Returns ``f + g*h`` where ``f, g, h`` are in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_add_mul(x**2 - 1, x - 2, x + 2)
    2*x**2 - 5

    """
    return dup_add(f, dup_mul(g, h, K), K)


def dmp_add_mul(f: dmp[Er], g: dmp[Er], h: dmp[Er], u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Returns ``f + g*h`` where ``f, g, h`` are in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dmp_add_mul(x**2 + y, x, x + 2)
    2*x**2 + 2*x + y

    """
    return dmp_add(f, dmp_mul(g, h, u, K), u, K)


def dup_sub_mul(f: dup[Er], g: dup[Er], h: dup[Er], K: Domain[Er]) -> dup[Er]:
    """
    Returns ``f - g*h`` where ``f, g, h`` are in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_sub_mul(x**2 - 1, x - 2, x + 2)
    3

    """
    return dup_sub(f, dup_mul(g, h, K), K)


def dmp_sub_mul(f: dmp[Er], g: dmp[Er], h: dmp[Er], u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Returns ``f - g*h`` where ``f, g, h`` are in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dmp_sub_mul(x**2 + y, x, x + 2)
    -2*x + y

    """
    return dmp_sub(f, dmp_mul(g, h, u, K), u, K)


def dup_mul(f: dup[Er], g: dup[Er], K: Domain[Er]) -> dup[Er]:
    """
    Multiply dense polynomials in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_mul(x - 2, x + 2)
    x**2 - 4

    """
    if f == g:
        return dup_sqr(f, K)

    if not (f and g):
        return []

    df = dup_degree(f)
    dg = dup_degree(g)

    n = max(df, dg) + 1

    if n < 100 or not K.is_Exact:
        h: list[Er] = []

        for i in range(0, df + dg + 1):
            coeff = K.zero

            for j in range(max(0, i - dg), min(df, i) + 1):
                coeff += f[j]*g[i - j]

            h.append(coeff)

        return dup_strip(h, K)
    else:
        # Use Karatsuba's algorithm (divide and conquer), see e.g.:
        # Joris van der Hoeven, Relax But Don't Be Too Lazy,
        # J. Symbolic Computation, 11 (2002), section 3.1.1.
        n2 = n//2

        fl, gl = dup_slice(f, 0, n2, K), dup_slice(g, 0, n2, K)

        fh = dup_rshift(dup_slice(f, n2, n, K), n2, K)
        gh = dup_rshift(dup_slice(g, n2, n, K), n2, K)

        lo, hi = dup_mul(fl, gl, K), dup_mul(fh, gh, K)

        mid = dup_mul(dup_add(fl, fh, K), dup_add(gl, gh, K), K)
        mid = dup_sub(mid, dup_add(lo, hi, K), K)

        return dup_add(dup_add(lo, dup_lshift(mid, n2, K), K),
                       dup_lshift(hi, 2*n2, K), K)


def dmp_mul(f: dmp[Er], g: dmp[Er], u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Multiply dense polynomials in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dmp_mul(x*y + 1, x)
    x**2*y + x

    """
    if not u:
        return _dmp(dup_mul(_dup(f), _dup(g), K))

    if f == g:
        return dmp_sqr(f, u, K)

    df = dmp_degree(f, u)

    if df < 0:
        return f

    dg = dmp_degree(g, u)

    if dg < 0:
        return g

    h: list[dmp[Er]] = []
    v = u - 1

    for i in range(0, df + dg + 1):
        coeff = dmp_zero(v, K)

        for j in range(max(0, i - dg), min(df, i) + 1):
            coeff = dmp_add(coeff, dmp_mul(f[j], g[i - j], v, K), v, K)

        h.append(coeff)

    return dmp_strip(h, u, K)


def _dup_series_mul_base(f: dup[Er], g: dup[Er], n: int, K: Domain[Er]) -> dup[Er]:
    """
    Helper function for dup_series_mul that uses a naive O(n^2) algorithm
    """
    f = f[::-1]
    g = g[::-1]

    h = [K.zero] * n
    for i in range(min(len(f), n)):
        for j in range(min(len(g), n - i)):
            h[i + j] += f[i] * g[j]

    return dup_reverse(h, K)


def _dup_series_mul_karatsuba(f: dup[Er], g: dup[Er], prec: int, K: Domain[Er]) -> dup[Er]:
    """
    Helper function for dup_series_mul that uses Karatsuba's algorithm
    """
    df = dup_degree(f)
    dg = dup_degree(g)

    n = max(df, dg) + 1
    n2 =  n // 2

    fl = dup_slice(f, 0, n2, K)
    gl = dup_slice(g, 0, n2, K)

    fh = dup_rshift(dup_slice(f, n2, n, K), n2, K)
    gh = dup_rshift(dup_slice(g, n2, n, K), n2, K)

    lo = dup_series_mul(fl, gl, prec, K)

    hi = []
    mid = []
    if n2 < prec:
        hi = dup_series_mul(fh, gh, prec - n2, K)

        a = dup_add(fl, fh, K)
        b = dup_add(gl, gh, K)

        mid = dup_series_mul(a, b, prec - n2, K)
        mid = dup_sub(mid, dup_add(lo, hi, K), K)

    res = lo
    if mid:
        res = dup_add(res, dup_lshift(mid, n2, K), K)
    if hi:
        res = dup_add(res, dup_lshift(hi, 2 * n2, K), K)

    res = dup_truncate(res, prec, K)
    return dup_strip(res, K)


def dup_series_mul(f: dup[Er], g: dup[Er], n: int, K: Domain[Er]) -> dup[Er]:
    """
    Multiply dense polynomials in ``K[[x]]`` modulo ``x**n``.

    Examples
    ========
    >>> from sympy import ZZ
    >>> from sympy.polys.densearith import dup_series_mul
    >>> from sympy.polys.densebasic import dup_from_list, dup_print
    >>> f = dup_from_list([1, 0, 1], ZZ)
    >>> g = dup_from_list([2, 1], ZZ)
    >>> h = dup_series_mul(f, g, 4, ZZ)
    >>> dup_print(h, 'x')
    2*x**3 + x**2 + 2*x + 1

    """
    if not f or not g or n <= 0:
        return []

    if min(len(f), len(g)) < 100:
        return _dup_series_mul_base(f, g, n, K)
    else:
        return _dup_series_mul_karatsuba(f, g, n, K)


def dup_sqr(f: dup[Er], K: Domain[Er]) -> dup[Er]:
    """
    Square dense polynomials in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_sqr(x**2 + 1)
    x**4 + 2*x**2 + 1

    """
    h: list[Er] = []
    df = len(f) - 1

    for i in range(0, 2*df + 1):
        c = K.zero

        jmin = max(0, i - df)
        jmax = min(i, df)

        n = jmax - jmin + 1

        jmax = jmin + n // 2 - 1

        for j in range(jmin, jmax + 1):
            c += f[j]*f[i - j]

        c += c

        if n & 1:
            elem = f[jmax + 1]
            c += elem**2

        h.append(c)

    return dup_strip(h, K)


def dmp_sqr(f: dmp[Er], u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Square dense polynomials in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dmp_sqr(x**2 + x*y + y**2)
    x**4 + 2*x**3*y + 3*x**2*y**2 + 2*x*y**3 + y**4

    """
    if not u:
        return _dmp(dup_sqr(_dup(f), K))

    df = dmp_degree(f, u)

    if df < 0:
        return f

    h: list[dmp[Er]] = []
    v = u - 1

    for i in range(0, 2*df + 1):
        c = dmp_zero(v, K)

        jmin = max(0, i - df)
        jmax = min(i, df)

        n = jmax - jmin + 1

        jmax = jmin + n // 2 - 1

        for j in range(jmin, jmax + 1):
            c = dmp_add(c, dmp_mul(f[j], f[i - j], v, K), v, K)

        c = dmp_mul_ground(c, K(2), v, K)

        if n & 1:
            elem = dmp_sqr(f[jmax + 1], v, K)
            c = dmp_add(c, elem, v, K)

        h.append(c)

    return dmp_strip(h, u, K)


def dup_series_sqr(f: dup[Er], n: int, K: Domain[Er]) -> dup[Er]:
    """
    Square dense polynomials in ``K[[x]]`` modulo ``x**n``.

    Examples
    ========
    >>> from sympy import ZZ
    >>> from sympy.polys.densearith import dup_series_sqr
    >>> from sympy.polys.densebasic import dup_from_list, dup_print
    >>> f = dup_from_list([1, 2, 3], ZZ)
    >>> p = dup_series_sqr(f, 5, ZZ)
    >>> dup_print(p, 'x')
    x**4 + 4*x**3 + 10*x**2 + 12*x + 9

    """
    if not f or n<=0:
        return []

    df = dup_degree(f)

    if df > 100:
        # Use Karatsuba's algorithm for larger polynomials
        return dup_series_mul(f, f, n, K)


    h: list[Er] = []

    for i in range(n):
        c = K.zero

        j_min = max(0, i - df)
        j_max = (i + 1) // 2

        for j in range(j_min, j_max):
            c += f[df - j] * f[df - (i - j)]

        c += c

        if i % 2 == 0:
            j = i // 2
            if j <= df:
                c += f[df - j]**2

        h.append(c)

    return dup_strip(h[::-1], K)


def dup_pow(f: dup[Er], n: int, K: Domain[Er]) -> dup[Er]:
    """
    Raise ``f`` to the ``n``-th power in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_pow(x - 2, 3)
    x**3 - 6*x**2 + 12*x - 8

    """
    if not n:
        return [K.one]
    if n < 0:
        raise ValueError("Cannot raise polynomial to a negative power")
    if n == 1 or not f or f == [K.one]:
        return f

    g = [K.one]

    while True:
        n, m = n//2, n

        if m % 2:
            g = dup_mul(g, f, K)

            if not n:
                break

        f = dup_sqr(f, K)

    return g


def dmp_pow(f: dmp[Er], n: int, u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Raise ``f`` to the ``n``-th power in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dmp_pow(x*y + 1, 3)
    x**3*y**3 + 3*x**2*y**2 + 3*x*y + 1

    """
    if not u:
        return _dmp(dup_pow(_dup(f), n, K))

    if not n:
        return dmp_one(u, K)
    if n < 0:
        raise ValueError("Cannot raise polynomial to a negative power")
    if n == 1 or dmp_zero_p(f, u) or dmp_one_p(f, u, K):
        return f

    g = dmp_one(u, K)

    while True:
        n, m = n//2, n

        if m & 1:
            g = dmp_mul(g, f, u, K)

            if not n:
                break

        f = dmp_sqr(f, u, K)

    return g


def _dup_recurse_pow(f: dup[Er], exp: int, prec: int, K: Domain[Er]) -> dup[Er]:
    """
    Helper function for dup_series_pow.
    """
    if exp == 1:
        return dup_truncate(f, prec, K)

    q, r = divmod(exp, 2)
    half = _dup_recurse_pow(f, q, prec, K)
    square = dup_series_mul(half, half, prec, K)

    if r == 0:
        return square
    return dup_series_mul(square, f, prec, K)


def dup_series_pow(f: dup[Er], n: int, prec: int, K: Domain[Er]) -> dup[Er]:
    """
    Raise polynomial ``f`` to power ``n`` modulo ``x**prec`` in ``K[[x]]``.

    Examples
    ========
    >>> from sympy import ZZ
    >>> from sympy.polys.densearith import dup_series_pow
    >>> from sympy.polys.densebasic import dup_from_list, dup_print
    >>> f = dup_from_list([1, 2, 3], ZZ)
    >>> p = dup_series_pow(f, 2, 5, ZZ)
    >>> dup_print(p, 'x')
    x**4 + 4*x**3 + 10*x**2 + 12*x + 9

    """
    if not f:
        return []
    if n < 0:
        raise ValueError("Negative exponent not supported")
    if n == 0:
        return [K.one]
    if n == 1 or f == [K.one]:
        return dup_truncate(f, prec, K)

    return _dup_recurse_pow(f, n, prec, K)


def dup_pdiv(f: dup[Er], g: dup[Er], K: Domain[Er]) -> tuple[dup[Er], dup[Er]]:
    """
    Polynomial pseudo-division in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_pdiv(x**2 + 1, 2*x - 4)
    (2*x + 4, 20)

    """
    df = dup_degree(f)
    dg = dup_degree(g)

    q: dup[Er] = []
    r = f
    dr = df

    if not g:
        raise ZeroDivisionError("polynomial division")
    elif df < dg:
        return q, r

    N = df - dg + 1
    lc_g = dup_LC(g, K)

    while True:
        lc_r = dup_LC(r, K)
        j, N = dr - dg, N - 1

        Q = dup_mul_ground(q, lc_g, K)
        q = dup_add_term(Q, lc_r, j, K)

        R = dup_mul_ground(r, lc_g, K)
        G = dup_mul_term(g, lc_r, j, K)
        r = dup_sub(R, G, K)

        _dr, dr = dr, dup_degree(r)

        if dr < dg:
            break
        elif not (dr < _dr):
            raise PolynomialDivisionFailed(f, g, K)

    c = lc_g**N

    q = dup_mul_ground(q, c, K)
    r = dup_mul_ground(r, c, K)

    return q, r


def dup_prem(f: dup[Er], g: dup[Er], K: Domain[Er]) -> dup[Er]:
    """
    Polynomial pseudo-remainder in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_prem(x**2 + 1, 2*x - 4)
    20

    """
    df = dup_degree(f)
    dg = dup_degree(g)

    r, dr = f, df

    if not g:
        raise ZeroDivisionError("polynomial division")
    elif df < dg:
        return r

    N = df - dg + 1
    lc_g = dup_LC(g, K)

    while True:
        lc_r = dup_LC(r, K)
        j, N = dr - dg, N - 1

        R = dup_mul_ground(r, lc_g, K)
        G = dup_mul_term(g, lc_r, j, K)
        r = dup_sub(R, G, K)

        _dr, dr = dr, dup_degree(r)

        if dr < dg:
            break
        elif not (dr < _dr):
            raise PolynomialDivisionFailed(f, g, K)

    return dup_mul_ground(r, lc_g**N, K)


def dup_pquo(f: dup[Er], g: dup[Er], K: Domain[Er]) -> dup[Er]:
    """
    Polynomial exact pseudo-quotient in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_pquo(x**2 - 1, 2*x - 2)
    2*x + 2

    >>> R.dup_pquo(x**2 + 1, 2*x - 4)
    2*x + 4

    """
    return dup_pdiv(f, g, K)[0]


def dup_pexquo(f: dup[Er], g: dup[Er], K: Domain[Er]) -> dup[Er]:
    """
    Polynomial pseudo-quotient in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_pexquo(x**2 - 1, 2*x - 2)
    2*x + 2

    >>> R.dup_pexquo(x**2 + 1, 2*x - 4)
    Traceback (most recent call last):
    ...
    ExactQuotientFailed: [2, -4] does not divide [1, 0, 1]

    """
    q, r = dup_pdiv(f, g, K)

    if not r:
        return q
    else:
        raise ExactQuotientFailed(f, g)


def dmp_pdiv(f: dmp[Er], g: dmp[Er], u: int, K: Domain[Er]) -> tuple[dmp[Er], dmp[Er]]:
    """
    Polynomial pseudo-division in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dmp_pdiv(x**2 + x*y, 2*x + 2)
    (2*x + 2*y - 2, -4*y + 4)

    """
    if not u:
        uq, ur = dup_pdiv(_dup(f), _dup(g), K)
        return _dmp(uq), _dmp(ur)

    df = dmp_degree(f, u)
    dg = dmp_degree(g, u)

    if dg < 0:
        raise ZeroDivisionError("polynomial division")

    q, r, dr = dmp_zero(u, K), f, df

    if df < dg:
        return q, r

    N = df - dg + 1
    lc_g = dmp_LC(g, K)

    while True:
        lc_r = dmp_LC(r, K)
        j, N = dr - dg, N - 1

        Q = dmp_mul_term(q, lc_g, 0, u, K)
        q = dmp_add_term(Q, lc_r, j, u, K)

        R = dmp_mul_term(r, lc_g, 0, u, K)
        G = dmp_mul_term(g, lc_r, j, u, K)
        r = dmp_sub(R, G, u, K)

        _dr, dr = dr, dmp_degree(r, u)

        if dr < dg:
            break
        elif not (dr < _dr):
            raise PolynomialDivisionFailed(f, g, K)

    c = dmp_pow(lc_g, N, u - 1, K)

    q = dmp_mul_term(q, c, 0, u, K)
    r = dmp_mul_term(r, c, 0, u, K)

    return q, r


def dmp_prem(f: dmp[Er], g: dmp[Er], u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Polynomial pseudo-remainder in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dmp_prem(x**2 + x*y, 2*x + 2)
    -4*y + 4

    """
    if not u:
        return _dmp(dup_prem(_dup(f), _dup(g), K))

    df = dmp_degree(f, u)
    dg = dmp_degree(g, u)

    if dg < 0:
        raise ZeroDivisionError("polynomial division")

    r, dr = f, df

    if df < dg:
        return r

    N = df - dg + 1
    lc_g = dmp_LC(g, K)

    while True:
        lc_r = dmp_LC(r, K)
        j, N = dr - dg, N - 1

        R = dmp_mul_term(r, lc_g, 0, u, K)
        G = dmp_mul_term(g, lc_r, j, u, K)
        r = dmp_sub(R, G, u, K)

        _dr, dr = dr, dmp_degree(r, u)

        if dr < dg:
            break
        elif not (dr < _dr):
            raise PolynomialDivisionFailed(f, g, K)

    c = dmp_pow(lc_g, N, u - 1, K)

    return dmp_mul_term(r, c, 0, u, K)


def dmp_pquo(f: dmp[Er], g: dmp[Er], u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Polynomial exact pseudo-quotient in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> f = x**2 + x*y
    >>> g = 2*x + 2*y
    >>> h = 2*x + 2

    >>> R.dmp_pquo(f, g)
    2*x

    >>> R.dmp_pquo(f, h)
    2*x + 2*y - 2

    """
    return dmp_pdiv(f, g, u, K)[0]


def dmp_pexquo(f: dmp[Er], g: dmp[Er], u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Polynomial pseudo-quotient in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> f = x**2 + x*y
    >>> g = 2*x + 2*y
    >>> h = 2*x + 2

    >>> R.dmp_pexquo(f, g)
    2*x

    >>> R.dmp_pexquo(f, h)
    Traceback (most recent call last):
    ...
    ExactQuotientFailed: [[2], [2]] does not divide [[1], [1, 0], []]

    """
    q, r = dmp_pdiv(f, g, u, K)

    if dmp_zero_p(r, u):
        return q
    else:
        raise ExactQuotientFailed(f, g)


def dup_rr_div(
    f: dup[Eeuclid], g: dup[Eeuclid], K: Domain[Eeuclid]
) -> tuple[dup[Eeuclid], dup[Eeuclid]]:
    """
    Univariate division with remainder over a ring.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_rr_div(x**2 + 1, 2*x - 4)
    (0, x**2 + 1)

    """
    df = dup_degree(f)
    dg = dup_degree(g)

    q: dup[Eeuclid] = []
    r = f
    dr = df

    if not g:
        raise ZeroDivisionError("polynomial division")
    elif df < dg:
        return q, r

    lc_g = dup_LC(g, K)

    while True:
        lc_r = dup_LC(r, K)

        if lc_r % lc_g:
            break

        c = K.exquo(lc_r, lc_g)
        j = dr - dg

        q = dup_add_term(q, c, j, K)
        h = dup_mul_term(g, c, j, K)
        r = dup_sub(r, h, K)

        _dr, dr = dr, dup_degree(r)

        if dr < dg:
            break
        elif not (dr < _dr):
            raise PolynomialDivisionFailed(f, g, K)

    return q, r


def dmp_rr_div(
        f: dmp[Eeuclid], g: dmp[Eeuclid], u: int, K: Domain[Eeuclid]
) -> tuple[dmp[Eeuclid], dmp[Eeuclid]]:
    """
    Multivariate division with remainder over a ring.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dmp_rr_div(x**2 + x*y, 2*x + 2)
    (0, x**2 + x*y)

    """
    if not u:
        uq, ur = dup_rr_div(_dup(f), _dup(g), K)
        return _dmp(uq), _dmp(ur)

    df = dmp_degree(f, u)
    dg = dmp_degree(g, u)

    if dg < 0:
        raise ZeroDivisionError("polynomial division")

    q, r, dr = dmp_zero(u, K), f, df

    if df < dg:
        return q, r

    lc_g, v = dmp_LC(g, K), u - 1

    while True:
        lc_r = dmp_LC(r, K)
        c, R = dmp_rr_div(lc_r, lc_g, v, K)

        if not dmp_zero_p(R, v):
            break

        j = dr - dg

        q = dmp_add_term(q, c, j, u, K)
        h = dmp_mul_term(g, c, j, u, K)
        r = dmp_sub(r, h, u, K)

        _dr, dr = dr, dmp_degree(r, u)

        if dr < dg:
            break
        elif not (dr < _dr):
            raise PolynomialDivisionFailed(f, g, K)

    return q, r


def dup_ff_div(f: dup[Ef], g: dup[Ef], K: Domain[Ef]) -> tuple[dup[Ef], dup[Ef]]:
    """
    Polynomial division with remainder over a field.

    Examples
    ========

    >>> from sympy.polys import ring, QQ
    >>> R, x = ring("x", QQ)

    >>> R.dup_ff_div(x**2 + 1, 2*x - 4)
    (1/2*x + 1, 5)

    """
    df = dup_degree(f)
    dg = dup_degree(g)

    q: dup[Ef] = []
    r = f
    dr = df

    if not g:
        raise ZeroDivisionError("polynomial division")
    elif df < dg:
        return q, r

    lc_g = dup_LC(g, K)

    while True:
        lc_r = dup_LC(r, K)

        c = K.exquo(lc_r, lc_g)
        j = dr - dg

        q = dup_add_term(q, c, j, K)
        h = dup_mul_term(g, c, j, K)
        r = dup_sub(r, h, K)

        _dr, dr = dr, dup_degree(r)

        if dr < dg:
            break
        elif dr == _dr and not K.is_Exact:
            # remove leading term created by rounding error
            r = dup_strip(r[1:], K)
            dr = dup_degree(r)
            if dr < dg:
                break
        elif not (dr < _dr):
            raise PolynomialDivisionFailed(f, g, K)

    return q, r


def dmp_ff_div(
    f: dmp[Ef], g: dmp[Ef], u: int, K: Domain[Ef]
) -> tuple[dmp[Ef], dmp[Ef]]:
    """
    Polynomial division with remainder over a field.

    Examples
    ========

    >>> from sympy.polys import ring, QQ
    >>> R, x,y = ring("x,y", QQ)

    >>> R.dmp_ff_div(x**2 + x*y, 2*x + 2)
    (1/2*x + 1/2*y - 1/2, -y + 1)

    """
    if not u:
        uq, ur = dup_ff_div(_dup(f), _dup(g), K)
        return _dmp(uq), _dmp(ur)

    df = dmp_degree(f, u)
    dg = dmp_degree(g, u)

    if dg < 0:
        raise ZeroDivisionError("polynomial division")

    q, r, dr = dmp_zero(u, K), f, df

    if df < dg:
        return q, r

    lc_g, v = dmp_LC(g, K), u - 1

    while True:
        lc_r = dmp_LC(r, K)
        c, R = dmp_ff_div(lc_r, lc_g, v, K)

        if not dmp_zero_p(R, v):
            break

        j = dr - dg

        q = dmp_add_term(q, c, j, u, K)
        h = dmp_mul_term(g, c, j, u, K)
        r = dmp_sub(r, h, u, K)

        _dr, dr = dr, dmp_degree(r, u)

        if dr < dg:
            break
        elif not (dr < _dr):
            raise PolynomialDivisionFailed(f, g, K)

    return q, r


def dup_div(f: dup[Er], g: dup[Er], K: Domain[Er]) -> tuple[dup[Er], dup[Er]]:
    """
    Polynomial division with remainder in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ, QQ

    >>> R, x = ring("x", ZZ)
    >>> R.dup_div(x**2 + 1, 2*x - 4)
    (0, x**2 + 1)

    >>> R, x = ring("x", QQ)
    >>> R.dup_div(x**2 + 1, 2*x - 4)
    (1/2*x + 1, 5)

    """
    if K.is_Field:
        return dup_ff_div(f, g, K) # type: ignore
    else:
        return dup_rr_div(f, g, K) # type: ignore


def dup_rem(f: dup[Er], g: dup[Er], K: Domain[Er]) -> dup[Er]:
    """
    Returns polynomial remainder in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ, QQ

    >>> R, x = ring("x", ZZ)
    >>> R.dup_rem(x**2 + 1, 2*x - 4)
    x**2 + 1

    >>> R, x = ring("x", QQ)
    >>> R.dup_rem(x**2 + 1, 2*x - 4)
    5

    """
    return dup_div(f, g, K)[1]


def dup_quo(f: dup[Er], g: dup[Er], K: Domain[Er]) -> dup[Er]:
    """
    Returns polynomial quotient in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ, QQ

    >>> R, x = ring("x", ZZ)
    >>> R.dup_quo(x**2 + 1, 2*x - 4)
    0

    >>> R, x = ring("x", QQ)
    >>> R.dup_quo(x**2 + 1, 2*x - 4)
    1/2*x + 1

    """
    return dup_div(f, g, K)[0]


def dup_exquo(f: dup[Er], g: dup[Er], K: Domain[Er]) -> dup[Er]:
    """
    Returns exact polynomial quotient in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_exquo(x**2 - 1, x - 1)
    x + 1

    >>> R.dup_exquo(x**2 + 1, 2*x - 4)
    Traceback (most recent call last):
    ...
    ExactQuotientFailed: [2, -4] does not divide [1, 0, 1]

    """
    q, r = dup_div(f, g, K)

    if not r:
        return q
    else:
        raise ExactQuotientFailed(f, g)


def dmp_div(
    f: dmp[Er], g: dmp[Er], u: int, K: Domain[Er]
) -> tuple[dmp[Er], dmp[Er]]:
    """
    Polynomial division with remainder in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ, QQ

    >>> R, x,y = ring("x,y", ZZ)
    >>> R.dmp_div(x**2 + x*y, 2*x + 2)
    (0, x**2 + x*y)

    >>> R, x,y = ring("x,y", QQ)
    >>> R.dmp_div(x**2 + x*y, 2*x + 2)
    (1/2*x + 1/2*y - 1/2, -y + 1)

    """
    if K.is_Field:
        return dmp_ff_div(f, g, u, K) # type: ignore
    else:
        return dmp_rr_div(f, g, u, K) # type: ignore


def dmp_rem(f: dmp[Er], g: dmp[Er], u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Returns polynomial remainder in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ, QQ

    >>> R, x,y = ring("x,y", ZZ)
    >>> R.dmp_rem(x**2 + x*y, 2*x + 2)
    x**2 + x*y

    >>> R, x,y = ring("x,y", QQ)
    >>> R.dmp_rem(x**2 + x*y, 2*x + 2)
    -y + 1

    """
    return dmp_div(f, g, u, K)[1]


def dmp_quo(f: dmp[Er], g: dmp[Er], u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Returns exact polynomial quotient in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ, QQ

    >>> R, x,y = ring("x,y", ZZ)
    >>> R.dmp_quo(x**2 + x*y, 2*x + 2)
    0

    >>> R, x,y = ring("x,y", QQ)
    >>> R.dmp_quo(x**2 + x*y, 2*x + 2)
    1/2*x + 1/2*y - 1/2

    """
    return dmp_div(f, g, u, K)[0]


def dmp_exquo(f: dmp[Er], g: dmp[Er], u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Returns polynomial quotient in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> f = x**2 + x*y
    >>> g = x + y
    >>> h = 2*x + 2

    >>> R.dmp_exquo(f, g)
    x

    >>> R.dmp_exquo(f, h)
    Traceback (most recent call last):
    ...
    ExactQuotientFailed: [[2], [2]] does not divide [[1], [1, 0], []]

    """
    q, r = dmp_div(f, g, u, K)

    if dmp_zero_p(r, u):
        return q
    else:
        raise ExactQuotientFailed(f, g)


def dup_max_norm(f: dup[Eordered], K: Domain[Eordered]) -> Eordered:
    """
    Returns maximum norm of a polynomial in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_max_norm(-x**2 + 2*x - 3)
    3

    """
    if not f:
        return K.zero
    else:
        return max(dup_abs(f, K))


def dmp_max_norm(f: dmp[Eordered], u: int, K: Domain[Eordered]) -> Eordered:
    """
    Returns maximum norm of a polynomial in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dmp_max_norm(2*x*y - x - 3)
    3

    """
    if not u:
        return dup_max_norm(_dup(f), K)

    v = u - 1

    return max(dmp_max_norm(c, v, K) for c in f)


def dup_l1_norm(f: dup[Eabs], K: Domain[Eabs]) -> Eabs:
    """
    Returns l1 norm of a polynomial in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_l1_norm(2*x**3 - 3*x**2 + 1)
    6

    """
    if not f:
        return K.zero
    else:
        return K.sum(dup_abs(f, K))


def dmp_l1_norm(f: dmp[Eabs], u: int, K: Domain[Eabs]) -> Eabs:
    """
    Returns l1 norm of a polynomial in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dmp_l1_norm(2*x*y - x - 3)
    6

    """
    if not u:
        return dup_l1_norm(_dup(f), K)

    v = u - 1

    return K.sum(dmp_l1_norm(c, v, K) for c in f)


def dup_l2_norm_squared(f: dup[Er], K: Domain[Er]) -> Er:
    """
    Returns squared l2 norm of a polynomial in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_l2_norm_squared(2*x**3 - 3*x**2 + 1)
    14

    """
    return K.sum([coeff**2 for coeff in f])


def dmp_l2_norm_squared(f: dmp[Er], u: int, K: Domain[Er]) -> Er:
    """
    Returns squared l2 norm of a polynomial in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dmp_l2_norm_squared(2*x*y - x - 3)
    14

    """
    if not u:
        return dup_l2_norm_squared(_dup(f), K)

    v = u - 1

    return K.sum(dmp_l2_norm_squared(c, v, K) for c in f)


def dup_expand(polys: list[dup[Er]], K: Domain[Er]) -> dup[Er]:
    """
    Multiply together several polynomials in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_expand([x**2 - 1, x, 2])
    2*x**3 - 2*x

    """
    if not polys:
        return [K.one]

    f = polys[0]

    for g in polys[1:]:
        f = dup_mul(f, g, K)

    return f


def dmp_expand(polys: list[dmp[Er]], u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Multiply together several polynomials in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dmp_expand([x**2 + y**2, x + 1])
    x**3 + x**2 + x*y**2 + y**2

    """
    if not polys:
        return dmp_one(u, K)

    f = polys[0]

    for g in polys[1:]:
        f = dmp_mul(f, g, u, K)

    return f

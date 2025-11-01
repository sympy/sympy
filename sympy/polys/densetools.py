"""Advanced tools for dense recursive polynomials in ``K[x]`` or ``K[X]``. """

from __future__ import annotations

from sympy.polys.domains.domain import Domain, Er, Eg, Ef, Eeuclid
from sympy.polys.domains.field import Field
from sympy.polys.densearith import (
    dup_add_term, dmp_add_term,
    dup_lshift, dup_rshift,
    dup_add, dmp_add,
    dup_sub, dmp_sub,
    dup_mul, dmp_mul, dup_series_mul,
    dup_sqr,
    dup_div,
    dup_series_pow,
    dup_rem, dmp_rem,
    dup_mul_ground, dmp_mul_ground,
    dup_quo_ground, dmp_quo_ground,
    dup_exquo_ground, dmp_exquo_ground,
)
from sympy.polys.densebasic import (
    dup, dmp, _dup, _dmp, _dmp2, _ground_dmp,
    dup_strip, dmp_strip, dup_truncate,
    dup_convert, dmp_convert,
    dup_degree, dmp_degree,
    dmp_to_dict,
    dmp_from_dict,
    dup_LC, dmp_LC, dmp_ground_LC,
    dup_TC, dmp_TC,
    dmp_zero, dmp_ground,
    dmp_zero_p,
    dup_to_raw_dict, dup_from_raw_dict,
    dmp_to_raw_dict,
    dmp_zeros,
    dmp_include,
    dup_nth,
)
from sympy.polys.polyerrors import (
    MultivariatePolynomialError,
    DomainError
)

from math import ceil as _ceil, log2 as _log2, sqrt

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from sympy.external.gmpy import MPQ
    from sympy.polys.domains.ringextension import RingExtension
    from sympy.polys.domains.algebraicfield import AlgebraicField, Alg


def dup_integrate(f: dup[Ef], m: int, K: Field[Ef]) -> dup[Ef]:
    """
    Computes the indefinite integral of ``f`` in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, QQ
    >>> R, x = ring("x", QQ)

    >>> R.dup_integrate(x**2 + 2*x, 1)
    1/3*x**3 + x**2
    >>> R.dup_integrate(x**2 + 2*x, 2)
    1/12*x**4 + 1/3*x**3

    """
    if m <= 0 or not f:
        return f

    g = [K.zero]*m

    for i, c in enumerate(reversed(f)):
        n = i + 1

        for j in range(1, m):
            n *= i + j + 1

        g.insert(0, K.exquo(c, K(n)))

    return g


def dmp_integrate(f: dmp[Ef], m: int, u: int, K: Field[Ef]) -> dmp[Ef]:
    """
    Computes the indefinite integral of ``f`` in ``x_0`` in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, QQ
    >>> R, x,y = ring("x,y", QQ)

    >>> R.dmp_integrate(x + 2*y, 1)
    1/2*x**2 + 2*x*y
    >>> R.dmp_integrate(x + 2*y, 2)
    1/6*x**3 + x**2*y

    """
    if not u:
        return _dmp(dup_integrate(_dup(f), m, K))

    if m <= 0 or dmp_zero_p(f, u):
        return f

    g, v = dmp_zeros(m, u - 1, K), u - 1

    for i, c in enumerate(reversed(f)):
        n = i + 1

        for j in range(1, m):
            n *= i + j + 1

        g.insert(0, dmp_quo_ground(c, K(n), v, K))

    return g


def _rec_integrate_in(g: dmp[Ef], m: int, v: int, i: int, j: int, K: Field[Ef]) -> dmp[Ef]:
    """Recursive helper for :func:`dmp_integrate_in`."""
    if i == j:
        return dmp_integrate(g, m, v, K)

    w, i = v - 1, i + 1

    return dmp_strip([ _rec_integrate_in(c, m, w, i, j, K) for c in g ], v, K)


def dmp_integrate_in(f: dmp[Ef], m: int, j: int, u: int, K: Field[Ef]) -> dmp[Ef]:
    """
    Computes the indefinite integral of ``f`` in ``x_j`` in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, QQ
    >>> R, x,y = ring("x,y", QQ)

    >>> R.dmp_integrate_in(x + 2*y, 1, 0)
    1/2*x**2 + 2*x*y
    >>> R.dmp_integrate_in(x + 2*y, 1, 1)
    x*y + y**2

    """
    if j < 0 or j > u:
        raise IndexError("0 <= j <= u expected, got u = %d, j = %d" % (u, j))

    return _rec_integrate_in(f, m, u, 0, j, K)


def dup_diff(f: dup[Er], m: int, K: Domain[Er]) -> dup[Er]:
    """
    ``m``-th order derivative of a polynomial in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_diff(x**3 + 2*x**2 + 3*x + 4, 1)
    3*x**2 + 4*x + 3
    >>> R.dup_diff(x**3 + 2*x**2 + 3*x + 4, 2)
    6*x + 4

    """
    if m <= 0:
        return f

    n = dup_degree(f)

    if n < m:
        return []

    deriv: list[Er] = []

    if m == 1:
        for coeff in f[:-m]:
            deriv.append(K(n)*coeff)
            n -= 1
    else:
        for coeff in f[:-m]:
            k = n

            for i in range(n - 1, n - m, -1):
                k *= i

            deriv.append(K(k)*coeff)
            n -= 1

    return dup_strip(deriv, K)


def dmp_diff(f: dmp[Er], m: int, u: int, K: Domain[Er]) -> dmp[Er]:
    """
    ``m``-th order derivative in ``x_0`` of a polynomial in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> f = x*y**2 + 2*x*y + 3*x + 2*y**2 + 3*y + 1

    >>> R.dmp_diff(f, 1)
    y**2 + 2*y + 3
    >>> R.dmp_diff(f, 2)
    0

    """
    if not u:
        return _dmp(dup_diff(_dup(f), m, K))
    if m <= 0:
        return f

    n = dmp_degree(f, u)

    if n < m:
        return dmp_zero(u, K)

    deriv: list[dmp[Er]] = []
    v = u - 1

    if m == 1:
        for coeff in f[:-m]:
            deriv.append(dmp_mul_ground(coeff, K(n), v, K))
            n -= 1
    else:
        for coeff in f[:-m]:
            k = n

            for i in range(n - 1, n - m, -1):
                k *= i

            deriv.append(dmp_mul_ground(coeff, K(k), v, K))
            n -= 1

    return dmp_strip(deriv, u, K)


def _rec_diff_in(g: dmp[Er], m: int, v: int, i: int, j: int, K: Domain[Er]) -> dmp[Er]:
    """Recursive helper for :func:`dmp_diff_in`."""
    if i == j:
        return dmp_diff(g, m, v, K)

    w, i = v - 1, i + 1

    return dmp_strip([ _rec_diff_in(c, m, w, i, j, K) for c in g ], v, K)


def dmp_diff_in(f: dmp[Er], m: int, j: int, u: int, K: Domain[Er]) -> dmp[Er]:
    """
    ``m``-th order derivative in ``x_j`` of a polynomial in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> f = x*y**2 + 2*x*y + 3*x + 2*y**2 + 3*y + 1

    >>> R.dmp_diff_in(f, 1, 0)
    y**2 + 2*y + 3
    >>> R.dmp_diff_in(f, 1, 1)
    2*x*y + 2*x + 4*y + 3

    """
    if j < 0 or j > u:
        raise IndexError("0 <= j <= %s expected, got %s" % (u, j))

    return _rec_diff_in(f, m, u, 0, j, K)


def dup_eval(f: dup[Er], a: Er, K: Domain[Er]) -> Er:
    """
    Evaluate a polynomial at ``x = a`` in ``K[x]`` using Horner scheme.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_eval(x**2 + 2*x + 3, 2)
    11

    """
    if not a:
        return K.convert(dup_TC(f, K))

    result = K.zero

    for c in f:
        result *= a
        result += c

    return result


def dmp_eval(f: dmp[Er], a: Er, u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Evaluate a polynomial at ``x_0 = a`` in ``K[X]`` using the Horner scheme.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dmp_eval(2*x*y + 3*x + y + 2, 2)
    5*y + 8

    """
    if not u:
        return _ground_dmp(dup_eval(_dup(f), a, K))

    if not a:
        return dmp_TC(f, K)

    result, v = dmp_LC(f, K), u - 1

    for coeff in f[1:]:
        result = dmp_mul_ground(result, a, v, K)
        result = dmp_add(result, coeff, v, K)

    return result


def _rec_eval_in(g: dmp[Er], a: Er, v: int, i: int, j: int, K: Domain[Er]) -> dmp[Er]:
    """Recursive helper for :func:`dmp_eval_in`."""
    if i == j:
        return dmp_eval(g, a, v, K)

    v, i = v - 1, i + 1

    return dmp_strip([ _rec_eval_in(c, a, v, i, j, K) for c in g ], v, K)


def dmp_eval_in(f: dmp[Er], a: Er, j: int, u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Evaluate a polynomial at ``x_j = a`` in ``K[X]`` using the Horner scheme.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> f = 2*x*y + 3*x + y + 2

    >>> R.dmp_eval_in(f, 2, 0)
    5*y + 8
    >>> R.dmp_eval_in(f, 2, 1)
    7*x + 4

    """
    if j < 0 or j > u:
        raise IndexError("0 <= j <= %s expected, got %s" % (u, j))

    return _rec_eval_in(f, a, u, 0, j, K)


def _rec_eval_tail(g: dmp[Er], i: int, A: list[Er], u: int, K: Domain[Er]) -> dmp[Er]:
    """Recursive helper for :func:`dmp_eval_tail`."""
    if i == u:
        return _ground_dmp(dup_eval(_dup(g), A[-1], K))
    else:
        h = [ _rec_eval_tail(c, i + 1, A, u, K) for c in g ]

        if i < u - len(A) + 1:
            return h
        else:
            return _ground_dmp(dup_eval(_dup(h), A[-u + i - 1], K))


def dmp_eval_tail(f: dmp[Er], A: list[Er], u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Evaluate a polynomial at ``x_j = a_j, ...`` in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> f = 2*x*y + 3*x + y + 2

    >>> R.dmp_eval_tail(f, [2])
    7*x + 4
    >>> R.dmp_eval_tail(f, [2, 2])
    18

    """
    if not A:
        return f

    if dmp_zero_p(f, u):
        return dmp_zero(u - len(A), K)

    e = _rec_eval_tail(f, 0, A, u, K)

    if u == len(A) - 1:
        return e
    else:
        return dmp_strip(e, u - len(A), K)


def _rec_diff_eval(g: dmp[Er], m: int, a: Er, v: int, i: int, j: int, K: Domain[Er]) -> dmp[Er]:
    """Recursive helper for :func:`dmp_diff_eval`."""
    if i == j:
        return dmp_eval(dmp_diff(g, m, v, K), a, v, K)

    v, i = v - 1, i + 1

    return dmp_strip([ _rec_diff_eval(c, m, a, v, i, j, K) for c in g ], v, K)


def dmp_diff_eval_in(f: dmp[Er], m: int, a: Er, j: int, u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Differentiate and evaluate a polynomial in ``x_j`` at ``a`` in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> f = x*y**2 + 2*x*y + 3*x + 2*y**2 + 3*y + 1

    >>> R.dmp_diff_eval_in(f, 1, 2, 0)
    y**2 + 2*y + 3
    >>> R.dmp_diff_eval_in(f, 1, 2, 1)
    6*x + 11

    """
    if j > u:
        raise IndexError("-%s <= j < %s expected, got %s" % (u, u, j))
    if not j:
        return dmp_eval(dmp_diff(f, m, u, K), a, u, K)

    return _rec_diff_eval(f, m, a, u, 0, j, K)


def dup_trunc(f: dup[Eeuclid], p: Eeuclid, K: Domain[Eeuclid]) -> dup[Eeuclid]:
    """
    Reduce a ``K[x]`` polynomial modulo a constant ``p`` in ``K``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_trunc(2*x**3 + 3*x**2 + 5*x + 7, ZZ(3))
    -x**3 - x + 1

    """
    if K.is_ZZ:
        g: list[Eeuclid] = []

        for c in f:
            c = c % p

            if c > p // 2: # type: ignore
                g.append(c - p)
            else:
                g.append(c)
    elif K.is_FiniteField:
        # XXX: python-flint's nmod does not support %
        pi = int(p) # type: ignore
        g = [ K(int(c) % pi) for c in f ] # type: ignore
    else:
        g = [ c % p for c in f ]

    return dup_strip(g, K)


def dmp_trunc(f: dmp[Er], p: dmp[Er], u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Reduce a ``K[X]`` polynomial modulo a polynomial ``p`` in ``K[Y]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> f = 3*x**2*y + 8*x**2 + 5*x*y + 6*x + 2*y + 3
    >>> g = (y - 1).drop(x)

    >>> R.dmp_trunc(f, g)
    11*x**2 + 11*x + 5

    """
    return dmp_strip([ dmp_rem(c, p, u - 1, K) for c in f ], u, K)


def dmp_ground_trunc(f: dmp[Eeuclid], p: Eeuclid, u: int, K: Domain[Eeuclid]) -> dmp[Eeuclid]:
    """
    Reduce a ``K[X]`` polynomial modulo a constant ``p`` in ``K``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> f = 3*x**2*y + 8*x**2 + 5*x*y + 6*x + 2*y + 3

    >>> R.dmp_ground_trunc(f, ZZ(3))
    -x**2 - x*y - y

    """
    if not u:
        return _dmp(dup_trunc(_dup(f), p, K))

    v = u - 1

    return dmp_strip([ dmp_ground_trunc(c, p, v, K) for c in f ], u, K)


def dup_monic(f: dup[Ef], K: Field[Ef]) -> dup[Ef]:
    """
    Divide all coefficients by ``LC(f)`` in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ, QQ

    >>> R, x = ring("x", ZZ)
    >>> R.dup_monic(3*x**2 + 6*x + 9)
    x**2 + 2*x + 3

    >>> R, x = ring("x", QQ)
    >>> R.dup_monic(3*x**2 + 4*x + 2)
    x**2 + 4/3*x + 2/3

    """
    if not f:
        return f

    lc = dup_LC(f, K)

    if K.is_one(lc):
        return f
    else:
        return dup_exquo_ground(f, lc, K)


def dmp_ground_monic(f: dmp[Ef], u: int, K: Field[Ef]) -> dmp[Ef]:
    """
    Divide all coefficients by ``LC(f)`` in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ, QQ

    >>> R, x,y = ring("x,y", ZZ)
    >>> f = 3*x**2*y + 6*x**2 + 3*x*y + 9*y + 3

    >>> R.dmp_ground_monic(f)
    x**2*y + 2*x**2 + x*y + 3*y + 1

    >>> R, x,y = ring("x,y", QQ)
    >>> f = 3*x**2*y + 8*x**2 + 5*x*y + 6*x + 2*y + 3

    >>> R.dmp_ground_monic(f)
    x**2*y + 8/3*x**2 + 5/3*x*y + 2*x + 2/3*y + 1

    """
    if not u:
        return _dmp(dup_monic(_dup(f), K))

    if dmp_zero_p(f, u):
        return f

    lc = dmp_ground_LC(f, u, K)

    if K.is_one(lc):
        return f
    else:
        return dmp_exquo_ground(f, lc, u, K)


def dup_content(f: dup[Er], K: Domain[Er]) -> Er:
    """
    Compute the GCD of coefficients of ``f`` in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ, QQ

    >>> R, x = ring("x", ZZ)
    >>> f = 6*x**2 + 8*x + 12

    >>> R.dup_content(f)
    2

    >>> R, x = ring("x", QQ)
    >>> f = 6*x**2 + 8*x + 12

    >>> R.dup_content(f)
    2

    """
    from sympy.polys.domains import QQ

    if not f:
        return K.zero

    cont = K.zero

    if K == QQ:
        for c in f:
            cont = K.gcd(cont, c)
    else:
        for c in f:
            cont = K.gcd(cont, c)

            if K.is_one(cont):
                break

    return cont


def dmp_ground_content(f: dmp[Er], u: int, K: Domain[Er]) -> Er:
    """
    Compute the GCD of coefficients of ``f`` in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ, QQ

    >>> R, x,y = ring("x,y", ZZ)
    >>> f = 2*x*y + 6*x + 4*y + 12

    >>> R.dmp_ground_content(f)
    2

    >>> R, x,y = ring("x,y", QQ)
    >>> f = 2*x*y + 6*x + 4*y + 12

    >>> R.dmp_ground_content(f)
    2

    """
    from sympy.polys.domains import QQ

    if not u:
        return dup_content(_dup(f), K)

    if dmp_zero_p(f, u):
        return K.zero

    cont, v = K.zero, u - 1

    if K == QQ:
        for c in f:
            cont = K.gcd(cont, dmp_ground_content(c, v, K))
    else:
        for c in f:
            cont = K.gcd(cont, dmp_ground_content(c, v, K))

            if K.is_one(cont):
                break

    return cont


def dup_primitive(f: dup[Er], K: Domain[Er]) -> tuple[Er, dup[Er]]:
    """
    Compute content and the primitive form of ``f`` in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ, QQ

    >>> R, x = ring("x", ZZ)
    >>> f = 6*x**2 + 8*x + 12

    >>> R.dup_primitive(f)
    (2, 3*x**2 + 4*x + 6)

    >>> R, x = ring("x", QQ)
    >>> f = 6*x**2 + 8*x + 12

    >>> R.dup_primitive(f)
    (2, 3*x**2 + 4*x + 6)

    """
    if not f:
        return K.zero, f

    cont = dup_content(f, K)

    if K.is_one(cont):
        return cont, f
    else:
        return cont, dup_quo_ground(f, cont, K)


def dmp_ground_primitive(f: dmp[Er], u: int, K: Domain[Er]) -> tuple[Er, dmp[Er]]:
    """
    Compute content and the primitive form of ``f`` in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ, QQ

    >>> R, x,y = ring("x,y", ZZ)
    >>> f = 2*x*y + 6*x + 4*y + 12

    >>> R.dmp_ground_primitive(f)
    (2, x*y + 3*x + 2*y + 6)

    >>> R, x,y = ring("x,y", QQ)
    >>> f = 2*x*y + 6*x + 4*y + 12

    >>> R.dmp_ground_primitive(f)
    (2, x*y + 3*x + 2*y + 6)

    """
    if not u:
        cont, fu = dup_primitive(_dup(f), K)
        return cont, _dmp(fu)

    if dmp_zero_p(f, u):
        return K.zero, f

    cont = dmp_ground_content(f, u, K)

    if K.is_one(cont):
        return cont, f
    else:
        return cont, dmp_quo_ground(f, cont, u, K)


def dup_extract(f: dup[Er], g: dup[Er], K: Domain[Er]) -> tuple[Er, dup[Er], dup[Er]]:
    """
    Extract common content from a pair of polynomials in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_extract(6*x**2 + 12*x + 18, 4*x**2 + 8*x + 12)
    (2, 3*x**2 + 6*x + 9, 2*x**2 + 4*x + 6)

    """
    fc = dup_content(f, K)
    gc = dup_content(g, K)

    gcd = K.gcd(fc, gc)

    if not K.is_one(gcd):
        f = dup_quo_ground(f, gcd, K)
        g = dup_quo_ground(g, gcd, K)

    return gcd, f, g


def dmp_ground_extract(
    f: dmp[Er], g: dmp[Er], u: int, K: Domain[Er]
) -> tuple[Er, dmp[Er], dmp[Er]]:
    """
    Extract common content from a pair of polynomials in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dmp_ground_extract(6*x*y + 12*x + 18, 4*x*y + 8*x + 12)
    (2, 3*x*y + 6*x + 9, 2*x*y + 4*x + 6)

    """
    fc = dmp_ground_content(f, u, K)
    gc = dmp_ground_content(g, u, K)

    gcd = K.gcd(fc, gc)

    if not K.is_one(gcd):
        f = dmp_quo_ground(f, gcd, u, K)
        g = dmp_quo_ground(g, gcd, u, K)

    return gcd, f, g


def dup_real_imag(f: dup[Er], K: Domain[Er]) -> tuple[dmp[Er], dmp[Er]]:
    """
    Find ``f1`` and ``f2``, such that ``f(x+I*y) = f1(x,y) + f2(x,y)*I``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dup_real_imag(x**3 + x**2 + x + 1)
    (x**3 + x**2 - 3*x*y**2 + x - y**2 + 1, 3*x**2*y + 2*x*y - y**3 + y)

    >>> from sympy.abc import x, y, z
    >>> from sympy import I
    >>> (z**3 + z**2 + z + 1).subs(z, x+I*y).expand().collect(I)
    x**3 + x**2 - 3*x*y**2 + x - y**2 + I*(3*x**2*y + 2*x*y - y**3 + y) + 1

    """
    if not K.is_ZZ and not K.is_QQ:
        raise DomainError("computing real and imaginary parts is not supported over %s" % K)

    f1 = dmp_zero(1, K)
    f2 = dmp_zero(1, K)

    if not f:
        return f1, f2

    g = _dmp2([[[K.one, K.zero]], [[K.one], []]])
    h = dmp_ground(f[0], 2, K)

    for c in f[1:]:
        h = dmp_mul(h, g, 2, K)
        h = dmp_add_term(h, dmp_ground(c, 1, K), 0, 2, K)

    H = dmp_to_raw_dict(h, 2, K)

    for k, h in H.items():
        m = k % 4

        if not m:
            f1 = dmp_add(f1, h, 1, K)
        elif m == 1:
            f2 = dmp_add(f2, h, 1, K)
        elif m == 2:
            f1 = dmp_sub(f1, h, 1, K)
        else:
            f2 = dmp_sub(f2, h, 1, K)

    return f1, f2


def dup_mirror(f: dup[Er], K: Domain[Er]) -> dup[Er]:
    """
    Evaluate efficiently the composition ``f(-x)`` in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_mirror(x**3 + 2*x**2 - 4*x + 2)
    -x**3 + 2*x**2 + 4*x + 2

    """
    f = list(f)

    for i in range(len(f) - 2, -1, -2):
        f[i] = -f[i]

    return f


def dup_scale(f: dup[Er], a: Er, K: Domain[Er]) -> dup[Er]:
    """
    Evaluate efficiently composition ``f(a*x)`` in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_scale(x**2 - 2*x + 1, ZZ(2))
    4*x**2 - 4*x + 1

    """
    f, n, b = list(f), len(f) - 1, a

    for i in range(n - 1, -1, -1):
        f[i], b = b*f[i], b*a

    return f


def dup_shift(f: dup[Er], a: Er, K: Domain[Er]) -> dup[Er]:
    """
    Evaluate efficiently Taylor shift ``f(x + a)`` in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_shift(x**2 - 2*x + 1, ZZ(2))
    x**2 + 2*x + 1

    """
    f, n = list(f), len(f) - 1

    for i in range(n, 0, -1):
        for j in range(0, i):
            f[j + 1] += a*f[j]

    return f


def dmp_shift(f: dmp[Er], a: list[Er], u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Evaluate efficiently Taylor shift ``f(X + A)`` in ``K[X]``.

    Examples
    ========

    >>> from sympy import symbols, ring, ZZ
    >>> x, y = symbols('x y')
    >>> R, _, _ = ring([x, y], ZZ)

    >>> p = x**2*y + 2*x*y + 3*x + 4*y + 5

    >>> R.dmp_shift(R(p), [ZZ(1), ZZ(2)])
    x**2*y + 2*x**2 + 4*x*y + 11*x + 7*y + 22

    >>> p.subs({x: x + 1, y: y + 2}).expand()
    x**2*y + 2*x**2 + 4*x*y + 11*x + 7*y + 22
    """
    if not u:
        return _dmp(dup_shift(_dup(f), a[0], K))

    if dmp_zero_p(f, u):
        return f

    a0, a1 = a[0], a[1:]

    if any(a1):
        f = [ dmp_shift(c, a1, u-1, K) for c in f ]
    else:
        f = list(f)

    if a0:
        n = len(f) - 1

        for i in range(n, 0, -1):
            for j in range(0, i):
                afj = dmp_mul_ground(f[j], a0, u-1, K)
                f[j + 1] = dmp_add(f[j + 1], afj, u-1, K)

    return dmp_strip(f, u, K)


def dup_transform(f: dup[Er], p: dup[Er], q: dup[Er], K: Domain[Er]) -> dup[Er]:
    """
    Evaluate functional transformation ``q**n * f(p/q)`` in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_transform(x**2 - 2*x + 1, x**2 + 1, x - 1)
    x**4 - 2*x**3 + 5*x**2 - 4*x + 4

    """
    if not f:
        return []

    n = len(f) - 1
    h, Q = [f[0]], [[K.one]]

    for i in range(0, n):
        Q.append(dup_mul(Q[-1], q, K))

    for c, q in zip(f[1:], Q[1:]):
        h = dup_mul(h, p, K)
        q = dup_mul_ground(q, c, K)
        h = dup_add(h, q, K)

    return h


def dup_compose(f: dup[Er], g: dup[Er], K: Domain[Er]) -> dup[Er]:
    """
    Evaluate functional composition ``f(g)`` in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_compose(x**2 + x, x - 1)
    x**2 - x

    """
    if len(g) <= 1:
        return dup_strip([dup_eval(f, dup_LC(g, K), K)], K)

    if not f:
        return []

    h = [f[0]]

    for c in f[1:]:
        h = dup_mul(h, g, K)
        h = dup_add_term(h, c, 0, K)

    return h


def dmp_compose(f: dmp[Er], g: dmp[Er], u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Evaluate functional composition ``f(g)`` in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x,y = ring("x,y", ZZ)

    >>> R.dmp_compose(x*y + 2*x + y, y)
    y**2 + 3*y

    """
    if not u:
        return _dmp(dup_compose(_dup(f), _dup(g), K))

    if dmp_zero_p(f, u):
        return f

    h = [f[0]]

    for c in f[1:]:
        h = dmp_mul(h, g, u, K)
        h = dmp_add_term(h, c, 0, u, K)

    return h


def _dup_series_compose(f: dup[Er], g: dup[Er], n: int, K: Domain[Er]) -> dup[Er]:
    """
    Helper function for dup_series_compose using divide and conquer.
    """
    if len(f) == 1:
        return [f[0]]

    m = len(f) // 2
    f_high = f[:-m]
    f_low = f[-m:]

    comp0 = _dup_series_compose(f_low, g, n, K)
    comp1 = _dup_series_compose(f_high, g, n, K)

    g_power = dup_series_pow(g, m, n, K)
    high_term = dup_series_mul(comp1, g_power, n, K)

    result = dup_add(high_term, comp0, K)
    return dup_truncate(result, n, K)


def dup_series_compose(f: dup[Er], g: dup[Er], n: int, K: Domain[Er]) -> dup[Er]:
    """
    Compute ``f(g(x))`` mod ``x**n`` using divide and conquer composition.

    Examples
    ========
    >>> from sympy import ZZ
    >>> from sympy.polys.densetools import dup_series_compose
    >>> from sympy.polys.densebasic import dup_from_list, dup_print
    >>> f = dup_from_list([1, 1, 1], ZZ)
    >>> g = dup_from_list([1, 1], ZZ)
    >>> comp = dup_series_compose(f, g, 3, ZZ)
    >>> dup_print(comp, 'x')
    x**2 + 3*x + 3

    """
    f = dup_truncate(f, n, K)
    g = dup_truncate(g, n, K)

    if len(g) <= 1:
        return dup_strip([dup_eval(f, dup_LC(g, K), K)], K)

    if not f:
        return []

    return _dup_series_compose(f, g, n, K)


def _dup_right_decompose(f: dup[Er], s: int, K: Domain[Er]) -> dup[Er]:
    """Helper function for :func:`_dup_decompose`."""
    n = len(f) - 1
    lc = dup_LC(f, K)

    fd: dict[int, Er] = dup_to_raw_dict(f, K)
    g = { s: K.one }

    r = n // s

    for i in range(1, s):
        coeff = K.zero

        for j in range(0, i):
            if not n + j - i in fd:
                continue

            if not s - j in g:
                continue

            fc, gc = fd[n + j - i], g[s - j]
            coeff += (i - r*j)*fc*gc

        g[s - i] = K.quo(coeff, i*r*lc)

    return dup_from_raw_dict(g, K)


def _dup_left_decompose(f: dup[Er], h: dup[Er], K: Domain[Er]) -> dup[Er] | None:
    """Helper function for :func:`_dup_decompose`."""
    g: dict[int, Er] = {}
    i = 0

    while f:
        q, r = dup_div(f, h, K)

        if dup_degree(r) > 0:
            return None
        else:
            g[i] = dup_LC(r, K)
            f, i = q, i + 1

    return dup_from_raw_dict(g, K)


def _dup_decompose(f: dup[Er], K: Domain[Er]) -> tuple[dup[Er], dup[Er]] | None:
    """Helper function for :func:`dup_decompose`."""
    df = len(f) - 1

    for s in range(2, df):
        if df % s != 0:
            continue

        h = _dup_right_decompose(f, s, K)
        g = _dup_left_decompose(f, h, K)

        if g is not None:
            return g, h

    return None


def dup_decompose(f: dup[Er], K: Domain[Er]) -> list[dup[Er]]:
    """
    Computes functional decomposition of ``f`` in ``K[x]``.

    Given a univariate polynomial ``f`` with coefficients in a field of
    characteristic zero, returns list ``[f_1, f_2, ..., f_n]``, where::

              f = f_1 o f_2 o ... f_n = f_1(f_2(... f_n))

    and ``f_2, ..., f_n`` are monic and homogeneous polynomials of at
    least second degree.

    Unlike factorization, complete functional decompositions of
    polynomials are not unique, consider examples:

    1. ``f o g = f(x + b) o (g - b)``
    2. ``x**n o x**m = x**m o x**n``
    3. ``T_n o T_m = T_m o T_n``

    where ``T_n`` and ``T_m`` are Chebyshev polynomials.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_decompose(x**4 - 2*x**3 + x**2)
    [x**2, x**2 - x]

    References
    ==========

    .. [1] [Kozen89]_

    """
    F: list[dup[Er]] = []

    while True:
        result = _dup_decompose(f, K)

        if result is not None:
            f, h = result
            F = [h] + F
        else:
            break

    return [f] + F


def dmp_alg_inject(
    f: dmp[Er],
    u: int,
    K: RingExtension[Er, Eg],
) -> tuple[dmp[Eg], int, Domain[Eg]]:
    """
    Convert polynomial from ``K(a)[X]`` to ``K[a,X]``.

    Examples
    ========

    >>> from sympy.polys.densetools import dmp_alg_inject
    >>> from sympy import QQ, sqrt

    >>> K = QQ.algebraic_field(sqrt(2))

    >>> p = [K.from_sympy(sqrt(2)), K.zero, K.one]
    >>> P, lev, dom = dmp_alg_inject(p, 0, K)
    >>> P
    [[1, 0, 0], [1]]
    >>> lev
    1
    >>> dom
    QQ

    """
    if not (K.is_Algebraic or K.is_GaussianRing or K.is_GaussianField):
        raise DomainError('computation can be done only in an algebraic domain')

    fd: dict[tuple[int, ...], Er] = dmp_to_dict(f, u, K)
    h: dict[tuple[int, ...], Eg] = {}

    for f_monom, g in fd.items():
        for g_monom, c in K.to_dict(g).items():
            h[g_monom + f_monom] = c

    F = dmp_from_dict(h, u + 1, K.dom)

    return F, u + 1, K.dom


def dmp_lift(f: dmp[Alg], u: int, K: AlgebraicField) -> dmp[MPQ]:
    """
    Convert algebraic coefficients to integers in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys import ring, QQ
    >>> from sympy import I

    >>> K = QQ.algebraic_field(I)
    >>> R, x = ring("x", K)

    >>> f = x**2 + K([QQ(1), QQ(0)])*x + K([QQ(2), QQ(0)])

    >>> R.dmp_lift(f)
    x**4 + x**2 + 4*x + 4

    """
    # Circular import. Probably dmp_lift should be moved to euclidtools
    from .euclidtools import dmp_resultant

    F, v, K2 = dmp_alg_inject(f, u, K)

    p_a = K.mod.to_list()
    P_A = dmp_include(p_a, list(range(1, v + 1)), 0, K2)

    return dmp_resultant(F, P_A, v, K2) # type: ignore


def dup_sign_variations(f, K):
    """
    Compute the number of sign variations of ``f`` in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys import ring, ZZ
    >>> R, x = ring("x", ZZ)

    >>> R.dup_sign_variations(x**4 - x**2 - x + 1)
    2

    """
    def is_negative_sympy(a):
        if not a:
            # XXX: requires zero equivalence testing in the domain
            return False
        else:
            # XXX: This is inefficient. It should not be necessary to use a
            # symbolic expression here at least for algebraic fields. If the
            # domain elements can be numerically evaluated to real values with
            # precision then this should work. We first need to rule out zero
            # elements though.
            return bool(K.to_sympy(a) < 0)

    # XXX: There should be a way to check for real numeric domains and
    # Domain.is_negative should be fixed to handle all real numeric domains.
    # It should not be necessary to special case all these different domains
    # in this otherwise generic function.
    if K.is_ZZ or K.is_QQ or K.is_RR:
        is_negative = K.is_negative
    elif K.is_AlgebraicField and K.ext.is_comparable:
        is_negative = is_negative_sympy
    elif ((K.is_PolynomialRing or K.is_FractionField) and len(K.symbols) == 1 and
          (K.dom.is_ZZ or K.dom.is_QQ or K.is_AlgebraicField) and
          K.symbols[0].is_transcendental and K.symbols[0].is_comparable):
        # We can handle a polynomial ring like QQ[E] if there is a single
        # transcendental generator because then zero equivalence is assured.
        is_negative = is_negative_sympy
    else:
        raise DomainError("sign variation counting not supported over %s" % K)

    prev, k = K.zero, 0

    for coeff in f:
        if is_negative(coeff*prev):
            k += 1

        if coeff:
            prev = coeff

    return k


def dup_clear_denoms(f, K0, K1=None, convert=False):
    """
    Clear denominators, i.e. transform ``K_0`` to ``K_1``.

    Examples
    ========

    >>> from sympy.polys import ring, QQ
    >>> R, x = ring("x", QQ)

    >>> f = QQ(1,2)*x + QQ(1,3)

    >>> R.dup_clear_denoms(f, convert=False)
    (6, 3*x + 2)
    >>> R.dup_clear_denoms(f, convert=True)
    (6, 3*x + 2)

    """
    if K1 is None:
        if K0.has_assoc_Ring:
            K1 = K0.get_ring()
        else:
            K1 = K0

    common = K1.one

    for c in f:
        common = K1.lcm(common, K0.denom(c))

    if K1.is_one(common):
        if not convert:
            return common, f
        else:
            return common, dup_convert(f, K0, K1)

    # Use quo rather than exquo to handle inexact domains by discarding the
    # remainder.
    f = [K0.numer(c)*K1.quo(common, K0.denom(c)) for c in f]

    if not convert:
        return common, dup_convert(f, K1, K0)
    else:
        return common, f


def _rec_clear_denoms(g, v, K0, K1):
    """Recursive helper for :func:`dmp_clear_denoms`."""
    common = K1.one

    if not v:
        for c in g:
            common = K1.lcm(common, K0.denom(c))
    else:
        w = v - 1

        for c in g:
            common = K1.lcm(common, _rec_clear_denoms(c, w, K0, K1))

    return common


def dmp_clear_denoms(f, u, K0, K1=None, convert=False):
    """
    Clear denominators, i.e. transform ``K_0`` to ``K_1``.

    Examples
    ========

    >>> from sympy.polys import ring, QQ
    >>> R, x,y = ring("x,y", QQ)

    >>> f = QQ(1,2)*x + QQ(1,3)*y + 1

    >>> R.dmp_clear_denoms(f, convert=False)
    (6, 3*x + 2*y + 6)
    >>> R.dmp_clear_denoms(f, convert=True)
    (6, 3*x + 2*y + 6)

    """
    if not u:
        return dup_clear_denoms(f, K0, K1, convert=convert)

    if K1 is None:
        if K0.has_assoc_Ring:
            K1 = K0.get_ring()
        else:
            K1 = K0

    common = _rec_clear_denoms(f, u, K0, K1)

    if not K1.is_one(common):
        f = dmp_mul_ground(f, common, u, K0)

    if not convert:
        return common, f
    else:
        return common, dmp_convert(f, u, K0, K1)


def dup_revert(f: dup[Er], n: int, K: Domain[Er]) -> dup[Er]:
    """
    Compute ``f**(-1)`` mod ``x**n`` using Newton iteration.

    This function computes first ``2**n`` terms of a polynomial that
    is a result of inversion of a polynomial modulo ``x**n``. This is
    useful to efficiently compute series expansion of ``1/f``.

    Examples
    ========

    >>> from sympy.polys import ring, QQ
    >>> R, x = ring("x", QQ)

    >>> f = -QQ(1,720)*x**6 + QQ(1,24)*x**4 - QQ(1,2)*x**2 + 1

    >>> R.dup_revert(f, 8)
    61/720*x**6 + 5/24*x**4 + 1/2*x**2 + 1

    """
    g = [K.revert(dup_TC(f, K))]
    h = [K.one, K.zero, K.zero]

    N = int(_ceil(_log2(n)))

    for i in range(1, N + 1):
        a = dup_mul_ground(g, K(2), K)
        b = dup_mul(f, dup_sqr(g, K), K)
        g = dup_rem(dup_sub(a, b, K), h, K)
        h = dup_lshift(h, dup_degree(h), K)

    return g


def dmp_revert(f: dmp[Er], n: int, u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Compute ``f**(-1)`` mod ``x**n`` using Newton iteration.

    Examples
    ========

    >>> from sympy.polys import ring, QQ
    >>> R, x,y = ring("x,y", QQ)

    """
    if not u:
        return _dmp(dup_revert(_dup(f), n, K))
    else:
        raise MultivariatePolynomialError(f, n)


def _dup_series_reversion_small(f: dup[Er], n: int, K: Domain[Er]) -> dup[Er]:
    """
    Helper function for :func:`dup_series_reversion`.
    ``n`` should be less than or equal to 4.
    """
    if n < 1 or n > 4:
        raise ValueError("Only n <= 4 supported")

    f = dup_truncate(f, n, K)
    f = [K.zero] * (4 - len(f)) + f
    a, b, c, d = f

    cinv = K.revert(c)
    g = [K.zero] * n

    if n >= 2:
        g[-2] = cinv

    if n >= 3:
        g[-3] = -b * cinv ** 3

    if n >= 4:
        g[-4] = (2 * b ** 2 - a * c) * cinv ** 5

    return dup_strip(g, K)


def dup_series_reversion(f: dup[Er], n: int, K: Domain[Er]) -> dup[Er]:
    r"""
    Computes the compositional inverse of f using fast lagrange inversion.
    The result is computed modulo x**n.

    Examples
    ========
    >>> from sympy.polys import QQ
    >>> from sympy.polys.densetools import dup_series_reversion
    >>> from sympy.polys.densebasic import dup_from_list, dup_print
    >>> f = dup_from_list([QQ(1, 3), QQ(1, 4), QQ(1, 5), QQ(1, 6), QQ(1, 7), 0], QQ)
    >>> rev = dup_series_reversion(f, 7, QQ)
    >>> dup_print(rev, 'x')
    5528444159/32400*x**6 + 467419477/16200*x**5 - 789929/216*x**4 + 40817/90*x**3 - 343/6*x**2 + 7*x

    References
    ==========

    .. [1] Johansson, F. A fast algorithm for reversion of power series.
        https://arxiv.org/abs/1403.4676

    """
    if not f:
        return []

    if f[-1] != K.zero:
        raise ValueError("f must have zero constant term")
    if n<=4:
        return _dup_series_reversion_small(f, n, K)

    # Step 0: h = x / f mod x**(n-1)
    f1 = dup_rshift(f, 1, K)
    f1_inv = dup_revert(f1, n, K)
    h = dup_truncate(f1_inv, n, K)

    m = _ceil(sqrt(n - 1))

    # Precompute powers of h: h, h^2, ..., h^m mod x**(n-1)
    H = [h]
    for _ in range(1, m):
        h_prev = H[-1]
        h_next = dup_series_mul(h_prev, h, n, K)
        H.append(h_next)

    # Initialize g with n zeros
    g = [K.zero] * n

    # First block: compute g[i] for i = 1 to m - 1
    for i in range(1, m):
        coeff = dup_nth(H[i-1], i - 1, K)  # coeff = [x**(i-1)](h^i)
        # (1/i) * [x**(i-1)](h^i)
        g[-(i + 1)] = K.quo(coeff, K(i))

    # t = h^m
    t = H[m - 1]

    # Loop over blocks of size m
    for i in range(m, n, m):
        # g[i] = (1 / i) * [x**(i-1)](t)
        coeff = dup_nth(t, i - 1, K)
        g[-(i + 1)] = K.quo(coeff, K(i))

        # Now fill g[i + j] for 1 <= j < m
        for j in range(1, m):
            if i + j >= n:
                break

            s = K.zero
            for k in range(0, i + j):
                c1 = dup_nth(t, k, K)
                c2 = dup_nth(H[j-1], i + j - k - 1, K)
                s += c1 * c2

            g[-(i + j + 1)] = K.quo(s, K(i + j))

        # t = t * h^m mod x**(n-1)
        t = dup_series_mul(t, H[m - 1], n, K)

    return dup_strip(g, K)

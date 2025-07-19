"""Basic tools for dense recursive polynomials in ``K[x]`` or ``K[X]``."""

from __future__ import annotations


from typing import TYPE_CHECKING, TypeVar, Iterable, Callable, Any

from sympy.core import igcd
from sympy.core.expr import Expr
from sympy.polys.domains.domain import Domain, Er, Es, Eg
from sympy.polys.domains.polynomialring import PolynomialRing
from sympy.polys.monomials import monomial_min, monomial_ldiv
from sympy.polys.orderings import monomial_key, MonomialOrder

import random


_T = TypeVar("_T")
dup: TypeAlias = "list[_T]"
dmp: TypeAlias = "list[dmp[_T]]"
dup_tup: TypeAlias = "tuple[_T, ...]"
dmp_tup: TypeAlias = "tuple[dmp_tup[_T], ...]"
monom: TypeAlias = "tuple[int, ...]"


# The _dup and _dmp functions do not do anything but are needed so that a type
# checker can understand the conversion between the two types.
#
# A dup is a list of domain elements. A dmp is a list of lists of domain
# elements of arbitrary depth.


if TYPE_CHECKING:
    from typing import TypeAlias
    from sympy.polys.rings import PolyElement
    from sympy.polys.domains.algebraicfield import AlgebraicField
    from sympy.polys.polyclasses import ANP
    Epa = TypeVar("Epa", PolyElement, ANP)

    def _dup(p: dmp[_T], /) -> dup[_T]: ...
    def _dmp(p: dup[_T], /) -> dmp[_T]: ...
    def _dmp_tup(p: tuple[_T, ...], /) -> dmp_tup[_T]: ...
    def _idup(ps: tuple[dmp[_T], ...], /) -> tuple[dup[_T], ...]: ...
    def _idmp(ps: tuple[dup[_T], ...], /) -> tuple[dmp[_T], ...]: ...
else:

    def _dup(p, /):
        return p

    def _dmp(p, /):
        return p

    def _dmp_tup(p, /):
        return p

    def _idup(ps, /):
        return ps

    def _idmp(ps, /):
        return ps


# XXX: This causes lots of type: ignore. It would be better just to use -1:
ninf = float("-inf")


def dup_LC(f: dup[Er], K: Domain[Er]) -> Er:
    """
    Return the leading coefficient of ``f``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dup_LC

    >>> dup_LC([1, 2, 3], ZZ)
    1

    """
    if not f:
        return K.zero
    else:
        return f[0]


def dup_TC(f: dup[Er], K: Domain[Er]) -> Er:
    """
    Return the trailing coefficient of ``f``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dup_TC

    >>> dup_TC([1, 2, 3], ZZ)
    3

    """
    if not f:
        return K.zero
    else:
        return f[-1]


def dmp_LC(f: dmp[Er], K: Domain[Er]) -> dmp[Er]:
    """
    Return the leading coefficient of ``f``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_LC

    >>> f = ZZ.map([[1], [2, 3]])

    >>> dmp_LC(f, ZZ)
    [1]

    """
    if not f:
        # XXX: Remove this. It should not be needed since a zero dmp is
        # represented like [[[]]].
        return K.zero # type: ignore
    else:
        return f[0]


def dmp_TC(f: dmp[Er], K: Domain[Er]) -> dmp[Er]:
    """
    Return the trailing coefficient of ``f``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_TC

    >>> f = ZZ.map([[1], [2, 3]])

    >>> dmp_TC(f, ZZ)
    [2, 3]

    """
    if not f:
        return K.zero # type: ignore
    else:
        return f[-1]


poly_LC = dup_LC
poly_TC = dup_TC


def dmp_ground_LC(f: dmp[Er], u: int, K: Domain[Er]) -> Er:
    """
    Return the ground leading coefficient.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_ground_LC

    >>> f = ZZ.map([[[1], [2, 3]]])

    >>> dmp_ground_LC(f, 2, ZZ)
    1

    """
    while u:
        f = dmp_LC(f, K)
        u -= 1

    return dup_LC(_dup(f), K)


def dmp_ground_TC(f: dmp[Er], u: int, K: Domain[Er]) -> Er:
    """
    Return the ground trailing coefficient.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_ground_TC

    >>> f = ZZ.map([[[1], [2, 3]]])

    >>> dmp_ground_TC(f, 2, ZZ)
    3

    """
    while u:
        f = dmp_TC(f, K)
        u -= 1

    return dup_TC(_dup(f), K)


def dmp_true_LT(f: dmp[Er], u: int, K: Domain[Er]) -> tuple[monom, Er]:
    """
    Return the leading term ``c * x_1**n_1 ... x_k**n_k``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_true_LT

    >>> f = ZZ.map([[4], [2, 0], [3, 0, 0]])

    >>> dmp_true_LT(f, 1, ZZ)
    ((2, 0), 4)

    """
    monom = []

    while u:
        monom.append(len(f) - 1)
        f, u = f[0], u - 1

    if not f:
        monom.append(0)
    else:
        monom.append(len(f) - 1)

    return tuple(monom), dup_LC(_dup(f), K)


def dup_degree(f: dup[Er]) -> int | float:
    """
    Return the leading degree of ``f`` in ``K[x]``.

    Note that the degree of 0 is negative infinity (``float('-inf')``).

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dup_degree

    >>> f = ZZ.map([1, 2, 0, 3])

    >>> dup_degree(f)
    3

    """
    if not f:
        return ninf
    return len(f) - 1


def dmp_degree(f: dmp[Er], u: int) -> int | float:
    """
    Return the leading degree of ``f`` in ``x_0`` in ``K[X]``.

    Note that the degree of 0 is negative infinity (``float('-inf')``).

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_degree

    >>> dmp_degree([[[]]], 2)
    -inf

    >>> f = ZZ.map([[2], [1, 2, 3]])

    >>> dmp_degree(f, 1)
    1

    """
    if dmp_zero_p(f, u):
        return ninf
    else:
        return len(f) - 1


def _rec_degree_in(g: dmp[Er], v: int, i: int, j: int) -> int | float:
    """Recursive helper function for :func:`dmp_degree_in`."""
    if i == j:
        return dmp_degree(g, v)

    v, i = v - 1, i + 1

    return max(_rec_degree_in(c, v, i, j) for c in g)


def dmp_degree_in(f: dmp[Er], j: int, u: int) -> int | float:
    """
    Return the leading degree of ``f`` in ``x_j`` in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_degree_in

    >>> f = ZZ.map([[2], [1, 2, 3]])

    >>> dmp_degree_in(f, 0, 1)
    1
    >>> dmp_degree_in(f, 1, 1)
    2

    """
    if not j:
        return dmp_degree(f, u)
    if j < 0 or j > u:
        raise IndexError("0 <= j <= %s expected, got %s" % (u, j))

    return _rec_degree_in(f, u, 0, j)


def _rec_degree_list(g: dmp[Er], v: int, i: int, degs: list[int | float]) -> None:
    """Recursive helper for :func:`dmp_degree_list`."""
    degs[i] = max(degs[i], dmp_degree(g, v))

    if v > 0:
        v, i = v - 1, i + 1

        for c in g:
            _rec_degree_list(c, v, i, degs)


def dmp_degree_list(f: dmp[Er], u: int) -> tuple[int | float, ...]:
    """
    Return a list of degrees of ``f`` in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_degree_list

    >>> f = ZZ.map([[1], [1, 2, 3]])

    >>> dmp_degree_list(f, 1)
    (1, 2)

    """
    degs = [ninf] * (u + 1)
    _rec_degree_list(f, u, 0, degs)
    return tuple(degs)


def dup_strip(f: dup[Er]) -> dup[Er]:
    """
    Remove leading zeros from ``f`` in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys.densebasic import dup_strip

    >>> dup_strip([0, 0, 1, 2, 3, 0])
    [1, 2, 3, 0]

    """
    if not f or f[0]:
        return f

    i = 0

    for cf in f:
        if cf:
            break
        else:
            i += 1

    return f[i:]


def dmp_strip(f: dmp[Er], u: int) -> dmp[Er]:
    """
    Remove leading zeros from ``f`` in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys.densebasic import dmp_strip

    >>> dmp_strip([[], [0, 1, 2], [1]], 1)
    [[0, 1, 2], [1]]

    """
    if not u:
        return _dmp(dup_strip(_dup(f)))

    if dmp_zero_p(f, u):
        return f

    i, v = 0, u - 1

    for c in f:
        if not dmp_zero_p(c, v):
            break
        else:
            i += 1

    if i == len(f):
        return dmp_zero(u)
    else:
        return f[i:]


def _rec_validate(
    f: dmp[Er], g: dmp[Er] | Er, i: int, K: Domain[Er] | None
) -> set[int]:
    """Recursive helper for :func:`dmp_validate`."""
    if not isinstance(g, list):
        if K is not None and not K.of_type(g):
            raise TypeError("%s in %s in not of type %s" % (g, f, K.dtype))

        return {i - 1}
    elif not g:
        return {i}
    else:
        levels = set()

        for c in g:
            levels |= _rec_validate(f, c, i + 1, K)

        return levels


def _rec_strip(g: dmp[Er], v: int) -> dmp[Er]:
    """Recursive helper for :func:`_rec_strip`."""
    if not v:
        return _dmp(dup_strip(_dup(g)))

    w = v - 1

    return dmp_strip([_rec_strip(c, w) for c in g], v)


def dmp_validate(f: dmp[Er], K: Domain[Er] | None = None) -> tuple[dmp[Er], int]:
    """
    Return the number of levels in ``f`` and recursively strip it.

    Examples
    ========

    >>> from sympy.polys.densebasic import dmp_validate

    >>> dmp_validate([[], [0, 1, 2], [1]])
    ([[1, 2], [1]], 1)

    >>> dmp_validate([[1], 1])
    Traceback (most recent call last):
    ...
    ValueError: invalid data structure for a multivariate polynomial

    """
    levels = _rec_validate(f, f, 0, K)

    u = levels.pop()

    if not levels:
        return _rec_strip(f, u), u
    else:
        raise ValueError("invalid data structure for a multivariate polynomial")


def dup_reverse(f: dup[Er]) -> dup[Er]:
    """
    Compute ``x**n * f(1/x)``, i.e.: reverse ``f`` in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dup_reverse

    >>> f = ZZ.map([1, 2, 3, 0])

    >>> dup_reverse(f)
    [3, 2, 1]

    """
    return dup_strip(list(reversed(f)))


def dup_copy(f: dup[Er]) -> dup[Er]:
    """
    Create a new copy of a polynomial ``f`` in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dup_copy

    >>> f = ZZ.map([1, 2, 3, 0])

    >>> dup_copy([1, 2, 3, 0])
    [1, 2, 3, 0]

    """
    return list(f)


def dmp_copy(f: dmp[Er], u: int) -> dmp[Er]:
    """
    Create a new copy of a polynomial ``f`` in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_copy

    >>> f = ZZ.map([[1], [1, 2]])

    >>> dmp_copy(f, 1)
    [[1], [1, 2]]

    """
    if not u:
        return list(f)

    v = u - 1

    return [dmp_copy(c, v) for c in f]


def dup_to_tuple(f: dup[Er]) -> tuple[Er, ...]:
    """
    Convert `f` into a tuple.

    This is needed for hashing. This is similar to dup_copy().

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dup_copy

    >>> f = ZZ.map([1, 2, 3, 0])

    >>> dup_copy([1, 2, 3, 0])
    [1, 2, 3, 0]

    """
    return tuple(f)


def dmp_to_tuple(f: dmp[Er], u: int) -> dmp_tup[Er]:
    """
    Convert `f` into a nested tuple of tuples.

    This is needed for hashing.  This is similar to dmp_copy().

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_to_tuple

    >>> f = ZZ.map([[1], [1, 2]])

    >>> dmp_to_tuple(f, 1)
    ((1,), (1, 2))

    """
    if not u:
        return _dmp_tup(tuple(_dup(f)))
    v = u - 1

    return tuple(dmp_to_tuple(c, v) for c in f)


def dup_normal(f: dup[Er], K: Domain[Er]) -> dup[Er]:
    """
    Normalize univariate polynomial in the given domain.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dup_normal

    >>> dup_normal([0, 1, 2, 3], ZZ)
    [1, 2, 3]

    """
    return dup_strip([K.normal(c) for c in f])


def dmp_normal(f: dmp[Er], u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Normalize a multivariate polynomial in the given domain.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_normal

    >>> dmp_normal([[], [0, 1, 2]], 1, ZZ)
    [[1, 2]]

    """
    if not u:
        return _dmp(dup_normal(_dup(f), K))

    v = u - 1

    return dmp_strip([dmp_normal(c, v, K) for c in f], u)


def dup_convert(f: dup[Er], K0: Domain[Er] | None, K1: Domain[Es]) -> dup[Es]:
    """
    Convert the ground domain of ``f`` from ``K0`` to ``K1``.

    Examples
    ========

    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dup_convert

    >>> R, x = ring("x", ZZ)

    >>> dup_convert([R(1), R(2)], R.to_domain(), ZZ)
    [1, 2]
    >>> dup_convert([ZZ(1), ZZ(2)], ZZ, R.to_domain())
    [1, 2]

    """
    if K0 is not None and K0 == K1:
        return f  # type: ignore
    else:
        return dup_strip([K1.convert(c, K0) for c in f])


def dmp_convert(f: dmp[Er], u: int, K0: Domain[Er] | None, K1: Domain[Es]) -> dmp[Es]:
    """
    Convert the ground domain of ``f`` from ``K0`` to ``K1``.

    Examples
    ========

    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_convert

    >>> R, x = ring("x", ZZ)

    >>> dmp_convert([[R(1)], [R(2)]], 1, R.to_domain(), ZZ)
    [[1], [2]]
    >>> dmp_convert([[ZZ(1)], [ZZ(2)]], 1, ZZ, R.to_domain())
    [[1], [2]]

    """
    if not u:
        return _dmp(dup_convert(_dup(f), K0, K1))
    if K0 is not None and K0 == K1:
        return f  # type: ignore

    v = u - 1

    return dmp_strip([dmp_convert(c, v, K0, K1) for c in f], u)


def dup_from_sympy(f: list[Expr], K: Domain[Er]) -> dup[Er]:
    """
    Convert the ground domain of ``f`` from SymPy to ``K``.

    Examples
    ========

    >>> from sympy import S
    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dup_from_sympy

    >>> dup_from_sympy([S(1), S(2)], ZZ) == [ZZ(1), ZZ(2)]
    True

    """
    return dup_strip([K.from_sympy(c) for c in f])


def dmp_from_sympy(f: dmp[Expr], u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Convert the ground domain of ``f`` from SymPy to ``K``.

    Examples
    ========

    >>> from sympy import S
    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_from_sympy

    >>> dmp_from_sympy([[S(1)], [S(2)]], 1, ZZ) == [[ZZ(1)], [ZZ(2)]]
    True

    """
    if not u:
        return _dmp(dup_from_sympy(_dup(f), K))

    v = u - 1

    return dmp_strip([dmp_from_sympy(c, v, K) for c in f], u)


def dup_nth(f: dup[Er], n: int, K: Domain[Er]) -> Er:
    """
    Return the ``n``-th coefficient of ``f`` in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dup_nth

    >>> f = ZZ.map([1, 2, 3])

    >>> dup_nth(f, 0, ZZ)
    3
    >>> dup_nth(f, 4, ZZ)
    0

    """
    if n < 0:
        raise IndexError("'n' must be non-negative, got %i" % n)
    elif n >= len(f):
        return K.zero
    else:
        return f[dup_degree(f) - n]  # type: ignore


def dmp_nth(f: dmp[Er], n: int, u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Return the ``n``-th coefficient of ``f`` in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_nth

    >>> f = ZZ.map([[1], [2], [3]])

    >>> dmp_nth(f, 0, 1, ZZ)
    [3]
    >>> dmp_nth(f, 4, 1, ZZ)
    []

    """
    if n < 0:
        raise IndexError("'n' must be non-negative, got %i" % n)
    elif n >= len(f):
        return dmp_zero(u - 1)
    else:
        return f[dmp_degree(f, u) - n]  # type: ignore


def dmp_ground_nth(f: dmp[Er], N: Iterable[int], u: int, K: Domain[Er]) -> Er:
    """
    Return the ground ``n``-th coefficient of ``f`` in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_ground_nth

    >>> f = ZZ.map([[1], [2, 3]])

    >>> dmp_ground_nth(f, (0, 1), 1, ZZ)
    2

    """
    v = u

    for n in N:
        if n < 0:
            raise IndexError("`n` must be non-negative, got %i" % n)
        elif n >= len(f):
            return K.zero
        else:
            d = dmp_degree(f, v)
            if d == ninf:
                d = -1
            f, v = f[d - n], v - 1  # type: ignore

    return f  # type: ignore


def dmp_zero_p(f: dmp[Er], u: int) -> bool:
    """
    Return ``True`` if ``f`` is zero in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys.densebasic import dmp_zero_p

    >>> dmp_zero_p([[[[[]]]]], 4)
    True
    >>> dmp_zero_p([[[[[1]]]]], 4)
    False

    """
    while u:
        if len(f) != 1:
            return False

        f = f[0]
        u -= 1

    return not f


def dmp_zero(u: int) -> dmp:
    """
    Return a multivariate zero.

    Examples
    ========

    >>> from sympy.polys.densebasic import dmp_zero

    >>> dmp_zero(4)
    [[[[[]]]]]

    """
    r: dmp = []

    for i in range(u):
        r = [r]

    return r


def dmp_one_p(f: dmp[Er], u: int, K: Domain[Er]) -> bool:
    """
    Return ``True`` if ``f`` is one in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_one_p

    >>> dmp_one_p([[[ZZ(1)]]], 2, ZZ)
    True

    """
    return dmp_ground_p(f, K.one, u)


def dmp_one(u: int, K: Domain[Er]) -> dmp[Er] | Er:
    """
    Return a multivariate one over ``K``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_one

    >>> dmp_one(2, ZZ)
    [[[1]]]

    """
    return dmp_ground(K.one, u)


def dmp_ground_p(f: dmp[Er], c: Er | None, u: int) -> bool:
    """
    Return True if ``f`` is constant in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys.densebasic import dmp_ground_p

    >>> dmp_ground_p([[[3]]], 3, 2)
    True
    >>> dmp_ground_p([[[4]]], None, 2)
    True

    """
    if c is not None and not c:
        return dmp_zero_p(f, u)

    while u:
        if len(f) != 1:
            return False
        f = f[0]
        u -= 1

    if c is None:
        return len(f) <= 1
    else:
        return f == [c]


def dmp_ground(c: Er, u: int) -> dmp[Er] | Er:
    """
    Return a multivariate constant.

    Examples
    ========

    >>> from sympy.polys.densebasic import dmp_ground

    >>> dmp_ground(3, 5)
    [[[[[[3]]]]]]
    >>> dmp_ground(1, -1)
    1

    """
    if not c:
        return dmp_zero(u)

    if u < 0:
        return c

    f = _dmp([c])

    for i in range(u):
        f = [f]

    return f


def dmp_zeros(n: int, u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Return a list of multivariate zeros.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_zeros

    >>> dmp_zeros(3, 2, ZZ)
    [[[[]]], [[[]]], [[[]]]]
    >>> dmp_zeros(3, -1, ZZ)
    [0, 0, 0]

    """
    if not n:
        return []

    if u < 0:
        return _dmp([K.zero] * n)
    else:
        return [dmp_zero(u) for i in range(n)]


def dmp_grounds(c: Er, n: int, u: int) -> list[dmp[Er]] | list[Er]:
    """
    Return a list of multivariate constants.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_grounds

    >>> dmp_grounds(ZZ(4), 3, 2)
    [[[[4]]], [[[4]]], [[[4]]]]
    >>> dmp_grounds(ZZ(4), 3, -1)
    [4, 4, 4]

    """
    if not n:
        return []

    if u < 0:
        return _dmp([c]) * n
    else:
        return [dmp_ground(c, u) for i in range(n)]  # type: ignore


def dmp_negative_p(f: dmp[Er], u: int, K: Domain[Er]) -> bool:
    """
    Return ``True`` if ``LC(f)`` is negative.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_negative_p

    >>> dmp_negative_p([[ZZ(1)], [-ZZ(1)]], 1, ZZ)
    False
    >>> dmp_negative_p([[-ZZ(1)], [ZZ(1)]], 1, ZZ)
    True

    """
    return K.is_negative(dmp_ground_LC(f, u, K))


def dmp_positive_p(f: dmp[Er], u: int, K: Domain[Er]) -> bool:
    """
    Return ``True`` if ``LC(f)`` is positive.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_positive_p

    >>> dmp_positive_p([[ZZ(1)], [-ZZ(1)]], 1, ZZ)
    True
    >>> dmp_positive_p([[-ZZ(1)], [ZZ(1)]], 1, ZZ)
    False

    """
    return K.is_positive(dmp_ground_LC(f, u, K))


def dup_from_dict(f: dict[tuple[int], Er] | dict[int, Er], K: Domain[Er]) -> dup[Er]:
    """
    Create a ``K[x]`` polynomial from a ``dict``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dup_from_dict

    >>> dup_from_dict({(0,): ZZ(7), (2,): ZZ(5), (4,): ZZ(1)}, ZZ)
    [1, 0, 5, 0, 7]
    >>> dup_from_dict({}, ZZ)
    []

    """
    if not f:
        return []

    n, h = max(f.keys()), []

    if isinstance(n, int):
        for k in range(n, -1, -1):
            h.append(f.get(k, K.zero))  # type: ignore
    else:
        (n,) = n

        for k in range(n, -1, -1):
            h.append(f.get((k,), K.zero))  # type: ignore

    return dup_strip(h)


def dup_from_raw_dict(f: dict[int, Er], K: Domain[Er]) -> dup[Er]:
    """
    Create a ``K[x]`` polynomial from a raw ``dict``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dup_from_raw_dict

    >>> dup_from_raw_dict({0: ZZ(7), 2: ZZ(5), 4: ZZ(1)}, ZZ)
    [1, 0, 5, 0, 7]

    """
    if not f:
        return []

    n, h = max(f.keys()), []

    for k in range(n, -1, -1):
        h.append(f.get(k, K.zero))

    return dup_strip(h)


def dmp_from_dict(f: dict[tuple[int, ...], Er], u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Create a ``K[X]`` polynomial from a ``dict``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_from_dict

    >>> dmp_from_dict({(0, 0): ZZ(3), (0, 1): ZZ(2), (2, 1): ZZ(1)}, 1, ZZ)
    [[1, 0], [], [2, 3]]
    >>> dmp_from_dict({}, 0, ZZ)
    []

    """
    if not u:
        return _dmp(dup_from_dict(f, K))  # type: ignore
    if not f:
        return dmp_zero(u)

    coeffs: dict[int, dict[monom, Er]] = {}

    for monom, coeff in f.items():
        head, tail = monom[0], monom[1:]

        if head in coeffs:
            coeffs[head][tail] = coeff
        else:
            coeffs[head] = {tail: coeff}

    n, v, h = max(coeffs.keys()), u - 1, []

    for k in range(n, -1, -1):
        dcoeff = coeffs.get(k)

        if dcoeff is not None:
            h.append(dmp_from_dict(dcoeff, v, K))
        else:
            h.append(dmp_zero(v))

    return dmp_strip(h, u)


def dup_to_dict(
    f: dup[Er], K: Domain[Er] | None = None, zero: bool = False
) -> dict[tuple[int], Er]:
    """
    Convert ``K[x]`` polynomial to a ``dict``.

    Examples
    ========

    >>> from sympy.polys.densebasic import dup_to_dict

    >>> dup_to_dict([1, 0, 5, 0, 7])
    {(0,): 7, (2,): 5, (4,): 1}
    >>> dup_to_dict([])
    {}

    """
    if not f and zero:
        return {(0,): K.zero}  # type: ignore

    n, result = len(f) - 1, {}

    for k in range(0, n + 1):
        if f[n - k]:
            result[(k,)] = f[n - k]

    return result


def dup_to_raw_dict(
    f: dup[Er], K: Domain[Er] | None = None, zero: bool = False
) -> dict[int, Er]:
    """
    Convert a ``K[x]`` polynomial to a raw ``dict``.

    Examples
    ========

    >>> from sympy.polys.densebasic import dup_to_raw_dict

    >>> dup_to_raw_dict([1, 0, 5, 0, 7])
    {0: 7, 2: 5, 4: 1}

    """
    if not f and zero:
        return {0: K.zero}  # type: ignore

    n, result = len(f) - 1, {}

    for k in range(0, n + 1):
        if f[n - k]:
            result[k] = f[n - k]

    return result


def dmp_to_dict(
    f: dmp[Er], u: int, K: Domain[Er] | None = None, zero: bool = False
) -> dict[tuple[int, ...], Er]:
    """
    Convert a ``K[X]`` polynomial to a ``dict````.

    Examples
    ========

    >>> from sympy.polys.densebasic import dmp_to_dict

    >>> dmp_to_dict([[1, 0], [], [2, 3]], 1)
    {(0, 0): 3, (0, 1): 2, (2, 1): 1}
    >>> dmp_to_dict([], 0)
    {}

    """
    if not u:
        return dup_to_dict(_dup(f), K, zero=zero)  # type: ignore

    if dmp_zero_p(f, u) and zero:
        return {(0,) * (u + 1): K.zero}  # type: ignore

    n, v, result = dmp_degree(f, u), u - 1, {}

    if n == ninf:
        n = -1

    for k in range(0, n + 1):  # type: ignore
        h = dmp_to_dict(f[n - k], v)  # type: ignore

        for exp, coeff in h.items():
            result[(k,) + exp] = coeff

    return result


def dmp_swap(f: dmp[Er], i: int, j: int, u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Transform ``K[..x_i..x_j..]`` to ``K[..x_j..x_i..]``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_swap

    >>> f = ZZ.map([[[2], [1, 0]], []])

    >>> dmp_swap(f, 0, 1, 2, ZZ)
    [[[2], []], [[1, 0], []]]
    >>> dmp_swap(f, 1, 2, 2, ZZ)
    [[[1], [2, 0]], [[]]]
    >>> dmp_swap(f, 0, 2, 2, ZZ)
    [[[1, 0]], [[2, 0], []]]

    """
    if i < 0 or j < 0 or i > u or j > u:
        raise IndexError("0 <= i < j <= %s expected" % u)
    elif i == j:
        return f

    F: dict[monom, Er] = dmp_to_dict(f, u)
    H: dict[monom, Er] = {}

    for exp, coeff in F.items():
        H[exp[:i] + (exp[j],) + exp[i + 1 : j] + (exp[i],) + exp[j + 1 :]] = coeff

    return dmp_from_dict(H, u, K)


def dmp_permute(f: dmp[Er], P: list[int], u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Return a polynomial in ``K[x_{P(1)},..,x_{P(n)}]``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_permute

    >>> f = ZZ.map([[[2], [1, 0]], []])

    >>> dmp_permute(f, [1, 0, 2], 2, ZZ)
    [[[2], []], [[1, 0], []]]
    >>> dmp_permute(f, [1, 2, 0], 2, ZZ)
    [[[1], []], [[2, 0], []]]

    """
    F: dict[monom, Er] = dmp_to_dict(f, u)
    H: dict[monom, Er] = {}

    for exp, coeff in F.items():
        new_exp = [0] * len(exp)

        for e, p in zip(exp, P):
            new_exp[p] = e

        H[tuple(new_exp)] = coeff

    return dmp_from_dict(H, u, K)


def dmp_nest(f: dmp[Er], l: int, K: Domain[Er]) -> dmp[Er]:
    """
    Return a multivariate value nested ``l``-levels.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_nest

    >>> dmp_nest([[ZZ(1)]], 2, ZZ)
    [[[[1]]]]

    """
    if not isinstance(f, list):
        return dmp_ground(f, l)

    for i in range(l):
        f = [f]

    return f


def dmp_raise(f: dmp[Er], l: int, u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Return a multivariate polynomial raised ``l``-levels.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_raise

    >>> f = ZZ.map([[], [1, 2]])

    >>> dmp_raise(f, 2, 1, ZZ)
    [[[[]]], [[[1]], [[2]]]]

    """
    if not l:
        return f

    if not u:
        if not f:
            return dmp_zero(l)

        k = l - 1

        return [dmp_ground(c, k) for c in _dup(f)]  # type: ignore

    v = u - 1

    return [dmp_raise(cp, l, v, K) for cp in f]


def dup_deflate(f: dup[Er], K: Domain[Er]) -> tuple[int, dup[Er]]:
    """
    Map ``x**m`` to ``y`` in a polynomial in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dup_deflate

    >>> f = ZZ.map([1, 0, 0, 1, 0, 0, 1])

    >>> dup_deflate(f, ZZ)
    (3, [1, 1, 1])

    """
    if dup_degree(f) <= 0:
        return 1, f

    g = 0

    for i in range(len(f)):
        if not f[-i - 1]:
            continue

        g = igcd(g, i)

        if g == 1:
            return 1, f

    return g, f[::g]


def dmp_deflate(f: dmp[Er], u: int, K: Domain[Er]) -> tuple[tuple[int, ...], dmp[Er]]:
    """
    Map ``x_i**m_i`` to ``y_i`` in a polynomial in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_deflate

    >>> f = ZZ.map([[1, 0, 0, 2], [], [3, 0, 0, 4]])

    >>> dmp_deflate(f, 1, ZZ)
    ((2, 3), [[1, 2], [3, 4]])

    """
    if dmp_zero_p(f, u):
        return (1,) * (u + 1), f

    F: dict[monom, Er] = dmp_to_dict(f, u)
    B = [0] * (u + 1)

    for M in F.keys():
        for i, m in enumerate(M):
            B[i] = igcd(B[i], m)

    for i, b in enumerate(B):
        if not b:
            B[i] = 1

    Bt = tuple(B)

    if all(b == 1 for b in Bt):
        return Bt, f

    H = {}

    for A, coeff in F.items():
        N = [a // b for a, b in zip(A, Bt)]
        H[tuple(N)] = coeff

    return Bt, dmp_from_dict(H, u, K)


def dup_multi_deflate(
    polys: tuple[dup[Er], ...], K: Domain[Er]
) -> tuple[int, tuple[dup[Er], ...]]:
    """
    Map ``x**m`` to ``y`` in a set of polynomials in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dup_multi_deflate

    >>> f = ZZ.map([1, 0, 2, 0, 3])
    >>> g = ZZ.map([4, 0, 0])

    >>> dup_multi_deflate((f, g), ZZ)
    (2, ([1, 2, 3], [4, 0]))

    """
    G = 0

    for p in polys:
        if dup_degree(p) <= 0:
            return 1, polys

        g = 0

        for i in range(len(p)):
            if not p[-i - 1]:
                continue

            g = igcd(g, i)

            if g == 1:
                return 1, polys

        G = igcd(G, g)

    return G, tuple([p[::G] for p in polys])


def dmp_multi_deflate(
    polys: tuple[dmp[Er], ...], u: int, K: Domain[Er]
) -> tuple[tuple[int, ...], tuple[dmp[Er], ...]]:
    """
    Map ``x_i**m_i`` to ``y_i`` in a set of polynomials in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_multi_deflate

    >>> f = ZZ.map([[1, 0, 0, 2], [], [3, 0, 0, 4]])
    >>> g = ZZ.map([[1, 0, 2], [], [3, 0, 4]])

    >>> dmp_multi_deflate((f, g), 1, ZZ)
    ((2, 1), ([[1, 0, 0, 2], [3, 0, 0, 4]], [[1, 0, 2], [3, 0, 4]]))

    """
    if not u:
        M, H = dup_multi_deflate(_idup(polys), K)
        return (M,), _idmp(H)

    F, B = [], [0] * (u + 1)

    for p in polys:
        f: dict[monom, Er] = dmp_to_dict(p, u)

        if not dmp_zero_p(p, u):
            for m in f.keys():
                for i, e in enumerate(m):
                    B[i] = igcd(B[i], e)

        F.append(f)

    for i, b in enumerate(B):
        if not b:
            B[i] = 1

    Bt = tuple(B)

    if all(b == 1 for b in Bt):
        return Bt, polys

    H2: list[dmp[Er]] = []

    for f in F:
        h = {}

        for A, coeff in f.items():
            N = [a // b for a, b in zip(A, Bt)]
            h[tuple(N)] = coeff

        H2.append(dmp_from_dict(h, u, K))

    return Bt, tuple(H2)


def dup_inflate(f: dup[Er], m: int, K: Domain[Er]) -> dup[Er]:
    """
    Map ``y`` to ``x**m`` in a polynomial in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dup_inflate

    >>> f = ZZ.map([1, 1, 1])

    >>> dup_inflate(f, 3, ZZ)
    [1, 0, 0, 1, 0, 0, 1]

    """
    if m <= 0:
        raise IndexError("'m' must be positive, got %s" % m)
    if m == 1 or not f:
        return f

    result = [f[0]]

    for coeff in f[1:]:
        result.extend([K.zero] * (m - 1))
        result.append(coeff)

    return result


def _rec_inflate(g: dmp[Er], M: list[int], v: int, i: int, K: Domain[Er]) -> dmp[Er]:
    """Recursive helper for :func:`dmp_inflate`."""
    if not v:
        return _dmp(dup_inflate(_dup(g), M[i], K))
    if M[i] <= 0:
        raise IndexError("all M[i] must be positive, got %s" % M[i])

    w, j = v - 1, i + 1

    g = [_rec_inflate(c, M, w, j, K) for c in g]

    result = [g[0]]

    for coeff in g[1:]:
        for _ in range(1, M[i]):
            result.append(dmp_zero(w))

        result.append(coeff)

    return result


def dmp_inflate(f: dmp[Er], M: list[int], u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Map ``y_i`` to ``x_i**k_i`` in a polynomial in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_inflate

    >>> f = ZZ.map([[1, 2], [3, 4]])

    >>> dmp_inflate(f, (2, 3), 1, ZZ)
    [[1, 0, 0, 2], [], [3, 0, 0, 4]]

    """
    if not u:
        return _dmp(dup_inflate(_dup(f), M[0], K))

    if all(m == 1 for m in M):
        return f
    else:
        return _rec_inflate(f, M, u, 0, K)


def dmp_exclude(f: dmp[Er], u: int, K: Domain[Er]) -> tuple[list[int], dmp[Er], int]:
    """
    Exclude useless levels from ``f``.

    Return the levels excluded, the new excluded ``f``, and the new ``u``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_exclude

    >>> f = ZZ.map([[[1]], [[1], [2]]])

    >>> dmp_exclude(f, 2, ZZ)
    ([2], [[1], [1, 2]], 1)

    """
    if not u or dmp_ground_p(f, None, u):
        return [], f, u

    J: list[int] = []
    F: dict[monom, Er] = dmp_to_dict(f, u)

    for j in range(0, u + 1):
        for monom in F.keys():
            if monom[j]:
                break
        else:
            J.append(j)

    if not J:
        return [], f, u

    d = {}

    for monom, coeff in F.items():
        lmonom = list(monom)

        for j in reversed(J):
            del lmonom[j]

        d[tuple(lmonom)] = coeff

    u -= len(J)

    return J, dmp_from_dict(d, u, K), u


def dmp_include(f: dmp[Er], J: list[int], u: int, K: Domain[Er]) -> dmp[Er]:
    """
    Include useless levels in ``f``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_include

    >>> f = ZZ.map([[1], [1, 2]])

    >>> dmp_include(f, [2], 1, ZZ)
    [[[1]], [[1], [2]]]

    """
    if not J:
        return f

    F: dict[tuple[int, ...], Er] = dmp_to_dict(f, u)
    d: dict[tuple[int, ...], Er] = {}

    for monom, coeff in F.items():
        lmonom = list(monom)

        for j in J:
            lmonom.insert(j, 0)

        d[tuple(lmonom)] = coeff

    u += len(J)

    return dmp_from_dict(d, u, K)


# XXX: K could be a PolynomialRing or an AlgebraicField or ...
def dmp_inject(
    f: dmp[Epa]
    , u: int, K: PolynomialRing[Eg] | AlgebraicField[Epa, Eg]
    , front: bool = False
) -> tuple[dmp[Eg], int]:
    """
    Convert ``f`` from ``K[X][Y]`` to ``K[X,Y]``.

    Examples
    ========

    >>> from sympy.polys.rings import ring
    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_inject

    >>> R, x,y = ring("x,y", ZZ)

    >>> dmp_inject([R(1), x + 2], 0, R.to_domain())
    ([[[1]], [[1], [2]]], 2)
    >>> dmp_inject([R(1), x + 2], 0, R.to_domain(), front=True)
    ([[[1]], [[1, 2]]], 2)

    """
    d: dict[monom, Epa]
    h: dict[monom, Eg]

    d = dmp_to_dict(f, u)
    h = {}

    v = K.ngens - 1

    for f_monom, g in d.items():
        g = g.to_dict()

        for g_monom, c in g.items():
            if front:
                h[g_monom + f_monom] = c
            else:
                h[f_monom + g_monom] = c

    w = u + v + 1

    return dmp_from_dict(h, w, K.dom), w


def dmp_eject(
    f: dmp[Er], u: int, K: PolynomialRing[Er], front: bool = False
) -> dmp[PolyElement[Er]]:
    """
    Convert ``f`` from ``K[X,Y]`` to ``K[X][Y]``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_eject

    >>> dmp_eject([[[1]], [[1], [2]]], 2, ZZ['x', 'y'])
    [1, x + 2]

    """
    d: dict[monom, Er]
    h: dict[monom, dict[monom, Er]]

    d = dmp_to_dict(f, u)
    h = {}

    n = K.ngens
    v = u - K.ngens + 1

    for monom, c in d.items():
        if front:
            g_monom, f_monom = monom[:n], monom[n:]
        else:
            g_monom, f_monom = monom[-n:], monom[:-n]

        if f_monom in h:
            h[f_monom][g_monom] = c
        else:
            h[f_monom] = {g_monom: c}

    g = {monom: K(c) for monom, c in h.items()}

    return dmp_from_dict(g, v - 1, K)


def dup_terms_gcd(f: dup[Er], K: Domain[Er]) -> tuple[int, dup[Er]]:
    """
    Remove GCD of terms from ``f`` in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dup_terms_gcd

    >>> f = ZZ.map([1, 0, 1, 0, 0])

    >>> dup_terms_gcd(f, ZZ)
    (2, [1, 0, 1])

    """
    if dup_TC(f, K) or not f:
        return 0, f

    i = 0

    for c in reversed(f):
        if not c:
            i += 1
        else:
            break

    return i, f[:-i]


def dmp_terms_gcd(f: dmp[Er], u: int, K: Domain[Er]) -> tuple[tuple[int, ...], dmp[Er]]:
    """
    Remove GCD of terms from ``f`` in ``K[X]``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_terms_gcd

    >>> f = ZZ.map([[1, 0], [1, 0, 0], [], []])

    >>> dmp_terms_gcd(f, 1, ZZ)
    ((2, 1), [[1], [1, 0]])

    """
    if dmp_ground_TC(f, u, K) or dmp_zero_p(f, u):
        return (0,) * (u + 1), f

    F: dict[monom, Er] = dmp_to_dict(f, u)
    G = monomial_min(*list(F.keys()))

    if all(g == 0 for g in G):
        return G, f

    d: dict[monom, Er] = {}

    for monom, coeff in F.items():
        d[monomial_ldiv(monom, G)] = coeff

    return G, dmp_from_dict(d, u, K)


def _rec_list_terms(g: dmp[Er], v: int, monom: monom) -> list[tuple[monom, Er]]:
    """Recursive helper for :func:`dmp_list_terms`."""
    d: int = dmp_degree(g, v)  # type: ignore
    terms: list[tuple[tuple[int, ...], Er]] = []

    if not v:
        c: Er
        for i, c in enumerate(_dup(g)):
            if not c:
                continue

            terms.append((monom + (d - i,), c))
    else:
        w = v - 1

        for i, cp in enumerate(g):
            terms.extend(_rec_list_terms(cp, w, monom + (d - i,)))

    return terms


def dmp_list_terms(
    f: dmp[Er], u: int, K: Domain[Er], order: MonomialOrder | str | None = None
) -> list[tuple[monom, Er]]:
    """
    List all non-zero terms from ``f`` in the given order ``order``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_list_terms

    >>> f = ZZ.map([[1, 1], [2, 3]])

    >>> dmp_list_terms(f, 1, ZZ)
    [((1, 1), 1), ((1, 0), 1), ((0, 1), 2), ((0, 0), 3)]
    >>> dmp_list_terms(f, 1, ZZ, order='grevlex')
    [((1, 1), 1), ((1, 0), 1), ((0, 1), 2), ((0, 0), 3)]

    """
    terms: list[tuple[monom, Er]]

    def sort(terms, O):
        return sorted(terms, key=lambda term: O(term[0]), reverse=True)

    terms = _rec_list_terms(f, u, ())

    if not terms:
        return [((0,) * (u + 1), K.zero)]

    if order is None:
        return terms
    else:
        return sort(terms, monomial_key(order))


def dup_apply_pairs(
    f: dup[Er], g: dup[Er], h: Callable[..., Er], args: tuple[Any, ...], K: Domain[Er]
) -> dup[Er]:
    """
    Apply ``h`` to pairs of coefficients of ``f`` and ``g``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dup_apply_pairs

    >>> h = lambda x, y, z: 2*x + y - z

    >>> dup_apply_pairs([1, 2, 3], [3, 2, 1], h, (1,), ZZ)
    [4, 5, 6]

    """
    n, m = len(f), len(g)

    if n != m:
        if n > m:
            g = [K.zero] * (n - m) + g
        else:
            f = [K.zero] * (m - n) + f

    result: dup[Er] = []

    for a, b in zip(f, g):
        result.append(h(a, b, *args))

    return dup_strip(result)


def dmp_apply_pairs(
    f: dmp[Er],
    g: dmp[Er],
    h: Callable[..., Er],
    args: tuple[Any, ...],
    u: int,
    K: Domain[Er],
) -> dmp[Er]:
    """
    Apply ``h`` to pairs of coefficients of ``f`` and ``g``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dmp_apply_pairs

    >>> h = lambda x, y, z: 2*x + y - z

    >>> dmp_apply_pairs([[1], [2, 3]], [[3], [2, 1]], h, (1,), 1, ZZ)
    [[4], [5, 6]]

    """
    if not u:
        return _dmp(dup_apply_pairs(_dup(f), _dup(g), h, args, K))

    n, m, v = len(f), len(g), u - 1

    if n != m:
        if n > m:
            g = dmp_zeros(n - m, v, K) + g
        else:
            f = dmp_zeros(m - n, v, K) + f

    result = []

    for a, b in zip(f, g):
        result.append(dmp_apply_pairs(a, b, h, args, v, K))

    return dmp_strip(result, u)


def dup_slice(f: dup[Er], m: int, n: int, K: Domain[Er]) -> dup[Er]:
    """Take a continuous subsequence of terms of ``f`` in ``K[x]``."""
    k = len(f)

    if k >= m:
        M = k - m
    else:
        M = 0
    if k >= n:
        N = k - n
    else:
        N = 0

    f = f[N:M]

    while f and f[0] == K.zero:
        f.pop(0)

    if not f:
        return []
    else:
        return f + [K.zero] * m


def dmp_slice(f: dmp[Er], m: int, n: int, u: int, K: Domain[Er]) -> dmp[Er]:
    """Take a continuous subsequence of terms of ``f`` in ``K[X]``."""
    return dmp_slice_in(f, m, n, 0, u, K)


def dup_truncate(f, n, K):
    """
    Truncate ``f`` to the first ``n`` terms in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dup_truncate
    >>> from sympy.polys.densebasic import dup_from_list, dup_print
    >>> f = dup_from_list([1, 2, 3, 4, 5], ZZ)
    >>> t = dup_truncate(f, 3, ZZ)
    >>> dup_print(t, 'x')
    3*x**2 + 4*x + 5

    """
    return dup_slice(f, 0, n, K)


def dmp_slice_in(f: dmp[Er], m: int, n: int, j: int, u: int, K: Domain[Er]) -> dmp[Er]:
    """Take a continuous subsequence of terms of ``f`` in ``x_j`` in ``K[X]``."""
    if j < 0 or j > u:
        raise IndexError("-%s <= j < %s expected, got %s" % (u, u, j))

    if not u:
        return _dmp(dup_slice(_dup(f), m, n, K))

    d: dict[tuple[int, ...], Er] = dmp_to_dict(f, u)
    g: dict[tuple[int, ...], Er] = {}

    for monom, coeff in d.items():
        k = monom[j]

        if k < m or k >= n:
            monom = monom[:j] + (0,) + monom[j + 1 :]

        if monom in g:
            g[monom] += coeff
        else:
            g[monom] = coeff

    return dmp_from_dict(g, u, K)


def dup_random(n: int, a: int, b: int, K: Domain[Er]) -> dup[Er]:
    """
    Return a polynomial of degree ``n`` with coefficients in ``[a, b]``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dup_random

    >>> dup_random(3, -10, 10, ZZ) #doctest: +SKIP
    [-2, -8, 9, -4]

    """
    f = [K.convert(random.randint(a, b)) for _ in range(0, n + 1)]

    while not f[0]:
        f[0] = K.convert(random.randint(a, b))

    return f


def dup_from_list(f: list[Er], K: Domain[Er]) -> dup[Er]:
    """
    Create a ``K[x]`` polynomial from a list.

    Examples
    ========

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.densebasic import dup_from_list

    >>> p = dup_from_list([1, 0, 5, 0, 7], QQ); p
    [1, 0, 5, 0, 7]
    >>> QQ.of_type(1)
    False
    >>> QQ.of_type(p[0])
    True

    >>> dup_from_list([0, 0, 0, 2, 1, 0, 0], QQ)
    [2, 1, 0, 0]

    """
    if not f:
        return []

    return dup_strip([K.convert(c) for c in f])


def dup_print(f: dup[Er], sym: str) -> None:
    """
    Print a polynomial in ``K[x]``.

    Examples
    ========

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densebasic import dup_print, dup_from_list
    >>> f = dup_from_list([1, 0, 5, 0, 7], ZZ)
    >>> dup_print(f, 'x')
    x**4 + 5*x**2 + 7

    """
    if isinstance(sym, tuple):
        sym = sym[0]

    if not f:
        print("0")
        return

    deg = dup_degree(f)
    terms = []

    for i, coeff in enumerate(f):
        d = deg - i

        if coeff == 0:
            continue

        if d == 0:
            term = f"{coeff}"
        elif d == 1:
            if coeff == 1:
                term = f"{sym}"
            elif coeff == -1:
                term = f"-{sym}"
            else:
                term = f"{coeff}*{sym}"
        else:
            if coeff == 1:
                term = f"{sym}**{d}"
            elif coeff == -1:
                term = f"-{sym}**{d}"
            else:
                term = f"{coeff}*{sym}**{d}"

        terms.append(term)

    poly = " + ".join(terms).replace("+ -", "- ")
    print(poly)

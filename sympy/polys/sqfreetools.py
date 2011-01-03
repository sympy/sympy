"""Square--free decomposition algorithms and related tools. """

from sympy.polys.densebasic import (
    dup_strip,
    dup_LC, dmp_ground_LC,
    dmp_zero_p,
    dmp_ground,
    dup_degree, dmp_degree,
    dmp_raise, dmp_inject,
    dup_convert)

from sympy.polys.densearith import (
    dup_neg, dmp_neg,
    dup_sub, dmp_sub,
    dup_mul, dmp_mul,
    dup_exquo, dmp_exquo,
    dup_mul_ground, dmp_mul_ground)

from sympy.polys.densetools import (
    dup_diff, dmp_diff,
    dup_shift, dmp_compose,
    dup_monic, dmp_ground_monic,
    dup_primitive, dmp_ground_primitive)

from sympy.polys.euclidtools import (
    dup_inner_gcd, dmp_inner_gcd,
    dup_gcd, dmp_gcd,
    dmp_resultant)

from sympy.polys.galoistools import (
    gf_sqf_list, gf_sqf_part)

from sympy.polys.polyerrors import (
    MultivariatePolynomialError,
    DomainError)

from sympy.utilities import cythonized

def dup_sqf_p(f, K):
    """
    Return ``True`` if ``f`` is a square--free polynomial in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.sqfreetools import dup_sqf_p

    >>> dup_sqf_p([ZZ(1),-ZZ(2), ZZ(1)], ZZ)
    False
    >>> dup_sqf_p([ZZ(1), ZZ(0),-ZZ(1)], ZZ)
    True

    """
    if not f:
        return True
    else:
        return not dup_degree(dup_gcd(f, dup_diff(f, 1, K), K))

@cythonized("u")
def dmp_sqf_p(f, u, K):
    """
    Return ``True`` if ``f`` is a square--free polynomial in ``K[X]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.sqfreetools import dmp_sqf_p

    >>> f = ZZ.map([[1], [2, 0], [1, 0, 0]])
    >>> dmp_sqf_p(f, 1, ZZ)
    False

    >>> f = ZZ.map([[1], [], [1, 0, 0]])
    >>> dmp_sqf_p(f, 1, ZZ)
    True

    """
    if dmp_zero_p(f, u):
        return True
    else:
        return not dmp_degree(dmp_gcd(f, dmp_diff(f, 1, u, K), u, K), u)

@cythonized("s")
def dup_sqf_norm(f, K):
    """
    Square--free norm of ``f`` in ``K[x]``, useful over algebraic domains.

    Returns ``s``, ``f``, ``r``, such that ``g(x) = f(x-sa)`` and ``r(x) = Norm(g(x))``
    is a square-free polynomial over K, where ``a`` is the algebraic extension of ``K``.

    Example
    =======

    >>> from sympy import sqrt
    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.sqfreetools import dup_sqf_norm

    >>> K = QQ.algebraic_field(sqrt(3))

    >>> s, f, r = dup_sqf_norm([K(1), K(0), K(-2)], K)

    >>> s == 1
    True
    >>> f == [K(1), K([QQ(-2), QQ(0)]), K(1)]
    True
    >>> r == [1, 0, -10, 0, 1]
    True

    """
    if not K.is_Algebraic:
        raise DomainError("ground domain must be algebraic")

    s, g = 0, dmp_raise(K.mod.rep, 1, 0, K.dom)

    while True:
        h, _ = dmp_inject(f, 0, K, front=True)
        r = dmp_resultant(g, h, 1, K.dom)

        if dup_sqf_p(r, K.dom):
            break
        else:
            f, s = dup_shift(f, -K.unit, K), s+1

    return s, f, r

@cythonized("s,u")
def dmp_sqf_norm(f, u, K):
    """
    Square-free norm of ``f`` in ``K[X]``, useful over algebraic domains.

    Returns ``s``, ``f``, ``r``, such that ``g(x) = f(x-sa)`` and ``r(x) = Norm(g(x))``
    is a square--free polynomial over K, where ``a`` is the algebraic extension of ``K``.

    Example
    =======

    >>> from sympy import I
    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.sqfreetools import dmp_sqf_norm

    >>> K = QQ.algebraic_field(I)

    >>> s, f, r = dmp_sqf_norm([[K(1), K(0)], [K(1), K(0), K(0)]], 1, K)

    >>> s == 1
    True
    >>> f == [[K(1), K(0)], [K(1), K([QQ(-1), QQ(0)]), K(0)]]
    True
    >>> r == [[1, 0, 0], [2, 0, 0, 0], [1, 0, 1, 0, 0]]
    True

    """
    if not u:
        return dup_sqf_norm(f, K)

    if not K.is_Algebraic:
        raise DomainError("ground domain must be algebraic")

    g = dmp_raise(K.mod.rep, u+1, 0, K.dom)
    F = dmp_raise([K.one,-K.unit], u, 0, K)

    s = 0

    while True:
        h, _ = dmp_inject(f, u, K, front=True)
        r = dmp_resultant(g, h, u+1, K.dom)

        if dmp_sqf_p(r, u, K.dom):
            break
        else:
            f, s = dmp_compose(f, F, u, K), s+1

    return s, f, r

@cythonized("i")
def dup_gf_sqf_part(f, K):
    """Compute square--free part of ``f`` in ``GF(p)[x]``. """
    f = dup_convert(f, K, K.dom)
    g = gf_sqf_part(f, K.mod, K.dom)
    return dup_convert(g, K.dom, K)

def dmp_gf_sqf_part(f, K):
    """Compute square--free part of ``f`` in ``GF(p)[X]``. """
    raise DomainError('multivariate polynomials over %s' % K)

def dup_sqf_part(f, K):
    """
    Returns square--free part of a polynomial in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.sqfreetools import dup_sqf_part

    >>> dup_sqf_part([ZZ(1), ZZ(0), -ZZ(3), -ZZ(2)], ZZ)
    [1, -1, -2]

    """
    if not K.has_CharacteristicZero:
        return dup_gf_sqf_part(f, K)

    if not f:
        return f

    if K.is_negative(dup_LC(f, K)):
        f = dup_neg(f, K)

    gcd = dup_gcd(f, dup_diff(f, 1, K), K)
    sqf = dup_exquo(f, gcd, K)

    if K.has_Field or not K.is_Exact:
        return dup_monic(sqf, K)
    else:
        return dup_primitive(sqf, K)[1]

@cythonized("u")
def dmp_sqf_part(f, u, K):
    """
    Returns square--free part of a polynomial in ``K[X]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.sqfreetools import dmp_sqf_part

    >>> f = ZZ.map([[1], [2, 0], [1, 0, 0], []])

    >>> dmp_sqf_part(f, 1, ZZ)
    [[1], [1, 0], []]

    """
    if not u:
        return dup_sqf_part(f, K)

    if not K.has_CharacteristicZero:
        return dmp_gf_sqf_part(f, u, K)

    if dmp_zero_p(f, u):
        return f

    if K.is_negative(dmp_ground_LC(f, u, K)):
        f = dmp_neg(f, u, K)

    gcd = dmp_gcd(f, dmp_diff(f, 1, u, K), u, K)
    sqf = dmp_exquo(f, gcd, u, K)

    if K.has_Field or not K.is_Exact:
        return dmp_ground_monic(sqf, u, K)
    else:
        return dmp_ground_primitive(sqf, u, K)[1]

@cythonized("i")
def dup_gf_sqf_list(f, K, all=False):
    """Compute square--free decomposition of ``f`` in ``GF(p)[x]``. """
    f = dup_convert(f, K, K.dom)

    coeff, factors = gf_sqf_list(f, K.mod, K.dom, all=all)

    for i, (f, k) in enumerate(factors):
        factors[i] = (dup_convert(f, K.dom, K), k)

    return K.convert(coeff, K.dom), factors

def dmp_gf_sqf_list(f, u, K, all=False):
    """Compute square--free decomposition of ``f`` in ``GF(p)[X]``. """
    raise DomainError('multivariate polynomials over %s' % K)

@cythonized("i")
def dup_sqf_list(f, K, all=False):
    """
    Return square--free decomposition of a polynomial in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.sqfreetools import dup_sqf_list

    >>> f = ZZ.map([2, 16, 50, 76, 56, 16])

    >>> dup_sqf_list(f, ZZ)
    (2, [([1, 1], 2), ([1, 2], 3)])

    >>> dup_sqf_list(f, ZZ, all=True)
    (2, [([1], 1), ([1, 1], 2), ([1, 2], 3)])

    """
    if not K.has_CharacteristicZero:
        return dup_gf_sqf_list(f, K, all=all)

    if K.has_Field or not K.is_Exact:
        coeff = dup_LC(f, K)
        f = dup_monic(f, K)
    else:
        coeff, f = dup_primitive(f, K)

        if K.is_negative(dup_LC(f, K)):
            f = dup_neg(f, K)
            coeff = -coeff

    if dup_degree(f) <= 0:
        return coeff, []

    result, i = [], 1

    h = dup_diff(f, 1, K)
    g, p, q = dup_inner_gcd(f, h, K)

    while True:
        d = dup_diff(p, 1, K)
        h = dup_sub(q, d, K)

        if not h:
            result.append((p, i))
            break

        g, p, q = dup_inner_gcd(p, h, K)

        if all or dup_degree(g) > 0:
            result.append((g, i))

        i += 1

    return coeff, result

def dup_sqf_list_include(f, K, all=False):
    """
    Return square--free decomposition of a polynomial in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.sqfreetools import dup_sqf_list_include

    >>> f = ZZ.map([2, 16, 50, 76, 56, 16])

    >>> dup_sqf_list_include(f, ZZ)
    [([2, 2], 2), ([1, 2], 3)]

    >>> dup_sqf_list_include(f, ZZ, all=True)
    [([2], 1), ([1, 1], 2), ([1, 2], 3)]

    """
    coeff, factors = dup_sqf_list(f, K, all=all)

    if not factors:
        return [(dup_strip([coeff]), 1)]
    else:
        g = dup_mul_ground(factors[0][0], coeff, K)
        return [(g, factors[0][1])] + factors[1:]

@cythonized("u,i")
def dmp_sqf_list(f, u, K, all=False):
    """
    Return square--free decomposition of a polynomial in ``K[X]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.sqfreetools import dmp_sqf_list

    >>> f = ZZ.map([[1], [2, 0], [1, 0, 0], [], [], []])

    >>> dmp_sqf_list(f, 1, ZZ)
    (1, [([[1], [1, 0]], 2), ([[1], []], 3)])

    >>> dmp_sqf_list(f, 1, ZZ, all=True)
    (1, [([[1]], 1), ([[1], [1, 0]], 2), ([[1], []], 3)])

    """
    if not u:
        return dup_sqf_list(f, K, all=all)

    if not K.has_CharacteristicZero:
        return dmp_gf_sqf_list(f, u, K, all=all)

    if K.has_Field or not K.is_Exact:
        coeff = dmp_ground_LC(f, u, K)
        f = dmp_ground_monic(f, u, K)
    else:
        coeff, f = dmp_ground_primitive(f, u, K)

        if K.is_negative(dmp_ground_LC(f, u, K)):
            f = dmp_neg(f, u, K)
            coeff = -coeff

    if dmp_degree(f, u) <= 0:
        return coeff, []

    result, i = [], 1

    h = dmp_diff(f, 1, u, K)
    g, p, q = dmp_inner_gcd(f, h, u, K)

    while True:
        d = dmp_diff(p, 1, u, K)
        h = dmp_sub(q, d, u, K)

        if dmp_zero_p(h, u):
            result.append((p, i))
            break

        g, p, q = dmp_inner_gcd(p, h, u, K)

        if all or dmp_degree(g, u) > 0:
            result.append((g, i))

        i += 1

    return coeff, result

@cythonized("u")
def dmp_sqf_list_include(f, u, K, all=False):
    """
    Return square--free decomposition of a polynomial in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.sqfreetools import dmp_sqf_list_include

    >>> f = ZZ.map([[1], [2, 0], [1, 0, 0], [], [], []])

    >>> dmp_sqf_list_include(f, 1, ZZ)
    [([[1], [1, 0]], 2), ([[1], []], 3)]

    >>> dmp_sqf_list_include(f, 1, ZZ, all=True)
    [([[1]], 1), ([[1], [1, 0]], 2), ([[1], []], 3)]

    """
    if not u:
        return dup_sqf_list_include(f, K, all=all)

    coeff, factors = dmp_sqf_list(f, u, K, all=all)

    if not factors:
        return [(dmp_ground(coeff, u), 1)]
    else:
        g = dmp_mul_ground(factors[0][0], coeff, u, K)
        return [(g, factors[0][1])] + factors[1:]

def dup_gff_list(f, K):
    """
    Compute greatest factorial factorization of ``f`` in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.sqfreetools import dup_gff_list

    >>> f = ZZ.map([1, 2, -1, -2, 0, 0])

    >>> dup_gff_list(f, ZZ)
    [([1, 0], 1), ([1, 2], 4)]

    """
    if not f:
        raise ValueError("greatest factorial factorization doesn't exist for a zero polynomial")

    f = dup_monic(f, K)

    if not dup_degree(f):
        return []
    else:
        g = dup_gcd(f, dup_shift(f, K.one, K), K)
        H = dup_gff_list(g, K)

        for i, (h, k) in enumerate(H):
            g = dup_mul(g, dup_shift(h, -K(k), K), K)
            H[i] = (h, k + 1)

        f = dup_exquo(f, g, K)

        if not dup_degree(f):
            return H
        else:
            return [(f, 1)] + H

def dmp_gff_list(f, u, K):
    """
    Compute greatest factorial factorization of ``f`` in ``K[X]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.sqfreetools import dmp_gff_list

    """
    if not u:
        return dup_gff_list(f, K)
    else:
        raise MultivariatePolynomialError(f)

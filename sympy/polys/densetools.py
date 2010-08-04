"""Advanced tools for dense recursive polynomials in ``K[x]`` or ``K[X]``. """

from sympy.polys.densebasic import (
    dup_strip, dmp_strip,
    dup_convert, dmp_convert,
    dup_degree, dmp_degree, dmp_degree_in,
    dup_to_dict, dmp_to_dict,
    dup_from_dict, dmp_from_dict,
    dup_LC, dmp_LC, dmp_ground_LC,
    dup_TC, dmp_TC, dmp_ground_TC,
    dmp_zero, dmp_one, dmp_ground,
    dmp_zero_p, dmp_one_p,
    dmp_multi_deflate, dmp_inflate,
    dup_to_raw_dict, dup_from_raw_dict,
    dmp_raise, dmp_apply_pairs,
    dmp_inject, dmp_zeros,
    dup_terms_gcd
)

from sympy.polys.densearith import (
    dup_add_term, dmp_add_term,
    dup_mul_term, dmp_mul_term,
    dup_lshift, dup_rshift,
    dup_neg, dmp_neg,
    dup_add, dmp_add,
    dup_sub, dmp_sub,
    dup_mul, dmp_mul,
    dup_sqr, dmp_sqr,
    dup_pow, dmp_pow,
    dup_div, dmp_div,
    dup_rem, dmp_rem,
    dup_quo, dmp_quo,
    dup_exquo, dmp_exquo,
    dup_prem, dmp_prem,
    dup_expand, dmp_expand,
    dup_add_mul, dup_sub_mul,
    dup_mul_ground, dmp_mul_ground,
    dup_quo_ground, dmp_quo_ground,
    dup_exquo_ground, dmp_exquo_ground,
    dup_max_norm, dmp_max_norm
)

from sympy.polys.polyerrors import (
    HeuristicGCDFailed,
    HomomorphismFailed,
    RefinementFailed,
    NotInvertible,
    DomainError
)

from sympy.utilities import (
    cythonized, variations
)

from math import ceil, log

@cythonized("m,n,i,j")
def dup_integrate(f, m, K):
    """
    Computes indefinite integral of ``f`` in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.densetools import dup_integrate

    >>> dup_integrate([QQ(1), QQ(2), QQ(0)], 1, QQ)
    [1/3, 1/1, 0/1, 0/1]
    >>> dup_integrate([QQ(1), QQ(2), QQ(0)], 2, QQ)
    [1/12, 1/3, 0/1, 0/1, 0/1]

    """
    if m <= 0 or not f:
        return f

    g = [K.zero]*m

    for i, c in enumerate(reversed(f)):
        n = i+1

        for j in xrange(1, m):
            n *= i+j+1

        g.insert(0, K.exquo(c, K(n)))

    return g

@cythonized("m,u,v,n,i,j")
def dmp_integrate(f, m, u, K):
    """
    Computes indefinite integral of ``f`` in ``x_0`` in ``K[X]``.

    Example
    =======

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.densetools import dmp_integrate

    >>> dmp_integrate([[QQ(1)], [QQ(2), QQ(0)]], 1, 1, QQ)
    [[1/2], [2/1, 0/1], []]
    >>> dmp_integrate([[QQ(1)], [QQ(2), QQ(0)]], 2, 1, QQ)
    [[1/6], [1/1, 0/1], [], []]

    """
    if not u:
        return dup_integrate(f, m, K)

    if m <= 0 or dmp_zero_p(f, u):
        return f

    g, v = dmp_zeros(m, u-1, K), u-1

    for i, c in enumerate(reversed(f)):
        n = i+1

        for j in xrange(1, m):
            n *= i+j+1

        g.insert(0, dmp_quo_ground(c, K(n), v, K))

    return g

@cythonized("m,v,w,i,j")
def _rec_integrate_in(g, m, v, i, j, K):
    """Recursive helper for :func:`dmp_integrate_in`."""
    if i == j:
        return dmp_integrate(g, m, v, K)

    w, i = v-1, i+1

    return dmp_strip([ _rec_integrate_in(c, m, w, i, j, K) for c in g ], v)

@cythonized("m,j,u")
def dmp_integrate_in(f, m, j, u, K):
    """
    Computes indefinite integral of ``f`` in ``x_j`` in ``K[X]``.

    Example
    =======

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.densetools import dmp_integrate_in

    >>> dmp_integrate_in([[QQ(1)], [QQ(2), QQ(0)]], 1, 0, 1, QQ)
    [[1/2], [2/1, 0/1], []]
    >>> dmp_integrate_in([[QQ(1)], [QQ(2), QQ(0)]], 1, 1, 1, QQ)
    [[1/1, 0/1], [1/1, 0/1, 0/1]]

    """
    if j < 0 or j > u:
        raise IndexError("-%s <= j < %s expected, got %s" % (u, u, j))

    return _rec_integrate_in(f, m, u, 0, j, K)

@cythonized("m,n,k,i")
def dup_diff(f, m, K):
    """
    ``m``--th order derivative of a polynomial in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dup_diff

    >>> dup_diff([ZZ(1), ZZ(2), ZZ(3), ZZ(4)], 1, ZZ)
    [3, 4, 3]
    >>> dup_diff([ZZ(1), ZZ(2), ZZ(3), ZZ(4)], 2, ZZ)
    [6, 4]

    """
    if m <= 0:
        return f

    n = dup_degree(f)

    if n < m:
        return []

    deriv = []

    if m == 1:
        for coeff in f[:-m]:
            deriv.append(K(n)*coeff)
            n -= 1
    else:
        for coeff in f[:-m]:
            k = n

            for i in xrange(n-1, n-m, -1):
                k *= i

            deriv.append(K(k)*coeff)
            n -= 1

    return dup_strip(deriv)

@cythonized("u,v,m,n,k,i")
def dmp_diff(f, m, u, K):
    """
    ``m``-th order derivative in ``x_0`` of a polynomial in ``K[X]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dmp_diff

    >>> f = ZZ.map([[1, 2, 3], [2, 3, 1]])

    >>> dmp_diff(f, 1, 1, ZZ)
    [[1, 2, 3]]
    >>> dmp_diff(f, 2, 1, ZZ)
    [[]]

    """
    if not u:
        return dup_diff(f, m, K)
    if m <= 0:
        return f

    n = dmp_degree(f, u)

    if n < m:
        return dmp_zero(u)

    deriv, v = [], u-1

    if m == 1:
        for coeff in f[:-m]:
            deriv.append(dmp_mul_ground(coeff, K(n), v, K))
            n -= 1
    else:
        for coeff in f[:-m]:
            k = n

            for i in xrange(n-1, n-m, -1):
                k *= i

            deriv.append(dmp_mul_ground(coeff, K(k), v, K))
            n -= 1

    return dmp_strip(deriv, u)

@cythonized("m,v,w,i,j")
def _rec_diff_in(g, m, v, i, j, K):
    """Recursive helper for :func:`dmp_diff_in`."""
    if i == j:
        return dmp_diff(g, m, v, K)

    w, i = v-1, i+1

    return dmp_strip([ _rec_diff_in(c, m, w, i, j, K) for c in g ], v)

@cythonized("m,j,u")
def dmp_diff_in(f, m, j, u, K):
    """
    ``m``--th order derivative in ``x_j`` of a polynomial in ``K[X]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dmp_diff_in

    >>> f = ZZ.map([[1, 2, 3], [2, 3, 1]])

    >>> dmp_diff_in(f, 1, 0, 1, ZZ)
    [[1, 2, 3]]
    >>> dmp_diff_in(f, 1, 1, 1, ZZ)
    [[2, 2], [4, 3]]

    """
    if j < 0 or j > u:
        raise IndexError("-%s <= j < %s expected, got %s" % (u, u, j))

    return _rec_diff_in(f, m, u, 0, j, K)

def dup_eval(f, a, K):
    """
    Evaluate a polynomial at ``x = a`` in ``K[x]`` using Horner scheme.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dup_eval

    >>> dup_eval([ZZ(1), ZZ(2), ZZ(3)], 2, ZZ)
    11

    """
    if not a:
        return dup_TC(f, K)

    result = K.zero

    for c in f:
        result *= a
        result += c

    return result

@cythonized("u,v")
def dmp_eval(f, a, u, K):
    """
    Evaluate a polynomial at ``x_0 = a`` in ``K[X]`` using the Horner scheme.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dmp_eval

    >>> f = ZZ.map([[2, 3], [1, 2]])

    >>> dmp_eval(f, 2, 1, ZZ)
    [5, 8]

    """
    if not u:
        return dup_eval(f, a, K)

    if not a:
        return dmp_TC(f, K)

    result, v = dmp_LC(f, K), u-1

    for coeff in f[1:]:
        result = dmp_mul_ground(result, a, v, K)
        result = dmp_add(result, coeff, v, K)

    return result

@cythonized("v,i,j")
def _rec_eval_in(g, a, v, i, j, K):
    """Recursive helper for :func:`dmp_eval_in`."""
    if i == j:
        return dmp_eval(g, a, v, K)

    v, i = v-1, i+1

    return dmp_strip([ _rec_eval_in(c, a, v, i, j, K) for c in g ], v)

@cythonized("u")
def dmp_eval_in(f, a, j, u, K):
    """
    Evaluate a polynomial at ``x_j = a`` in ``K[X]`` using the Horner scheme.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dmp_eval_in

    >>> f = ZZ.map([[2, 3], [1, 2]])

    >>> dmp_eval_in(f, 2, 0, 1, ZZ)
    [5, 8]
    >>> dmp_eval_in(f, 2, 1, 1, ZZ)
    [7, 4]

    """
    if j < 0 or j > u:
        raise IndexError("-%s <= j < %s expected, got %s" % (u, u, j))

    return _rec_eval_in(f, a, u, 0, j, K)

@cythonized("i,u")
def _rec_eval_tail(g, i, A, u, K):
    """Recursive helper for :func:`dmp_eval_tail`."""
    if i == u:
        return dup_eval(g, A[-1], K)
    else:
        h = [ _rec_eval_tail(c, i+1, A, u, K) for c in g ]

        if i < u - len(A) + 1:
            return h
        else:
            return dup_eval(h, A[-u+i-1], K)

@cythonized("u")
def dmp_eval_tail(f, A, u, K):
    """
    Evaluate a polynomial at ``x_j = a_j, ...`` in ``K[X]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dmp_eval_tail

    >>> f = ZZ.map([[2, 3], [1, 2]])

    >>> dmp_eval_tail(f, (2, 2), 1, ZZ)
    18
    >>> dmp_eval_tail(f, (2,), 1, ZZ)
    [7, 4]

    """
    if not A:
        return f

    if dmp_zero_p(f, u):
        return dmp_zero(u - len(A))

    e = _rec_eval_tail(f, 0, A, u, K)

    if u == len(A)-1:
        return e
    else:
        return dmp_strip(e, u - len(A))

@cythonized("m,v,i,j")
def _rec_diff_eval(g, m, a, v, i, j, K):
    """Recursive helper for :func:`dmp_diff_eval`."""
    if i == j:
        return dmp_eval(dmp_diff(g, m, v, K), a, v, K)

    v, i = v-1, i+1

    return dmp_strip([ _rec_diff_eval(c, m, a, v, i, j, K) for c in g ], v)

@cythonized("m,j,u")
def dmp_diff_eval_in(f, m, a, j, u, K):
    """
    Differentiate and evaluate a polynomial in ``x_j`` at ``a`` in ``K[X]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dmp_diff_eval_in

    >>> f = ZZ.map([[1, 2, 3], [2, 3, 1]])

    >>> dmp_diff_eval_in(f, 1, 2, 0, 1, ZZ)
    [1, 2, 3]
    >>> dmp_diff_eval_in(f, 1, 2, 1, 1, ZZ)
    [6, 11]

    """
    if j > u:
        raise IndexError("-%s <= j < %s expected, got %s" % (u, u, j))
    if not j:
        return dmp_eval(dmp_diff(f, m, u, K), a, u, K)

    return _rec_diff_eval(f, m, a, u, 0, j, K)

def dup_trunc(f, p, K):
    """
    Reduce ``K[x]`` polynomial modulo a constant ``p`` in ``K``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dup_trunc

    >>> f = ZZ.map([2, 3, 5, 7])

    >>> dup_trunc(f, ZZ(3), ZZ)
    [-1, 0, -1, 1]

    """
    if K.is_ZZ:
        g = []

        for c in f:
            c = c % p

            if c > p // 2:
                g.append(c - p)
            else:
                g.append(c)
    else:
        g = [ c % p for c in f ]

    return dup_strip(g)

@cythonized("u")
def dmp_trunc(f, p, u, K):
    """
    Reduce ``K[X]`` polynomial modulo a polynomial ``p`` in ``K[Y]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dmp_trunc

    >>> f = ZZ.map([[3, 8], [5, 6], [2, 3]])
    >>> g = ZZ.map([1, -1])

    >>> dmp_trunc(f, g, 1, ZZ)
    [[11], [11], [5]]

    """
    return dmp_strip([ dmp_rem(c, p, u-1, K) for c in f ], u)

@cythonized("u,v")
def dmp_ground_trunc(f, p, u, K):
    """
    Reduce ``K[X]`` polynomial modulo a constant ``p`` in ``K``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dmp_ground_trunc

    >>> f = ZZ.map([[3, 8], [5, 6], [2, 3]])

    >>> dmp_ground_trunc(f, ZZ(3), 1, ZZ)
    [[-1], [-1, 0], [-1, 0]]

    """
    if not u:
        return dup_trunc(f, p, K)

    v = u-1

    return dmp_strip([ dmp_ground_trunc(c, p, v, K) for c in f ], u)

def dup_monic(f, K):
    """
    Divides all coefficients by ``LC(f)`` in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ, QQ
    >>> from sympy.polys.densetools import dup_monic

    >>> dup_monic([ZZ(3), ZZ(6), ZZ(9)], ZZ)
    [1, 2, 3]

    >>> dup_monic([QQ(3), QQ(4), QQ(2)], QQ)
    [1/1, 4/3, 2/3]

    """
    if not f:
        return f

    lc = dup_LC(f, K)

    if K.is_one(lc):
        return f
    else:
        return dup_quo_ground(f, lc, K)

@cythonized("u")
def dmp_ground_monic(f, u, K):
    """
    Divides all coefficients by ``LC(f)`` in ``K[X]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ, QQ
    >>> from sympy.polys.densetools import dmp_ground_monic

    >>> f = ZZ.map([[3, 6], [3, 0], [9, 3]])
    >>> g = QQ.map([[3, 8], [5, 6], [2, 3]])

    >>> dmp_ground_monic(f, 1, ZZ)
    [[1, 2], [1, 0], [3, 1]]

    >>> dmp_ground_monic(g, 1, QQ)
    [[1/1, 8/3], [5/3, 2/1], [2/3, 1/1]]

    """
    if not u:
        return dup_monic(f, K)

    if dmp_zero_p(f, u):
        return f

    lc = dmp_ground_LC(f, u, K)

    if K.is_one(lc):
        return f
    else:
        return dmp_quo_ground(f, lc, u, K)

def dup_rr_content(f, K):
    """
    Returns GCD of coefficients over a ring.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dup_rr_content

    >>> dup_rr_content([ZZ(6), ZZ(8), ZZ(12)], ZZ)
    2

    """
    cont = K.zero

    for c in f:
        cont = K.gcd(cont, c)

        if K.is_one(cont):
            break

    return cont

def dup_ff_content(f, K):
    """
    Returns GCD of coefficients over a field.

    Example
    =======

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.densetools import dup_ff_content

    >>> dup_ff_content([], QQ)
    0/1
    >>> dup_ff_content([QQ(1), QQ(1,2), QQ(2,3)], QQ)
    1/1

    """
    if not f:
        return K.zero
    else:
        return K.one

def dup_content(f, K):
    """
    Returns GCD of coefficients in `K[x]`.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dup_content

    >>> dup_content([ZZ(6), ZZ(8), ZZ(12)], ZZ)
    2

    """
    if K.has_Field or not K.is_Exact:
        return dup_ff_content(f, K)
    else:
        return dup_rr_content(f, K)

@cythonized("u,v")
def dmp_rr_ground_content(f, u, K):
    """
    Returns GCD of coefficients over a ring.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dmp_rr_ground_content

    >>> f = ZZ.map([[2, 6], [4, 12]])

    >>> dmp_rr_ground_content(f, 1, ZZ)
    2

    """
    if not u:
        return dup_rr_content(f, K)

    cont, v = K.zero, u-1

    for c in f:
        gc = dmp_rr_ground_content(c, v, K)
        cont = K.gcd(cont, gc)

        if K.is_one(cont):
            break

    return cont

@cythonized("u")
def dmp_ff_ground_content(f, u, K):
    """
    Returns GCD of coefficients over a field.

    Example
    =======

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.densetools import dmp_ff_ground_content

    >>> dmp_ff_ground_content([[]], 1, QQ)
    0/1
    >>> dmp_ff_ground_content([[QQ(2), QQ(1,2)], [QQ(4)]], 1, QQ)
    1/1

    """
    if dmp_zero_p(f, u):
        return K.zero
    else:
        return K.one

@cythonized("u")
def dmp_ground_content(f, u, K):
    """
    Returns GCD of coefficients in ``K[X]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dmp_ground_content

    >>> f = ZZ.map([[2, 6], [4, 12]])

    >>> dmp_ground_content(f, 1, ZZ)
    2

    """
    if not u:
        return dup_content(f, K)

    if K.has_Field or not K.is_Exact:
        return dmp_ff_ground_content(f, u, K)
    else:
        return dmp_rr_ground_content(f, u, K)

def dup_rr_primitive(f, K):
    """
    Returns content and a primitive polynomial over a ring.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dup_rr_primitive

    >>> dup_rr_primitive([ZZ(6), ZZ(8), ZZ(12)], ZZ)
    (2, [3, 4, 6])

    """
    cont = dup_content(f, K)

    if not f or K.is_one(cont):
        return cont, f
    else:
        return cont, dup_exquo_ground(f, cont, K)

def dup_ff_primitive(f, K):
    """
    Returns content and a primitive polynomial over a field.

    Example
    =======

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.densetools import dup_ff_primitive

    >>> dup_ff_primitive([QQ(1), QQ(1,2), QQ(2,3)], QQ)
    (1/1, [1/1, 1/2, 2/3])

    """
    if not f:
        return K.zero, f
    else:
        return K.one, f

def dup_primitive(f, K):
    """
    Returns content and a primitive polynomial in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dup_primitive

    >>> dup_primitive([ZZ(2), ZZ(6), ZZ(12)], ZZ)
    (2, [1, 3, 6])

    """
    if K.has_Field or not K.is_Exact:
        return dup_ff_primitive(f, K)
    else:
        return dup_rr_primitive(f, K)

@cythonized("u")
def dmp_rr_ground_primitive(f, u, K):
    """
    Compute content and the primitive from of ``f`` over a ring.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dmp_rr_ground_primitive

    >>> f = ZZ.map([[2, 6], [4, 12]])

    >>> dmp_rr_ground_primitive(f, 1, ZZ)
    (2, [[1, 3], [2, 6]])

    """
    cont = dmp_ground_content(f, u, K)

    if K.is_one(cont):
        return cont, f
    else:
        return cont, dmp_exquo_ground(f, cont, u, K)

@cythonized("u")
def dmp_ff_ground_primitive(f, u, K):
    """
    Compute content and the primitive from of ``f`` over a field.

    Example
    =======

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.densetools import dmp_ff_ground_primitive

    >>> f = [[QQ(2), QQ(1,2)], [QQ(4)]]

    >>> dmp_ff_ground_primitive(f, 1, QQ)
    (1/1, [[2/1, 1/2], [4/1]])

    """
    if dmp_zero_p(f, u):
        return K.zero, f
    else:
        return K.one, f

@cythonized("u")
def dmp_ground_primitive(f, u, K):
    """
    Compute content and the primitive form of ``f`` in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dmp_ground_primitive

    >>> f = ZZ.map([[2, 6], [4, 12]])

    >>> dmp_ground_primitive(f, 1, ZZ)
    (2, [[1, 3], [2, 6]])

    """
    if not u:
        return dup_primitive(f, K)

    if dmp_zero_p(f, u):
        return K.zero, f

    if K.has_Field or not K.is_Exact:
        return dmp_ff_ground_primitive(f, u, K)
    else:
        return dmp_rr_ground_primitive(f, u, K)

def dup_extract(f, g, K):
    """
    Extract common content from a pair of polynomials in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dup_extract

    >>> f = ZZ.map([6, 12, 18])
    >>> g = ZZ.map([4, 8, 12])

    >>> dup_extract(f, g, ZZ)
    (2, [3, 6, 9], [2, 4, 6])

    """
    fc = dup_content(f, K)
    gc = dup_content(g, K)

    gcd = K.gcd(fc, gc)

    if not K.is_one(gcd):
        f = dup_exquo_ground(f, gcd, K)
        g = dup_exquo_ground(g, gcd, K)

    return gcd, f, g

@cythonized("u")
def dmp_ground_extract(f, g, u, K):
    """
    Extract common content from a pair of polynomials in ``K[X]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dmp_ground_extract

    >>> f = ZZ.map([[6, 12], [18]])
    >>> g = ZZ.map([[4, 8], [12]])

    >>> dmp_ground_extract(f, g, 1, ZZ)
    (2, [[3, 6], [9]], [[2, 4], [6]])

    """
    fc = dmp_ground_content(f, u, K)
    gc = dmp_ground_content(g, u, K)

    gcd = K.gcd(fc, gc)

    if not K.is_one(gcd):
        f = dmp_exquo_ground(f, gcd, u, K)
        g = dmp_exquo_ground(g, gcd, u, K)

    return gcd, f, g

def dup_real_imag(f, K):
    """
    Return bivariate polynomials ``f1`` and ``f2``, such that ``f = f1 + f2*I``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dup_real_imag

    >>> dup_real_imag([ZZ(1), ZZ(1), ZZ(1), ZZ(1)], ZZ)
    ([[1], [1], [-3, 0, 1], [-1, 0, 1]], [[3, 0], [2, 0], [-1, 0, 1, 0]])

    """
    if not K.is_ZZ and not K.is_QQ:
        raise DomainError("computing real and imaginary parts is not supported over %s" % K)

    f1 = dmp_zero(1)
    f2 = dmp_zero(1)

    if not f:
        return f1, f2

    g = [[[K.one, K.zero]], [[K.one], []]]
    h = dmp_ground(f[0], 2)

    for c in f[1:]:
        h = dmp_mul(h, g, 2, K)
        h = dmp_add_term(h, dmp_ground(c, 1), 0, 2, K)

    H = dup_to_raw_dict(h)

    for k, h in H.iteritems():
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

@cythonized('i,n')
def dup_mirror(f, K):
    """
    Evaluate efficiently the composition ``f(-x)`` in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dup_mirror

    >>> dup_mirror([ZZ(1), ZZ(2), -ZZ(4), ZZ(2)], ZZ)
    [-1, 2, 4, 2]

    """
    f, n, a = list(f), dup_degree(f), -K.one

    for i in xrange(n-1, -1, -1):
        f[i], a = a*f[i], -a

    return f

@cythonized('i,n')
def dup_scale(f, a, K):
    """
    Evaluate efficiently composition ``f(a*x)`` in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dup_scale

    >>> dup_scale([ZZ(1), -ZZ(2), ZZ(1)], ZZ(2), ZZ)
    [4, -4, 1]

    """
    f, n, b = list(f), dup_degree(f), a

    for i in xrange(n-1, -1, -1):
        f[i], b = b*f[i], b*a

    return f

@cythonized('i,j,n')
def dup_shift(f, a, K):
    """
    Evaluate efficiently Taylor shift ``f(x + a)`` in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dup_shift

    >>> dup_shift([ZZ(1), -ZZ(2), ZZ(1)], ZZ(2), ZZ)
    [1, 2, 1]

    """
    f, n = list(f), dup_degree(f)

    for i in xrange(n, 0, -1):
        for j in xrange(0, i):
            f[j+1] += a*f[j]

    return f

@cythonized('i,n')
def dup_transform(f, p, q, K):
    """
    Evaluate functional transformation ``q**n * f(p/q)`` in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dup_transform

    >>> f = ZZ.map([1, -2, 1])
    >>> p = ZZ.map([1, 0, 1])
    >>> q = ZZ.map([1, -1])

    >>> dup_transform(f, p, q, ZZ)
    [1, -2, 5, -4, 4]

    """
    if not f:
        return []

    n = dup_degree(f)
    h, Q = [f[0]], [[K.one]]

    for i in xrange(0, n):
        Q.append(dup_mul(Q[-1], q, K))

    for c, q in zip(f[1:], Q[1:]):
        h = dup_mul(h, p, K)
        q = dup_mul_ground(q, c, K)
        h = dup_add(h, q, K)

    return h

def dup_compose(f, g, K):
    """
    Evaluate functional composition ``f(g)`` in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dup_compose

    >>> f = ZZ.map([1, 1, 0])
    >>> g = ZZ.map([1, -1])

    >>> dup_compose(f, g, ZZ)
    [1, -1, 0]

    """
    if len(g) <= 1:
        return dup_strip([dup_eval(f, dup_LC(g, K), K)])

    if not f:
        return []

    h = [f[0]]

    for c in f[1:]:
        h = dup_mul(h, g, K)
        h = dup_add_term(h, c, 0, K)

    return h

@cythonized("u")
def dmp_compose(f, g, u, K):
    """
    Evaluate functional composition ``f(g)`` in ``K[X]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dmp_compose

    >>> f = ZZ.map([[1, 2], [1, 0]])
    >>> g = ZZ.map([[1, 0]])

    >>> dmp_compose(f, g, 1, ZZ)
    [[1, 3, 0]]

    """
    if not u:
        return dup_compose(f, g, K)

    if dmp_zero_p(f, u):
        return f

    h = [f[0]]

    for c in f[1:]:
        h = dmp_mul(h, g, u, K)
        h = dmp_add_term(h, c, 0, u, K)

    return h

@cythonized("s,n,r,i,j")
def _dup_right_decompose(f, s, K):
    """Helper function for :func:`_dup_decompose`."""
    n = dup_degree(f)
    lc = dup_LC(f, K)

    f = dup_to_raw_dict(f)
    g = { s : K.one }

    r = n // s

    for i in xrange(1, s):
        coeff = K.zero

        for j in xrange(0, i):
            if not n+j-i in f:
                continue

            if not s-j in g:
                continue

            fc, gc = f[n+j-i], g[s-j]
            coeff += (i - r*j)*fc*gc

        g[s-i] = K.exquo(coeff, i*r*lc)

    return dup_from_raw_dict(g, K)

@cythonized("i")
def _dup_left_decompose(f, h, K):
    """Helper function for :func:`_dup_decompose`."""
    g, i = {}, 0

    while f:
        q, r = dup_div(f, h, K)

        if dup_degree(r) > 0:
            return None
        else:
            g[i] = dup_LC(r, K)
            f, i = q, i + 1

    return dup_from_raw_dict(g, K)

@cythonized("df,s")
def _dup_decompose(f, K):
    """Helper function for :func:`dup_decompose`."""
    df = dup_degree(f)

    for s in xrange(2, df):
        if df % s != 0:
            continue

        h = _dup_right_decompose(f, s, K)

        if h is not None:
            g = _dup_left_decompose(f, h, K)

            if g is not None:
                return g, h

    return None

def dup_decompose(f, K):
    """
    Computes functional decomposition of ``f`` in ``K[x]``.

    Given an univariate polynomial ``f`` with coefficients in a field of
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

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dup_decompose

    >>> f = ZZ.map([1, -2, 1, 0, 0])

    >>> dup_decompose(f, ZZ)
    [[1, 0, 0], [1, -1, 0]]

    References
    ==========

    .. [Kozen89] D. Kozen, S. Landau, Polynomial decomposition algorithms,
    Journal of Symbolic Computation 7 (1989), pp. 445-456

    """
    F = []

    while True:
        result = _dup_decompose(f, K)

        if result is not None:
            f, h = result
            F = [h] + F
        else:
            break

    return [f] + F

@cythonized("u")
def dmp_lift(f, u, K):
    """
    Convert algebraic coefficients to integers in ``K[X]``.

    Example
    =======

    >>> from sympy import I
    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.densetools import dmp_lift

    >>> K = QQ.algebraic_field(I)
    >>> f = [K(1), K([QQ(1), QQ(0)]), K([QQ(2), QQ(0)])]

    >>> dmp_lift(f, 0, K)
    [1/1, 0/1, 2/1, 0/1, 9/1, 0/1, -8/1, 0/1, 16/1]

    """
    if not K.is_Algebraic:
        raise DomainError('computation can be done only in an algebraic domain')

    F, monoms, polys = dmp_to_dict(f, u), [], []

    for monom, coeff in F.iteritems():
        if not coeff.is_ground:
            monoms.append(monom)

    perms = variations([-1, 1], len(monoms), repetition=True)

    for perm in perms:
        G = dict(F)

        for sign, monom in zip(perm, monoms):
            if sign == -1:
                G[monom] = -G[monom]

        polys.append(dmp_from_dict(G, u, K))

    return dmp_convert(dmp_expand(polys, u, K), u, K, K.dom)

def dup_sign_variations(f, K):
    """
    Compute the number of sign variations of ``f`` in ``K[x]``.

    Example
    =======

    >>> from sympy.polys.domains import ZZ
    >>> from sympy.polys.densetools import dup_sign_variations

    >>> f = ZZ.map([1, 0, -1, -1, 1])

    >>> dup_sign_variations(f, ZZ)
    2

    """
    prev, k = K.zero, 0

    for coeff in f:
        if coeff*prev < 0:
            k += 1

        if coeff:
            prev = coeff

    return k

def dup_clear_denoms(f, K0, K1=None, convert=False):
    """
    Clear denominators, i.e. transform ``K_0`` to ``K_1``.

    Example
    =======

    >>> from sympy.polys.domains import QQ, ZZ
    >>> from sympy.polys.densetools import dup_clear_denoms

    >>> f = [QQ(1,2), QQ(1,3)]
    >>> dup_clear_denoms(f, QQ, convert=False)
    (6, [3/1, 2/1])

    >>> f = [QQ(1,2), QQ(1,3)]
    >>> dup_clear_denoms(f, QQ, convert=True)
    (6, [3, 2])

    """
    if K1 is None:
        K1 = K0.get_ring()

    common = K1.one

    for c in f:
        common = K1.lcm(common, K0.denom(c))

    if not K1.is_one(common):
        f = dup_mul_ground(f, common, K0)

    if not convert:
        return common, f
    else:
        return common, dup_convert(f, K0, K1)

@cythonized("v,w")
def _rec_clear_denoms(g, v, K0, K1):
    """Recursive helper for :func:`dmp_clear_denoms`."""
    common = K1.one

    if not v:
        for c in g:
            common = K1.lcm(common, K0.denom(c))
    else:
        w = v-1

        for c in g:
            common = K1.lcm(common, _rec_clear_denoms(c, w, K0, K1))

    return common

@cythonized("u")
def dmp_clear_denoms(f, u, K0, K1=None, convert=False):
    """
    Clear denominators, i.e. transform ``K_0`` to ``K_1``.

    Example
    =======

    >>> from sympy.polys.domains import QQ, ZZ
    >>> from sympy.polys.densetools import dmp_clear_denoms

    >>> f = [[QQ(1,2)], [QQ(1,3), QQ(1)]]
    >>> dmp_clear_denoms(f, 1, QQ, convert=False)
    (6, [[3/1], [2/1, 6/1]])

    >>> f = [[QQ(1,2)], [QQ(1,3), QQ(1)]]
    >>> dmp_clear_denoms(f, 1, QQ, convert=True)
    (6, [[3], [2, 6]])

    """
    if not u:
        return dup_clear_denoms(f, K0, K1)

    if K1 is None:
        K1 = K0.get_ring()

    common = _rec_clear_denoms(f, u, K0, K1)

    if not K1.is_one(common):
        f = dmp_mul_ground(f, common, u, K0)

    if not convert:
        return common, f
    else:
        return common, dmp_convert(f, u, K0, K1)

@cythonized('i,n')
def dup_revert(f, n, K):
    """
    Compute ``f**(-1)`` mod ``x**n`` using Newton iteration.

    This function computes first ``2**n`` terms of a polynomial that
    is a result of inversion of a polynomial modulo ``x**n``. This is
    useful to efficiently compute series expansion of ``1/f``.

    Example
    =======

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.densetools import dup_revert

    >>> f = [-QQ(1,720), QQ(0), QQ(1,24), QQ(0), -QQ(1,2), QQ(0), QQ(1)]

    >>> dup_revert(f, 8, QQ)
    [61/720, 0/1, 5/24, 0/1, 1/2, 0/1, 1/1]

    """
    g = [K.revert(dup_TC(f, K))]
    h = [K.one, K.zero, K.zero]

    N = int(ceil(log(n, 2)))

    for i in xrange(1, N + 1):
        a = dup_mul_ground(g, K(2), K)
        b = dup_mul(f, dup_sqr(g, K), K)
        g = dup_rem(dup_sub(a, b, K), h, K)
        h = dup_lshift(h, dup_degree(h), K)

    return g

def dmp_revert(f, g, u, K):
    """
    Compute ``f**(-1)`` mod ``x**n`` using Newton iteration.

    Example
    =======

    >>> from sympy.polys.domains import QQ
    >>> from sympy.polys.densetools import dmp_revert

    """
    if not u:
        return dup_revert(f, g, K)
    else:
        raise MultivariatePolynomialError(f, g)

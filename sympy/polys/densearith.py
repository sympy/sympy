"""Arithmetics for dense recursive polynomials in `K[x]` or `K[X]`. """

from sympy.polys.densebasic import (
    dup_LC, dmp_LC,
    dup_degree, dmp_degree,
    dup_normal, dmp_normal,
    dup_strip, dmp_strip,
    dmp_zero_p, dmp_zero,
    dmp_one_p, dmp_one,
    dmp_zeros
)

from sympy.polys.polyerrors import (
    ExactQuotientFailed
)

from sympy.utilities import cythonized

@cythonized("i,n,m")
def dup_add_term(f, c, i, K):
    """
    Add `c*x**i` to `f` in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dup_add_term
    >>> dup_add_term([1, 0, -1], 2, 4, ZZ)
    [2, 0, 1, 0, -1]
    """
    if not c:
        return f

    n = len(f)
    m = n-i-1

    if i == n-1:
        return dup_strip([f[0]+c] + f[1:])
    else:
        if i >= n:
            return [c] + [K.zero]*(i-n) + f
        else:
            return f[:m] + [f[m]+c] + f[m+1:]

@cythonized("i,u,v,n,m")
def dmp_add_term(f, c, i, u, K):
    """
    Add `c(x_2..x_u)*x_0**i` to `f` in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dmp_add_term
    >>> dmp_add_term([[1, 0], [1]], [2], 2, 1, ZZ)
    [[2], [1, 0], [1]]
    """
    if not u:
        return dup_add_term(f, c, i, K)

    v = u-1

    if dmp_zero_p(c, v):
        return f

    n = len(f)
    m = n-i-1

    if i == n-1:
        return dmp_strip([dmp_add(f[0], c, v, K)] + f[1:], u)
    else:
        if i >= n:
            return [c] + dmp_zeros(i-n, v, K) + f
        else:
            return f[:m] + [dmp_add(f[m], c, v, K)] + f[m+1:]

@cythonized("i,n,m")
def dup_sub_term(f, c, i, K):
    """
    Subtract `c*x**i` from `f` in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dup_sub_term
    >>> dup_sub_term([2, 0, 1, 0, -1], 2, 4, ZZ)
    [1, 0, -1]
    """
    if not c:
        return f

    n = len(f)
    m = n-i-1

    if i == n-1:
        return dup_strip([f[0]-c] + f[1:])
    else:
        if i >= n:
            return [-c] + [K.zero]*(i-n) + f
        else:
            return f[:m] + [f[m]-c] + f[m+1:]

@cythonized("i,u,v,n,m")
def dmp_sub_term(f, c, i, u, K):
    """
    Subtract `c(x_2..x_u)*x_0**i` from `f` in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dmp_sub_term
    >>> dmp_sub_term([[2], [1, 0], [1]], [2], 2, 1, ZZ)
    [[1, 0], [1]]
    """
    if not u:
        return dup_add_term(f, -c, i, K)

    v = u-1

    if dmp_zero_p(c, v):
        return f

    n = len(f)
    m = n-i-1

    if i == n-1:
        return dmp_strip([dmp_sub(f[0], c, v, K)] + f[1:], u)
    else:
        if i >= n:
            return [dmp_neg(c, v, K)] + dmp_zeros(i-n, v, K) + f
        else:
            return f[:m] + [dmp_sub(f[m], c, v, K)] + f[m+1:]

@cythonized("i")
def dup_mul_term(f, c, i, K):
    """
    Multiply `f` by `c*x**i` in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dup_mul_term
    >>> dup_mul_term([1, 0, -1], 3, 2, ZZ)
    [3, 0, -3, 0, 0]
    """
    if not c or not f:
        return []
    else:
        return [ cf * c for cf in f ] + [K.zero]*i

@cythonized("i,u,v")
def dmp_mul_term(f, c, i, u, K):
    """
    Multiply `f` by `c(x_2..x_u)*x_0**i` in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dmp_mul_term
    >>> dmp_mul_term([[1, 0], [1], []], [3, 0], 2, 1, ZZ)
    [[3, 0, 0], [3, 0], [], [], []]
    """
    if not u:
        return dup_mul_term(f, c, i, K)

    v = u-1

    if dmp_zero_p(f, u):
        return f
    if dmp_zero_p(c, v):
        return dmp_zero(u)
    else:
        return [ dmp_mul(cf, c, v, K) for cf in f ] + dmp_zeros(i, v, K)

def dup_mul_ground(f, c, K):
    """
    Multiply `f` by a constant value in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dup_mul_ground
    >>> dup_mul_ground([1, 2, -1], 3, ZZ)
    [3, 6, -3]
    """
    if not c or not f:
        return []
    else:
        return [ cf * c for cf in f ]

@cythonized("u,v")
def dmp_mul_ground(f, c, u, K):
    """
    Multiply `f` by a constant value in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dmp_mul_ground
    >>> dmp_mul_ground([[2], [2, 0]], 3, 1, ZZ)
    [[6], [6, 0]]
    """
    if not u:
        return dup_mul_ground(f, c, K)

    v = u-1

    return [ dmp_mul_ground(cf, c, v, K) for cf in f ]

def dup_quo_ground(f, c, K):
    r"""
    Quotient by a constant in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import QQ
    >>> from sympy.polys.densearith import dup_quo_ground
    >>> dup_quo_ground([QQ(1), QQ(0), QQ(2)], QQ(2), QQ) == \
    ... [QQ(1, 2), QQ(0), QQ(1)]
    True
    """
    if not c:
        raise ZeroDivisionError('polynomial division')
    if not f:
        return f

    return [ K.quo(cf, c) for cf in f ]

@cythonized("u,v")
def dmp_quo_ground(f, c, u, K):
    r"""
    Quotient by a constant in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import QQ
    >>> from sympy.polys.densearith import dmp_quo_ground
    >>> dmp_quo_ground([[QQ(1), QQ(0)], [QQ(2)], []], QQ(2), 1, QQ) == \
    ... [[QQ(1, 2), QQ(0)], [QQ(1)], []]
    True
    """
    if not u:
        return dup_quo_ground(f, c, K)

    v = u-1

    return [ dmp_quo_ground(cf, c, v, K) for cf in f ]

def dup_exquo_ground(f, c, K):
    r"""
    Exact quotient by a constant in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ, QQ
    >>> from sympy.polys.densearith import dup_exquo_ground
    >>> dup_exquo_ground([3, 0, 2], 2, ZZ)
    [1, 0, 1]
    >>> dup_exquo_ground([QQ(3), QQ(0), QQ(2)], QQ(2), QQ) == \
    ... [QQ(3, 2), QQ(0), QQ(1)]
    True
    """
    if not c:
        raise ZeroDivisionError('polynomial division')
    if not f:
        return f

    if K.has_Field or not K.is_Exact:
        return [ K.quo(cf, c) for cf in f ]
    else:
        return [ cf // c for cf in f ]

@cythonized("u,v")
def dmp_exquo_ground(f, c, u, K):
    r"""
    Exact quotient by a constant in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ, QQ
    >>> from sympy.polys.densearith import dmp_exquo_ground
    >>> dmp_exquo_ground([[2, 0], [3], []], 2, 1, ZZ)
    [[1, 0], [1], []]
    >>> dmp_exquo_ground([[QQ(2), QQ(0)], [QQ(3)], []], QQ(2), 1, QQ) == \
    ... [[QQ(1), QQ(0)], [QQ(3, 2)], []]
    True
    """
    if not u:
        return dup_exquo_ground(f, c, K)

    v = u-1

    return [ dmp_exquo_ground(cf, c, v, K) for cf in f ]

@cythonized("n")
def dup_lshift(f, n, K):
    """
    Efficiently multiply `f` by `x**n` in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dup_lshift
    >>> dup_lshift([1, 0, 1], 2, ZZ)
    [1, 0, 1, 0, 0]
    """
    if not f:
        return f
    else:
        return f + [K.zero]*n

@cythonized("n")
def dup_rshift(f, n, K):
    """
    Efficiently divide `f` by `x**n` in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dup_rshift
    >>> dup_rshift([1, 0, 1, 0, 0], 2, ZZ)
    [1, 0, 1]
    >>> dup_rshift([1, 0, 1, 0, 2], 2, ZZ)
    [1, 0, 1]
    """
    return f[:-n]

def dup_abs(f, K):
    """
    Make all coefficients positive in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dup_abs
    >>> dup_abs([1, 0, -1], ZZ)
    [1, 0, 1]
    """
    return [ K.abs(coeff) for coeff in f ]

@cythonized("u,v")
def dmp_abs(f, u, K):
    """
    Make all coefficients positive in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dmp_abs
    >>> dmp_abs([[1, 0], [-1], []], 1, ZZ)
    [[1, 0], [1], []]
    """
    if not u:
        return dup_abs(f, K)

    v = u-1

    return [ dmp_abs(cf, v, K) for cf in f ]

def dup_neg(f, K):
    """
    Negate a polynomial in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dup_neg
    >>> dup_neg([1, 0, -1], ZZ)
    [-1, 0, 1]
    """
    return [ -coeff for coeff in f ]

@cythonized("u,v")
def dmp_neg(f, u, K):
    """
    Negate a polynomial in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dmp_neg
    >>> dmp_neg([[1, 0], [-1], []], 1, ZZ)
    [[-1, 0], [1], []]
    """
    if not u:
        return dup_neg(f, K)

    v = u-1

    return [ dmp_neg(cf, u-1, K) for cf in f ]

@cythonized("df,dg,k")
def dup_add(f, g, K):
    """
    Add dense polynomials in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dup_add
    >>> dup_add([1, 0, -1], [1, -2], ZZ)
    [1, 1, -3]
    """
    if not f:
        return g
    if not g:
        return f

    df = dup_degree(f)
    dg = dup_degree(g)

    if df == dg:
        return dup_strip([ a + b for a, b in zip(f, g) ])
    else:
        k = abs(df - dg)

        if df > dg:
            h, f = f[:k], f[k:]
        else:
            h, g = g[:k], g[k:]

        return h + [ a + b for a, b in zip(f, g) ]

@cythonized("u,v,df,dg,k")
def dmp_add(f, g, u, K):
    """
    Add dense polynomials in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dmp_add
    >>> dmp_add([[1], [], [1, 0]], [[1, 0], [1], []], 1, ZZ)
    [[1, 1], [1], [1, 0]]
    """
    if not u:
        return dup_add(f, g, K)

    df = dmp_degree(f, u)

    if df < 0:
        return g

    dg = dmp_degree(g, u)

    if dg < 0:
        return f

    v = u-1

    if df == dg:
        return dmp_strip([ dmp_add(a, b, v, K) for a, b in zip(f, g) ], u)
    else:
        k = abs(df - dg)

        if df > dg:
            h, f = f[:k], f[k:]
        else:
            h, g = g[:k], g[k:]

        return h + [ dmp_add(a, b, v, K) for a, b in zip(f, g) ]

@cythonized("df,dg,k")
def dup_sub(f, g, K):
    """
    Subtract dense polynomials in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dup_sub
    >>> dup_sub([1, 0, -1], [1, -2], ZZ)
    [1, -1, 1]
    """
    if not f:
        return dup_neg(g, K)
    if not g:
        return f

    df = dup_degree(f)
    dg = dup_degree(g)

    if df == dg:
        return dup_strip([ a - b for a, b in zip(f, g) ])
    else:
        k = abs(df - dg)

        if df > dg:
            h, f = f[:k], f[k:]
        else:
            h, g = dup_neg(g[:k], K), g[k:]

        return h + [ a - b for a, b in zip(f, g) ]

@cythonized("u,v,df,dg,k")
def dmp_sub(f, g, u, K):
    """
    Subtract dense polynomials in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dmp_sub
    >>> dmp_sub([[1], [], [1, 0]], [[1, 0], [1], []], 1, ZZ)
    [[-1, 1], [-1], [1, 0]]
    """
    if not u:
        return dup_sub(f, g, K)

    df = dmp_degree(f, u)

    if df < 0:
        return dmp_neg(g, u, K)

    dg = dmp_degree(g, u)

    if dg < 0:
        return f

    v = u-1

    if df == dg:
        return dmp_strip([ dmp_sub(a, b, v, K) for a, b in zip(f, g) ], u)
    else:
        k = abs(df - dg)

        if df > dg:
            h, f = f[:k], f[k:]
        else:
            h, g = dmp_neg(g[:k], u, K), g[k:]

        return h + [ dmp_sub(a, b, v, K) for a, b in zip(f, g) ]

def dup_add_mul(f, g, h, K):
    """
    Returns `f + g*h` where `f, g, h` are in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dup_add_mul
    >>> dup_add_mul([1, 0, -1], [1, -2], [1, 2], ZZ)
    [2, 0, -5]
    """
    return dup_add(f, dup_mul(g, h, K), K)

@cythonized("u")
def dmp_add_mul(f, g, h, u, K):
    """
    Returns `f + g*h` where `f, g, h` are in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dmp_add_mul
    >>> dmp_add_mul([[1], [], [1, 0]], [[1], []], [[1], [2]], 1, ZZ)
    [[2], [2], [1, 0]]
    """
    return dmp_add(f, dmp_mul(g, h, u, K), u, K)

def dup_sub_mul(f, g, h, K):
    """
    Returns `f - g*h` where `f, g, h` are in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dup_sub_mul
    >>> dup_sub_mul([1, 0, -1], [1, -2], [1, 2], ZZ)
    [3]
    """
    return dup_sub(f, dup_mul(g, h, K), K)

@cythonized("u")
def dmp_sub_mul(f, g, h, u, K):
    """
    Returns `f - g*h` where `f, g, h` are in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dmp_sub_mul
    >>> dmp_sub_mul([[1], [], [1, 0]], [[1], []], [[1], [2]], 1, ZZ)
    [[-2], [1, 0]]
    """
    return dmp_sub(f, dmp_mul(g, h, u, K), u, K)

@cythonized("df,dg,i,j")
def dup_mul(f, g, K):
    """
    Multiply dense polynomials in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dup_mul
    >>> dup_mul([1, -2], [1, 2], ZZ)
    [1, 0, -4]
    """
    if f == g:
        return dup_sqr(f, K)

    if not (f and g):
        return []

    df = dup_degree(f)
    dg = dup_degree(g)

    h = []

    for i in xrange(0, df+dg+1):
        coeff = K.zero

        for j in xrange(max(0, i-dg), min(df, i)+1):
            coeff += f[j]*g[i-j]

        h.append(coeff)

    return h

@cythonized("u,v,df,dg,i,j")
def dmp_mul(f, g, u, K):
    """
    Multiply dense polynomials in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dmp_mul
    >>> dmp_mul([[1, 0], [1]], [[1], []], 1, ZZ)
    [[1, 0], [1], []]
    """
    if not u:
        return dup_mul(f, g, K)

    if f == g:
        return dmp_sqr(f, u, K)

    df = dmp_degree(f, u)

    if df < 0:
        return f

    dg = dmp_degree(g, u)

    if dg < 0:
        return g

    h, v = [], u-1

    for i in xrange(0, df+dg+1):
        coeff = dmp_zero(v)

        for j in xrange(max(0, i-dg), min(df, i)+1):
            coeff = dmp_add(coeff, dmp_mul(f[j], g[i-j], v, K), v, K)

        h.append(coeff)

    return h

@cythonized("df,jmin,jmax,n,i,j")
def dup_sqr(f, K):
    """
    Square dense polynomials in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dup_sqr
    >>> dup_sqr([1, 0, 1], ZZ)
    [1, 0, 2, 0, 1]
    """
    df, h = dup_degree(f), []

    for i in xrange(0, 2*df+1):
        c = K.zero

        jmin = max(0, i-df)
        jmax = min(i, df)

        n = jmax - jmin + 1

        jmax = jmin + n // 2 - 1

        for j in xrange(jmin, jmax+1):
            c += f[j]*f[i-j]

        c += c

        if n & 1:
            elem = f[jmax+1]
            c += elem**2

        h.append(c)

    return h

@cythonized("u,v,df,jmin,jmax,n,i,j")
def dmp_sqr(f, u, K):
    """
    Square dense polynomials in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dmp_sqr
    >>> dmp_sqr([[1], [1, 0], [1, 0, 0]], 1, ZZ)
    [[1], [2, 0], [3, 0, 0], [2, 0, 0, 0], [1, 0, 0, 0, 0]]
    """
    if not u:
        return dup_sqr(f, K)

    df = dmp_degree(f, u)

    if df < 0:
        return f

    h, v = [], u-1

    for i in xrange(0, 2*df+1):
        c = dmp_zero(v)

        jmin = max(0, i-df)
        jmax = min(i, df)

        n = jmax - jmin + 1

        jmax = jmin + n // 2 - 1

        for j in xrange(jmin, jmax+1):
            c = dmp_add(c, dmp_mul(f[j], f[i-j], v, K), v, K)

        c = dmp_mul_ground(c, 2, v, K)

        if n & 1:
            elem = dmp_sqr(f[jmax+1], v, K)
            c = dmp_add(c, elem, v, K)

        h.append(c)

    return h

@cythonized("n,m")
def dup_pow(f, n, K):
    """
    Raise f to the n-th power in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dup_pow
    >>> dup_pow([1, -2], 3, ZZ)
    [1, -6, 12, -8]
    """
    if not n:
        return [K.one]
    if n < 0:
        raise ValueError("can't raise polynomial to a negative power")
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

@cythonized("u,n,m")
def dmp_pow(f, n, u, K):
    """
    Raise f to the n-th power in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dmp_pow
    >>> dmp_pow([[1, 0], [1]], 3, 1, ZZ)
    [[1, 0, 0, 0], [3, 0, 0], [3, 0], [1]]
    """
    if not u:
        return dup_pow(f, n, K)

    if not n:
        return dmp_one(u, K)
    if n < 0:
        raise ValueError("can't raise polynomial to a negative power")
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

@cythonized("df,dg,dr,N,j")
def dup_pdiv(f, g, K):
    """
    Polynomial pseudo-division in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dup_pdiv
    >>> dup_pdiv([1, 0, 1], [2, -4], ZZ)
    ([2, 4], [20])
    """
    df = dup_degree(f)
    dg = dup_degree(g)

    q, r = [], f

    if not g:
        raise ZeroDivisionError("polynomial division")
    elif df < dg:
        return q, r

    N = df - dg + 1
    lc_g = dup_LC(g, K)

    while True:
        dr = dup_degree(r)

        if dr < dg:
            break

        lc_r = dup_LC(r, K)
        j, N = dr-dg, N-1

        Q = dup_mul_ground(q, lc_g, K)
        q = dup_add_term(Q, lc_r, j, K)

        R = dup_mul_ground(r, lc_g, K)
        G = dup_mul_term(g, lc_r, j, K)
        r = dup_sub(R, G, K)

    c = lc_g**N

    q = dup_mul_ground(q, c, K)
    r = dup_mul_ground(r, c, K)

    return q, r

@cythonized("df,dg,dr,N,j")
def dup_prem(f, g, K):
    """
    Polynomial pseudo-remainder in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dup_prem
    >>> dup_prem([1, 0, 1], [2, -4], ZZ)
    [20]
    """
    df = dup_degree(f)
    dg = dup_degree(g)

    r = f

    if not g:
        raise ZeroDivisionError("polynomial division")
    elif df < dg:
        return r

    N = df - dg + 1
    lc_g = dup_LC(g, K)

    while True:
        dr = dup_degree(r)

        if dr < dg:
            break

        lc_r = dup_LC(r, K)
        j, N = dr-dg, N-1

        R = dup_mul_ground(r, lc_g, K)
        G = dup_mul_term(g, lc_r, j, K)
        r = dup_sub(R, G, K)

    return dup_mul_ground(r, lc_g**N, K)

def dup_pquo(f, g, K):
    """
    Polynomial pseudo-quotient in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dup_pquo
    >>> dup_pquo([1, 0, 1], [2, -4], ZZ)
    Traceback (most recent call last):
    ...
    ExactQuotientFailed: [2, -4] does not divide [1, 0, 1]
    >>> dup_pquo([1, 0, -1], [2, -2], ZZ)
    [2, 2]
    """
    q, r = dup_pdiv(f, g, K)

    if not r:
        return q
    else:
        raise ExactQuotientFailed('%s does not divide %s' % (g, f))

def dup_pexquo(f, g, K):
    """
    Polynomial exact pseudo-quotient in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dup_pexquo
    >>> dup_pexquo([1, 0, 1], [2, -4], ZZ)
    [2, 4]
    >>> dup_pexquo([1, 0, -1], [2, -2], ZZ)
    [2, 2]
    """
    return dup_pdiv(f, g, K)[0]

@cythonized("u,df,dg,dr,N,j")
def dmp_pdiv(f, g, u, K):
    """
    Polynomial pseudo-division in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dmp_pdiv
    >>> dmp_pdiv([[1], [1, 0], []], [[2], [2]], 1, ZZ)
    ([[2], [2, -2]], [[-4, 4]])
    """
    if not u:
        return dup_pdiv(f, g, K)

    df = dmp_degree(f, u)
    dg = dmp_degree(g, u)

    if dg < 0:
        raise ZeroDivisionError("polynomial division")

    q, r = dmp_zero(u), f

    if df < dg:
        return q, r

    N = df - dg + 1
    lc_g = dmp_LC(g, K)

    while True:
        dr = dmp_degree(r, u)

        if dr < dg:
            break

        lc_r = dmp_LC(r, K)
        j, N = dr-dg, N-1

        Q = dmp_mul_term(q, lc_g, 0, u, K)
        q = dmp_add_term(Q, lc_r, j, u, K)

        R = dmp_mul_term(r, lc_g, 0, u, K)
        G = dmp_mul_term(g, lc_r, j, u, K)
        r = dmp_sub(R, G, u, K)

    c = dmp_pow(lc_g, N, u-1, K)

    q = dmp_mul_term(q, c, 0, u, K)
    r = dmp_mul_term(r, c, 0, u, K)

    return q, r

@cythonized("u,df,dg,dr,N,j")
def dmp_prem(f, g, u, K):
    """
    Polynomial pseudo-remainder in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dmp_prem
    >>> dmp_prem([[1], [1, 0], []], [[2], [2]], 1, ZZ)
    [[-4, 4]]
    """
    if not u:
        return dup_prem(f, g, K)

    df = dmp_degree(f, u)
    dg = dmp_degree(g, u)

    if dg < 0:
        raise ZeroDivisionError("polynomial division")

    r = f

    if df < dg:
        return r

    N = df - dg + 1
    lc_g = dmp_LC(g, K)

    while True:
        dr = dmp_degree(r, u)

        if dr < dg:
            break

        lc_r = dmp_LC(r, K)
        j, N = dr-dg, N-1

        R = dmp_mul_term(r, lc_g, 0, u, K)
        G = dmp_mul_term(g, lc_r, j, u, K)
        r = dmp_sub(R, G, u, K)

    c = dmp_pow(lc_g, N, u-1, K)

    return dmp_mul_term(r, c, 0, u, K)

def dmp_pquo(f, g, u, K):
    """
    Polynomial pseudo-quotient in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dmp_pquo
    >>> dmp_pquo([[1], [1, 0], []], [[2], [2]], 1, ZZ)
    Traceback (most recent call last):
    ...
    ExactQuotientFailed: [[2], [2]] does not divide [[1], [1, 0], []]
    >>> dmp_pquo([[1], [1, 0], []], [[2], [2, 0]], 1, ZZ)
    [[2], []]
    """
    q, r = dmp_pdiv(f, g, u, K)

    if dmp_zero_p(r, u):
        return q
    else:
        raise ExactQuotientFailed('%s does not divide %s' % (g, f))

def dmp_pexquo(f, g, u, K):
    """
    Polynomial exact pseudo-quotient in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dmp_pexquo
    >>> dmp_pexquo([[1], [1, 0], []], [[2], [2]], 1, ZZ)
    [[2], [2, -2]]
    >>> dmp_pexquo([[1], [1, 0], []], [[2], [2, 0]], 1, ZZ)
    [[2], []]
    """
    return dmp_pdiv(f, g, u, K)[0]

@cythonized("df,dg,dr,j")
def dup_rr_div(f, g, K):
    """
    Univariate division with remainder over a ring.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dup_rr_div
    >>> dup_rr_div([1, 0, 1], [2, -4], ZZ)
    ([], [1, 0, 1])
    """
    df = dup_degree(f)
    dg = dup_degree(g)

    q, r = [], f

    if not g:
        raise ZeroDivisionError("polynomial division")
    elif df < dg:
        return q, r

    lc_g = dup_LC(g, K)

    while True:
        dr = dup_degree(r)

        if dr < dg:
            break

        lc_r = dup_LC(r, K)

        if lc_r % lc_g:
            break

        c = K.exquo(lc_r, lc_g)
        j = dr - dg

        q = dup_add_term(q, c, j, K)
        h = dup_mul_term(g, c, j, K)

        r = dup_sub(r, h, K)

    return q, r

@cythonized("u,df,dg,dr,j")
def dmp_rr_div(f, g, u, K):
    """
    Multivariate division with remainder over a ring.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dmp_rr_div
    >>> dmp_rr_div([[1], [1, 0], []], [[2], [2]], 1, ZZ)
    ([[]], [[1], [1, 0], []])
    """
    if not u:
        return dup_rr_div(f, g, K)

    df = dmp_degree(f, u)
    dg = dmp_degree(g, u)

    if dg < 0:
        raise ZeroDivisionError("polynomial division")

    q, r = dmp_zero(u), f

    if df < dg:
        return q, r

    lc_g, v = dmp_LC(g, K), u-1

    while True:
        dr = dmp_degree(r, u)

        if dr < dg:
            break

        lc_r = dmp_LC(r, K)

        c, R = dmp_rr_div(lc_r, lc_g, v, K)

        if not dmp_zero_p(R, v):
            break

        j = dr - dg

        q = dmp_add_term(q, c, j, u, K)
        h = dmp_mul_term(g, c, j, u, K)

        r = dmp_sub(r, h, u, K)

    return q, r

@cythonized("df,dg,dr,j")
def dup_ff_div(f, g, K):
    r"""Polynomial division with remainder over a field.

    Example
    =======
    >>> from sympy.polys.algebratools import QQ
    >>> from sympy.polys.densearith import dup_ff_div
    >>> dup_ff_div([QQ(1), QQ(0), QQ(1)], [QQ(2), QQ(-4)], QQ) == \
    ... ([QQ(1,2), QQ(1)], [QQ(5)])
    True
    """
    df = dup_degree(f)
    dg = dup_degree(g)

    q, r = [], f

    if not g:
        raise ZeroDivisionError("polynomial division")
    elif df < dg:
        return q, r

    lc_g = dup_LC(g, K)

    while True:
        dr = dup_degree(r)

        if dr < dg:
            break

        lc_r = dup_LC(r, K)

        c = K.exquo(lc_r, lc_g)
        j = dr - dg

        q = dup_add_term(q, c, j, K)
        h = dup_mul_term(g, c, j, K)

        r = dup_sub(r, h, K)

        if not K.is_Exact:
            r = dup_normal(r, K)

    return q, r

@cythonized("u,df,dg,dr,j")
def dmp_ff_div(f, g, u, K):
    r"""
    Polynomial division with remainder over a field.

    Example
    =======
    >>> from sympy.polys.algebratools import QQ
    >>> from sympy.polys.densearith import dmp_ff_div
    >>> dmp_ff_div([[QQ(1)], [QQ(1), QQ(0)], []], [[QQ(2)], [QQ(2)]],
    ... 1, QQ) == \
    ... ([[QQ(1, 2)], [QQ(1, 2), QQ(-1, 2)]], [[QQ(-1, 1), QQ(1)]])
    True
    """
    if not u:
        return dup_ff_div(f, g, K)

    df = dmp_degree(f, u)
    dg = dmp_degree(g, u)

    if dg < 0:
        raise ZeroDivisionError("polynomial division")

    q, r = dmp_zero(u), f

    if df < dg:
        return q, r

    lc_g, v = dmp_LC(g, K), u-1

    while True:
        dr = dmp_degree(r, u)

        if dr < dg:
            break

        lc_r = dmp_LC(r, K)

        c, R = dmp_ff_div(lc_r, lc_g, v, K)

        if not dmp_zero_p(R, v):
            break

        j = dr - dg

        q = dmp_add_term(q, c, j, u, K)
        h = dmp_mul_term(g, c, j, u, K)

        r = dmp_sub(r, h, u, K)

    return q, r

def dup_div(f, g, K):
    r"""Polynomial division with remainder in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ, QQ
    >>> from sympy.polys.densearith import dup_div
    >>> dup_div([1, 0, 1], [2, -4], ZZ)
    ([], [1, 0, 1])
    >>> dup_div([QQ(1), QQ(0), QQ(1)], [QQ(2), QQ(-4)], QQ) == \
    ... ([QQ(1,2), QQ(1)], [QQ(5)])
    True
    """
    if K.has_Field or not K.is_Exact:
        return dup_ff_div(f, g, K)
    else:
        return dup_rr_div(f, g, K)

def dup_rem(f, g, K):
    r"""
    Returns polynomial remainder in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ, QQ
    >>> from sympy.polys.densearith import dup_rem
    >>> dup_rem([1, 0, 1], [2, -4], ZZ)
    [1, 0, 1]
    >>> dup_rem([QQ(1), QQ(0), QQ(1)], [QQ(2), QQ(-4)], QQ) == \
    ... [QQ(5)]
    True
    """
    return dup_div(f, g, K)[1]

def dup_quo(f, g, K):
    """
    Returns polynomial quotient in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ, QQ
    >>> from sympy.polys.densearith import dup_quo
    >>> dup_quo([1, 0, 1], [2, -4], ZZ)
    Traceback (most recent call last):
    ...
    ExactQuotientFailed: [2, -4] does not divide [1, 0, 1]
    >>> dup_quo([1, 0, -1], [1, -1], ZZ)
    [1, 1]
    """
    q, r = dup_div(f, g, K)

    if not r:
        return q
    else:
        raise ExactQuotientFailed('%s does not divide %s' % (g, f))

def dup_exquo(f, g, K):
    r"""
    Returns exact polynomial quotient in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ, QQ
    >>> from sympy.polys.densearith import dup_exquo
    >>> dup_exquo([1, 0, 1], [2, -4], ZZ)
    []
    >>> dup_exquo([QQ(1), QQ(0), QQ(1)], [QQ(2), QQ(-4)], QQ) == \
    ... [QQ(1,2), QQ(1)]
    True
    """
    return dup_div(f, g, K)[0]

@cythonized("u")
def dmp_div(f, g, u, K):
    r"""
    Polynomial division with remainder in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ, QQ
    >>> from sympy.polys.densearith import dmp_div
    >>> dmp_div([[1], [1, 0], []], [[2], [2]], 1, ZZ)
    ([[]], [[1], [1, 0], []])
    >>> dmp_div([[QQ(1)], [QQ(1), QQ(0)], []], [[QQ(2)], [QQ(2)]],
    ... 1, QQ) == \
    ... ([[QQ(1, 2)], [QQ(1, 2), QQ(-1, 2)]], [[QQ(-1, 1), QQ(1)]])
    True
    """
    if K.has_Field or not K.is_Exact:
        return dmp_ff_div(f, g, u, K)
    else:
        return dmp_rr_div(f, g, u, K)

@cythonized("u")
def dmp_rem(f, g, u, K):
    r"""
    Returns polynomial remainder in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ, QQ
    >>> from sympy.polys.densearith import dmp_rem
    >>> dmp_rem([[1], [1, 0], []], [[2], [2]], 1, ZZ)
    [[1], [1, 0], []]
    >>> dmp_rem([[QQ(1)], [QQ(1), QQ(0)], []], [[QQ(2)], [QQ(2)]],
    ... 1, QQ) == \
    ... [[QQ(-1, 1), QQ(1)]]
    True
    """
    return dmp_div(f, g, u, K)[1]

@cythonized("u")
def dmp_quo(f, g, u, K):
    """
    Returns polynomial quotient in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ, QQ
    >>> from sympy.polys.densearith import dmp_quo
    >>> dmp_quo([[1], [1, 0], []], [[2], [2]], 1, ZZ)
    Traceback (most recent call last):
    ...
    ExactQuotientFailed: [[2], [2]] does not divide [[1], [1, 0], []]
    >>> dmp_quo([[1], [1, 0], []], [[1], [1, 0]], 1, ZZ)
    [[1], []]
    """
    q, r = dmp_div(f, g, u, K)

    if dmp_zero_p(r, u):
        return q
    else:
        raise ExactQuotientFailed('%s does not divide %s' % (g, f))

@cythonized("u")
def dmp_exquo(f, g, u, K):
    r"""
    Returns exact polynomial quotient in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ, QQ
    >>> from sympy.polys.densearith import dmp_exquo
    >>> dmp_exquo([[1], [1, 0], []], [[2], [2]], 1, ZZ)
    [[]]
    >>> dmp_exquo([[QQ(1)], [QQ(1), QQ(0)], []], [[QQ(2)], [QQ(2)]],
    ... 1, QQ) == \
    ... [[QQ(1, 2)], [QQ(1, 2), QQ(-1, 2)]]
    True
    """
    return dmp_div(f, g, u, K)[0]

def dup_max_norm(f, K):
    """
    Returns maximum norm of a polynomial in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dup_max_norm
    >>> dup_max_norm([-1, 2, -3], ZZ)
    3
    """
    if not f:
        return K.zero
    else:
        return max(dup_abs(f, K))

@cythonized("u,v")
def dmp_max_norm(f, u, K):
    """
    Returns maximum norm of a polynomial in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dmp_max_norm
    >>> dmp_max_norm([[2, -1], [-3]], 1, ZZ)
    3
    """
    if not u:
        return dup_max_norm(f, K)

    v = u-1

    return max([ dmp_max_norm(c, v, K) for c in f ])

def dup_l1_norm(f, K):
    """
    Returns l1 norm of a polynomial in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dup_l1_norm
    >>> dup_l1_norm([2, -3, 0, 1], ZZ)
    6
    """
    if not f:
        return K.zero
    else:
        return sum(dup_abs(f, K))

@cythonized("u,v")
def dmp_l1_norm(f, u, K):
    """
    Returns l1 norm of a polynomial in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dmp_l1_norm
    >>> dmp_l1_norm([[2, -1], [-3]], 1, ZZ)
    6
    """
    if not u:
        return dup_l1_norm(f, K)

    v = u-1

    return sum([ dmp_l1_norm(c, v, K) for c in f ])

def dup_expand(polys, K):
    """
    Multiply together several polynomials in `K[x]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dup_expand
    >>> dup_expand(([1, 0, -1], [1, 0], [2]), ZZ)
    [2, 0, -2, 0]
    """
    if not polys:
        return [K.one]

    f = polys[0]

    for g in polys[1:]:
        f = dup_mul(f, g, K)

    return f

@cythonized("u")
def dmp_expand(polys, u, K):
    """
    Multiply together several polynomials in `K[X]`.

    Example
    =======
    >>> from sympy.polys.algebratools import ZZ
    >>> from sympy.polys.densearith import dmp_expand
    >>> dmp_expand(([[1], [], [1, 0, 0]], [[1], [1]]), 1, ZZ)
    [[1], [1], [1, 0, 0], [1, 0, 0]]
    """
    if not polys:
        return dmp_one(u, K)

    f = polys[0]

    for g in polys[1:]:
        f = dmp_mul(f, g, u, K)

    return f


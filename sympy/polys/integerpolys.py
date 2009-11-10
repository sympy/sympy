"""Univariate and multivariate polynomials with coefficients in the integer ring. """

from sympy.polys.galoispolys import (
    gf_from_int_poly, gf_to_int_poly, gf_degree, gf_from_dict,
    gf_lshift, gf_add_mul, gf_mul, gf_div, gf_quo, gf_rem,
    gf_gcd, gf_gcdex, gf_sqf_p, gf_factor_sqf)

from sympy.ntheory import randprime, nextprime, isprime, factorint
from sympy.ntheory.modular import crt1, crt2
from sympy.utilities import any, all, subsets

from sympy.core.numbers import Integer

from math import floor, ceil, log, sqrt
from random import randint

from sympy.core.numbers import igcd, igcdex
from sympy.mpmath.libmpf import isqrt

INT_TYPE = int
INT_ZERO = 0
INT_ONE  = 1

from copy import deepcopy

def factorial(m):
    k = m

    while m > 1:
        m -= 1
        k *= m

    return k

class ExactQuotientFailed(Exception):
    pass

class HeuristicGCDFailed(Exception):
    pass

class ExtraneousFactors(Exception):
    pass

def poly_LC(f):
    """Returns leading coefficient of f. """
    if not f:
        return INT_ZERO
    else:
        return f[0]

def poly_TC(f):
    """Returns trailing coefficient of f. """
    if not f:
        return INT_ZERO
    else:
        return f[-1]

def poly_nth(f, n):
    """Returns n-th coefficient of f. """
    if n < 0 or n > len(f)-1:
        raise IndexError
    else:
        return f[zzx_degree(f)-n]

def poly_level(f):
    """Return the number of nested lists in f. """
    if poly_univariate_p(f):
        return 1
    else:
        return 1 + poly_level(poly_LC(f))

def poly_univariate_p(f):
    """Returns True if f is univariate. """
    if not f:
        return True
    else:
        return type(f[0]) is not list

def zzx_degree(f):
    """Returns leading degree of f in Z[x]. """
    return len(f) - 1

def zzX_degree(f):
    """Returns leading degree of f in x_1 in Z[X]. """
    if zzX_zero_p(f):
        return -1
    else:
        return len(f) - 1

def zzX_degree_for(f, k):
    """Returns leading degree of f in x_k in Z[X]. """
    if k < 0:
        k += poly_level(f) + 1

    if k == 1:
        return zzX_degree(f)

    def rec_degree(g, l):
        if l == k:
            return zzX_degree(g)
        else:
            return max([ rec_degree(coeff, l+1) for coeff in g ])

    return rec_degree(f, 1)

def zzX_degree_all(f):
    """Returns total degree of f in Z[X]. """
    degs = [-1]*poly_level(f)

    def rec_degree(g, l):
        degs[l-1] = max(degs[l-1], zzX_degree(g))

        if not poly_univariate_p(g):
            for coeff in g:
                rec_degree(coeff, l+1)

    rec_degree(f, 1)

    return tuple(degs)

def zzx_strip(f):
    """Remove leading zeros from f in Z[x]. """
    if not f or f[0]:
        return f

    k = 0

    for coeff in f:
        if coeff:
            break
        else:
            k += 1

    return f[k:]

def zzX_strip(f):
    """Remove leading zeros from f in Z[X]. """
    if poly_univariate_p(f):
        return zzx_strip(f)

    if zzX_zero_p(f):
        return f

    k = 0

    for coeff in f:
        if not zzX_zero_p(coeff):
            break
        else:
            k += 1

    if k == len(f):
        return zzX_zero_of(f)
    else:
        return f[k:]

def zzX_valid_p(f):
    """Returns True if f is a valid polynomial in Z[x]. """
    levels = []

    def rec_valid(g, l):
        if poly_univariate_p(g):
            levels.append(l)
            return zzx_strip(g) == g
        else:
            return zzX_strip(g) == g and \
                all([ rec_valid(h, l+1) for h in g ])

    return rec_valid(f, 1) and len(set(levels)) == 1

def zzX_zz_LC(f):
    """Returns integer leading coefficient. """
    if poly_univariate_p(f):
        return poly_LC(f)
    else:
        return zzX_zz_LC(poly_LC(f))

def zzX_zz_TC(f):
    """Returns integer trailing coefficient. """
    if poly_univariate_p(f):
        return poly_TC(f)
    else:
        return zzX_zz_TC(poly_TC(f))

def zzX_zero(l):
    """Returns multivariate zero. """
    if not l:
        return INT_ZERO
    elif l == 1:
        return []
    else:
        return [zzX_zero(l-1)]

def zzX_zero_of(f, d=0):
    """Returns multivariate zero of f. """
    return zzX_zero(poly_level(f)-d)

def zzX_const(l, c):
    """Returns multivariate constant. """
    if not c:
        return zzX_zero(l)
    else:
        if not l:
            return INT_TYPE(c)
        elif l == 1:
            return [INT_TYPE(c)]
        else:
            return [zzX_const(l-1, c)]

def zzX_const_of(f, c, d=0):
    """Returns multivariate constant of f. """
    return zzX_const(poly_level(f)-d, c)

def zzX_zeros_of(f, k, d=0):
    """Returns a list of multivariate zeros of f. """
    if poly_univariate_p(f):
        return [INT_ZERO]*k

    l = poly_level(f)-d

    if not k:
        return []
    else:
        return [ zzX_zero(l) for i in xrange(k) ]

def zzX_consts_of(f, c, k, d=0):
    """Returns a list of multivariate constants of f. """
    if poly_univariate_p(f):
        return [INT_TYPE(c)]*k

    l = poly_level(f)-d

    if not k:
        return []
    else:
        return [ zzX_const(l, c) for i in xrange(k) ]

def zzX_zero_p(f):
    """Returns True if f is zero in Z[X]. """
    if poly_univariate_p(f):
        return not f
    else:
        if len(f) == 1:
            return zzX_zero_p(f[0])
        else:
            return False

def zzx_one_p(f):
    """Returns True if f is one in Z[x]. """
    return f == [INT_ONE]

def zzX_one_p(f):
    """Returns True if f is one in Z[X]. """
    if poly_univariate_p(f):
        return zzx_one_p(f)
    else:
        if len(f) == 1:
            return zzX_one_p(f[0])
        else:
            return False

def zzX_value(l, f):
    """Returns multivariate value nested l-levels. """
    if type(f) is not list:
        return zzX_const(l, f)
    else:
        if not l:
            return f
        else:
            return [zzX_value(l-1, f)]

def zzX_lift(l, f):
    """Returns multivariate polynomial lifted l-levels. """
    if poly_univariate_p(f):
        if not f:
            return zzX_zero(l+1)
        else:
            return [ zzX_const(l, c) for c in f ]
    else:
        return [ zzX_lift(l, c) for c in f ]

def zzx_from_dict(f):
    """Create Z[x] polynomial from a dict. """
    if not f:
        return []

    n, h = max(f.iterkeys()), []

    for k in xrange(n, -1, -1):
        h.append(INT_TYPE(int(f.get(k, 0))))

    return zzx_strip(h)

def zzX_from_dict(f, l):
    """Create Z[X] polynomial from a dict. """
    if l == 1:
        return zzx_from_dict(f)
    elif not f:
        return zzX_zero(l)

    coeffs = {}

    for monom, coeff in f.iteritems():
        head, tail = monom[0], monom[1:]

        if len(tail) == 1:
            tail = tail[0]

        if head in coeffs:
            coeffs[head][tail] = INT_TYPE(int(coeff))
        else:
            coeffs[head] = { tail : INT_TYPE(int(coeff)) }

    n, h = max(coeffs.iterkeys()), []

    for k in xrange(n, -1, -1):
        coeff = coeffs.get(k)

        if coeff is not None:
            h.append(zzX_from_dict(coeff, l-1))
        else:
            h.append(zzX_zero(l-1))

    return zzX_strip(h)

def zzx_to_dict(f):
    """Convert Z[x] polynomial to a dict. """
    n, result = zzx_degree(f), {}

    for i in xrange(0, n+1):
        if f[n-i]:
            result[i] = f[n-i]

    return result

def zzX_to_dict(f):
    """Convert Z[X] polynomial to a dict. """
    if poly_univariate_p(f):
        return zzx_to_dict(f)

    n, result = zzX_degree(f), {}

    for i in xrange(0, n+1):
        h = zzX_to_dict(f[n-i])

        for exp, coeff in h.iteritems():
            if type(exp) is not tuple:
                exp = (exp,)

            result[(i,)+exp] = coeff

    return result

def zzx_from_poly(f):
    """Convert Poly instance to a recursive dense polynomial in Z[x]. """
    return zzx_from_dict(dict(zip([ m for (m,) in f.monoms ], f.coeffs)))

def zzX_from_poly(f):
    """Convert Poly instance to a recursive dense polynomial in Z[X]. """
    if f.is_univariate:
        return zzx_from_poly(f)
    else:
        return zzX_from_dict(dict(zip(f.monoms, f.coeffs)), len(f.symbols))

def zzx_to_poly(f, *symbols):
    """Convert recursive dense polynomial to a Poly in Z[x]. """
    from sympy.polys import Poly

    terms = {}

    for monom, coeff in zzx_to_dict(f).iteritems():
        terms[(monom,)] = Integer(int(coeff))

    return Poly(terms, *symbols)

def zzX_to_poly(f, *symbols):
    """Convert recursive dense polynomial to a Poly in Z[X]. """
    from sympy.polys import Poly

    terms = {}

    for monom, coeff in zzX_to_dict(f).iteritems():
        terms[monom] = Integer(int(coeff))

    return Poly(terms, *symbols)

def zzX_swap(f, i=1, j=2):
    """Transform Z[..x_i..x_j..] to Z[..x_j..x_i..]. """
    l = poly_level(f)

    if i < 1 or j < 1 or i > l or j > l:
        raise ValueError("1 <= i < j <= lev(f) expected")
    elif i == j:
        return f
    else:
        i, j = i-1, j-1

    F, H = zzX_to_dict(f), {}

    for exp, coeff in F.iteritems():
        H[exp[:i]    + (exp[j],) +
          exp[i+1:j] +
          (exp[i],)  + exp[j+1:]] = coeff

    return zzX_from_dict(H, l)

def zzx_abs(f):
    """Make all coefficients positive in Z[x]. """
    return [ abs(coeff) for coeff in f ]

def zzX_abs(f):
    """Make all coefficients positive in Z[X]. """
    if poly_univariate_p(f):
        return zzx_abs(f)
    else:
        return [ zzX_abs(coeff) for coeff in f ]

def zzx_neg(f):
    """Negate a polynomial in Z[x]. """
    return [ -coeff for coeff in f ]

def zzX_neg(f):
    """Negate a polynomial in Z[X]. """
    if poly_univariate_p(f):
        return zzx_neg(f)
    else:
        return [ zzX_neg(coeff) for coeff in f ]

def zzx_add_term(f, c, k=0):
    """Add c*x**k to f in Z[x]. """
    if not c:
        return f

    n = len(f)
    m = n-k-1

    if k == n-1:
        return zzx_strip([f[0]+c] + f[1:])
    else:
        if k >= n:
            return [c] + [INT_ZERO]*(k-n) + f
        else:
            return f[:m] + [f[m]+c] + f[m+1:]

def zzX_add_term(f, c, k=0):
    """Add c*x**k to f in Z[X]. """
    if poly_univariate_p(f):
        return zzx_add_term(f, c, k)

    if zzX_zero_p(c):
        return f

    n = len(f)
    m = n-k-1

    if k == n-1:
        return zzX_strip([zzX_add(f[0], c)] + f[1:])
    else:
        if k >= n:
            return [c] + zzX_zeros_of(f, k-n, 1) + f
        else:
            return f[:m] + [zzX_add(f[m], c)] + f[m+1:]

def zzx_sub_term(f, c, k=0):
    """Subtract c*x**k from f in Z[x]. """
    if not c:
        return f

    n = len(f)
    m = n-k-1

    if k == n-1:
        return zzx_strip([f[0]-c] + f[1:])
    else:
        if k >= n:
            return [-c] + [INT_ZERO]*(k-n) + f
        else:
            return f[:m] + [f[m]-c] + f[m+1:]

def zzX_sub_term(f, c, k=0):
    """Subtract c*x**k from f in Z[X]. """
    return zzX_add_term(f, zzX_neg(c), k)

def zzx_mul_term(f, c, k):
    """Multiply f by c*x**k in Z[x]. """
    if not c or not f:
        return []
    else:
        return [ c * coeff for coeff in f ] + [INT_ZERO]*k

def zzX_mul_term(f, c, k):
    """Multiply f by c*x**k in Z[X]. """
    if poly_univariate_p(f):
        return zzx_mul_term(f, c, k)
    elif zzX_zero_p(f):
        return f
    elif zzX_zero_p(c):
        return zzX_zero_of(f)
    else:
        return [ zzX_mul(c, coeff) for coeff in f ] + zzX_zeros_of(f, k, 1)

def zzx_mul_const(f, c):
    """Multiply f by constant value in Z[x]. """
    if not c or not f:
        return []
    else:
        return [ c * coeff for coeff in f ]

def zzX_mul_const(f, c):
    """Multiply f by constant value in Z[X]. """
    if poly_univariate_p(f):
        return zzx_mul_const(f, c)
    else:
        return [ zzX_mul_const(coeff, c) for coeff in f ]

def zzx_quo_const(f, c):
    """Exact quotient by a constant in Z[x]. """
    if not c:
        raise ZeroDivisionError('polynomial division')
    elif not f:
        return f
    else:
        h = []

        for coeff in f:
            if coeff % c:
                raise ExactQuotientFailed('%s does not divide %s' % (c, coeff))
            else:
                h.append(coeff // c)

        return h

def zzX_quo_const(f, c):
    """Exact quotient by a constant in Z[X]. """
    if poly_univariate_p(f):
        return zzx_quo_const(f, c)
    else:
        return [ zzX_quo_const(coeff, c) for coeff in f ]

def zzx_compose_term(f, k):
    """Map y -> x**k in a polynomial in Z[x]. """
    if k <= 0:
        raise ValueError("'k' must be positive, got %s" % k)
    if k == 1 or not f:
        return f

    result = [f[0]]

    for coeff in f[1:]:
        result.extend([0]*(k-1))
        result.append(coeff)

    return result

def zzX_compose_term(f, K):
    """Map y_i -> x_i**k_i in a polynomial in Z[X]. """
    def rec_compose(g, l):
        if poly_univariate_p(g):
            return zzx_compose_term(g, K[l])

        if K[l] <= 0:
            raise ValueError("All 'K[i]' must be positive, got %s" % K[l])

        g = [ rec_compose(c, l+1) for c in g ]
        result, L = [g[0]], poly_level(g) - 1

        for coeff in g[1:]:
            for i in xrange(1, K[l]):
                result.append(zzX_zero(L))

            result.append(coeff)

        return result

    if all([ k == 1 for k in K ]):
        return f
    else:
        return rec_compose(f, 0)

def zzx_reduce(f):
    """Map x**k -> y in a polynomial in Z[x]. """
    if zzx_degree(f) <= 0:
        return 1, f

    g = INT_ZERO

    for i in xrange(len(f)):
        if not f[-i-1]:
            continue

        g = igcd(g, i)

        if g == 1:
            return 1, f

    return g, f[::g]

def zzX_reduce(f):
    """Map x_i**k_i -> y_i in a polynomial in Z[X]. """
    if zzX_zero_p(f):
        return (1,)*poly_level(f), f

    F, H = zzX_to_dict(f), {}

    def ilgcd(M):
        g = 0

        for m in M:
            g = igcd(g, m)

            if g == 1:
                break

        return g or 1

    M = tuple(map(lambda *row: ilgcd(row), *F.keys()))

    if all([ b == 1 for b in M ]):
        return M, f

    for m, coeff in F.iteritems():
        N = [ a // b for a, b in zip(m, M) ]
        H[tuple(N)] = coeff

    return M, zzX_from_dict(H, len(M))

def zzx_multi_reduce(*polys):
    """Map x**k -> y in a set of polynomials in Z[x]. """
    G = INT_ZERO

    for p in polys:
        if zzx_degree(p) <= 0:
            return 1, polys

        g = INT_ZERO

        for i in xrange(len(p)):
            if not p[-i-1]:
                continue

            g = igcd(g, i)

            if g == 1:
                return 1, polys

        G = igcd(G, g)

        if G == 1:
            return 1, polys

    return G, tuple([ p[::G] for p in polys ])

def zzX_multi_reduce(*polys):
    """Map x_i**k_i -> y_i in a set of polynomials in Z[X]. """
    def ilgcd(M):
        g = 0

        for m in M:
            g = igcd(g, m)

            if g == 1:
                break

        return g or 1

    l = poly_level(polys[0])

    if l == 1:
        M, H = zzx_multi_reduce(*polys)
        return (M,), H

    F, M, H = [], [], []

    for p in polys:
        f = zzX_to_dict(p)

        if zzX_zero_p(p):
            m = (0,)*l
        else:
            m = map(lambda *row: ilgcd(row), *f.keys())

        F.append(f)
        M.append(m)

    M = tuple(map(lambda *row: ilgcd(row), *M))

    if all([ b == 1 for b in M ]):
        return M, polys

    for f in F:
        h = {}

        for m, coeff in f.iteritems():
            N = [ a // b for a, b in zip(m, M) ]
            h[tuple(N)] = coeff

        H.append(zzX_from_dict(h, len(m)))

    return M, tuple(H)

def zzx_add(f, g):
    """Add polynomials in Z[x]. """
    if not f:
        return g
    if not g:
        return f

    df = zzx_degree(f)
    dg = zzx_degree(g)

    if df == dg:
        return zzx_strip([ a + b for a, b in zip(f, g) ])
    else:
        k = abs(df - dg)

        if df > dg:
            h, f = f[:k], f[k:]
        else:
            h, g = g[:k], g[k:]

        return h + [ a + b for a, b in zip(f, g) ]

def zzX_add(f, g):
    """Add polynomials in Z[X]. """
    if poly_univariate_p(f):
        return zzx_add(f, g)

    if zzX_zero_p(f):
        return g
    if zzX_zero_p(g):
        return f

    df = zzX_degree(f)
    dg = zzX_degree(g)

    if df == dg:
        return zzX_strip([ zzX_add(a, b) for a, b in zip(f, g) ])
    else:
        k = abs(df - dg)

        if df > dg:
            h, f = f[:k], f[k:]
        else:
            h, g = g[:k], g[k:]

        return h + [ zzX_add(a, b) for a, b in zip(f, g) ]

def zzx_sub(f, g):
    """Subtract polynomials in Z[x]. """
    if not g:
        return f
    if not f:
        return zzx_neg(g)

    df = zzx_degree(f)
    dg = zzx_degree(g)

    if df == dg:
        return zzx_strip([ a - b for a, b in zip(f, g) ])
    else:
        k = abs(df - dg)

        if df > dg:
            h, f = f[:k], f[k:]
        else:
            h, g = zzx_neg(g[:k]), g[k:]

        return h + [ a - b for a, b in zip(f, g) ]

def zzX_sub(f, g):
    """Subtract polynomials in Z[x]. """
    if poly_univariate_p(f):
        return zzx_sub(f, g)

    if zzX_zero_p(g):
        return f
    if zzX_zero_p(f):
        return zzX_neg(g)

    df = zzX_degree(f)
    dg = zzX_degree(g)

    if df == dg:
        return zzX_strip([ zzX_sub(a, b) for a, b in zip(f, g) ])
    else:
        k = abs(df - dg)

        if df > dg:
            h, f = f[:k], f[k:]
        else:
            h, g = zzX_neg(g[:k]), g[k:]

        return h + [ zzX_sub(a, b) for a, b in zip(f, g) ]

def zzx_add_mul(f, g, h):
    """Returns f + g*h where f, g, h in Z[x]. """
    return zzx_add(f, zzx_mul(g, h))

def zzX_add_mul(f, g, h):
    """Returns f + g*h where f, g, h in Z[X]. """
    return zzX_add(f, zzX_mul(g, h))

def zzx_sub_mul(f, g, h):
    """Returns f - g*h where f, g, h in Z[x]. """
    return zzx_sub(f, zzx_mul(g, h))

def zzX_sub_mul(f, g, h):
    """Returns f - g*h where f, g, h in Z[X]. """
    return zzX_sub(f, zzX_mul(g, h))

def zzx_mul(f, g):
    """Multiply polynomials in Z[x]. """
    if f == g:
        return zzx_sqr(f)

    if not (f and g):
        return []

    df = zzx_degree(f)
    dg = zzx_degree(g)

    h = []

    for i in xrange(0, df+dg+1):
        coeff = 0

        for j in xrange(max(0, i-dg), min(df, i)+1):
            coeff += f[j]*g[i-j]

        h.append(coeff)

    return h

def zzX_mul(f, g):
    """Multiply polynomials in Z[X]. """
    if poly_univariate_p(f):
        return zzx_mul(f, g)

    if f == g:
        return zzX_sqr(f)

    if zzX_zero_p(f):
        return f
    if zzX_zero_p(g):
        return g

    df = zzX_degree(f)
    dg = zzX_degree(g)

    h, l = [], poly_level(f)-1

    for i in xrange(0, df+dg+1):
        coeff = zzX_zero(l)

        for j in xrange(max(0, i-dg), min(df, i)+1):
            coeff = zzX_add(coeff, zzX_mul(f[j], g[i-j]))

        h.append(coeff)

    return h

def zzx_sqr(f):
    """Square polynomials in Z[x]. """
    df, h = zzx_degree(f), []

    for i in xrange(0, 2*df+1):
        coeff = INT_ZERO

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

        h.append(coeff)

    return h

def zzX_sqr(f):
    """Square polynomials in Z[X]. """
    if poly_univariate_p(f):
        return zzx_sqr(f)

    if zzX_zero_p(f):
        return f

    df = zzX_degree(f)
    l = poly_level(f)-1

    h = []

    for i in xrange(0, 2*df+1):
        coeff = zzX_zero(l)

        jmin = max(0, i-df)
        jmax = min(i, df)

        n = jmax - jmin + 1

        jmax = jmin + n // 2 - 1

        for j in xrange(jmin, jmax+1):
            coeff = zzX_add(coeff, zzX_mul(f[j], f[i-j]))

        coeff = zzX_mul_const(coeff, 2)

        if n & 1:
            elem = zzX_sqr(f[jmax+1])
            coeff = zzX_add(coeff, elem)

        h.append(coeff)

    return h

def zzx_pow(f, n):
    """Raise f to the n-th power in Z[x]. """
    if not n:
        return [INT_ONE]
    if n == 1 or f == [] or f == [1]:
        return f

    g = [INT_ONE]

    while True:
        n, m = n//2, n

        if m & 1:
            g = zzx_mul(g, f)

            if n == 0:
                break

        f = zzx_sqr(f)

    return g

def zzX_pow(f, n):
    """Raise f to the n-th power in Z[X]. """
    if poly_univariate_p(f):
        return zzx_pow(f, n)
    if not n:
        return zzX_const_of(f, 1)
    if n == 1 or zzX_zero_p(f) or zzX_one_p(f):
        return f

    g = zzX_const_of(f, 1)

    while True:
        n, m = n//2, n

        if m & 1:
            g = zzX_mul(g, f)

            if n == 0:
                break

        f = zzX_sqr(f)

    return g

def zzx_expand(*polys):
    """Multiply together several polynomials in Z[x]. """
    f = polys[0]

    for g in polys[1:]:
        f = zzx_mul(f, g)

    return f

def zzX_expand(*polys):
    """Multiply together several polynomials in Z[X]. """
    f = polys[0]

    for g in polys[1:]:
        f = zzX_mul(f, g)

    return f

def zzx_div(f, g):
    """Returns quotient and remainder in Z[x]. """
    df = zzx_degree(f)
    dg = zzx_degree(g)

    if not g:
        raise ZeroDivisionError("polynomial division")
    elif df < dg:
        return [], f

    q, r = [], f

    while True:
        dr = zzx_degree(r)

        if dr < dg:
            break

        lc_r = poly_LC(r)
        lc_g = poly_LC(g)

        if lc_r % lc_g != 0:
            break

        c, k = lc_r // lc_g, dr - dg

        q = zzx_add_term(q, c, k)
        h = zzx_mul_term(g, c, k)
        r = zzx_sub(r, h)

    return q, r

def zzX_div(f, g):
    """Returns quotient and remainder in Z[X]. """
    if poly_univariate_p(f):
        return zzx_div(f, g)

    df = zzX_degree(f)
    dg = zzX_degree(g)

    if dg < 0:
        raise ZeroDivisionError("polynomial division")

    q, r = zzX_zero_of(f), f

    if df < dg:
        return q, r

    while True:
        dr = zzX_degree(r)

        if dr < dg:
            break

        lc_r = poly_LC(r)
        lc_g = poly_LC(g)

        c, R = zzX_div(lc_r, lc_g)

        if not zzX_zero_p(R):
            break

        k = dr - dg

        q = zzX_add_term(q, c, k)
        h = zzX_mul_term(g, c, k)
        r = zzX_sub(r, h)

    return q, r

def zzx_quo(f, g):
    """Returns polynomial remainder in Z[x]. """
    return zzx_div(f, g)[0]

def zzX_quo(f, g):
    """Returns polynomial remainder in Z[X]. """
    return zzX_div(f, g)[0]

def zzx_rem(f, g):
    """Returns polynomial remainder in Z[x]. """
    return zzx_div(f, g)[1]

def zzX_rem(f, g):
    """Returns polynomial remainder in Z[X]. """
    return zzX_div(f, g)[1]

def zzx_max_norm(f):
    """Returns maximum norm of a polynomial in Z[x]. """
    if not f:
        return INT_ZERO
    else:
        return max(zzx_abs(f))

def zzX_max_norm(f):
    """Returns maximum norm of a polynomial in Z[X]. """
    if poly_univariate_p(f):
        return zzx_max_norm(f)
    else:
        return max([ zzX_max_norm(coeff) for coeff in f ])

def zzx_l1_norm(f):
    """Returns l1 norm of a polynomial in Z[x]. """
    if not f:
        return INT_ZERO
    else:
        return sum(zzx_abs(f))

def zzX_l1_norm(f):
    """Returns l1 norm of a polynomial in Z[X]. """
    if poly_univariate_p(f):
        return zzx_l1_norm(f)
    else:
        return sum([ zzX_l1_norm(coeff) for coeff in f ])

def zzx_mignotte_bound(f):
    """Mignotte bound for univariate polynomials. """
    a = zzx_max_norm(f)
    b = abs(poly_LC(f))
    n = zzx_degree(f)

    return INT_TYPE(isqrt(n+1))*2**n*a*b

def zzX_mignotte_bound(f):
    """Mignotte bound for multivariate polynomials. """
    a = zzX_max_norm(f)
    b = abs(zzX_zz_LC(f))
    n = sum(zzX_degree_all(f))

    return INT_TYPE(isqrt(n+1))*2**n*a*b

def zzx_diff(f, m=1):
    """m-th order derivative of a polynomial in Z[x]. """
    n, deriv = zzx_degree(f), []

    n = zzX_degree(f)

    if n < m:
        return []

    deriv, c = [], INT_ONE

    for i in xrange(0, m):
        c, n = c*n, n-1

    for coeff in f[:-m]:
        deriv.append(coeff*c)
        c, n = n*(c//(n+m)), n-1

    return deriv

def zzX_diff(f, m=1):
    """m-th order derivative in x_0 of a polynomial in Z[X]. """
    if poly_univariate_p(f):
        return zzx_diff(f, m)

    n = zzX_degree(f)

    if n < m:
        return zzX_zero_of(f)

    deriv, c = [], INT_ONE

    for i in xrange(0, m):
        c, n = c*n, n-1

    for coeff in f[:-m]:
        h = zzX_mul_const(coeff, c)
        c, n = n*(c//(n+m)), n-1
        deriv.append(h)

    return deriv

def zzX_diff_for(f, k, m=1):
    """m-th order derivative in x_k of a polynomial in Z[X]. """
    if k < 0:
        k += poly_level(f) + 1

    if k == 1:
        return zzX_diff(f, m)

    def rec_diff(g, l):
        if l == k:
            return zzX_diff(g, m)
        else:
            return zzX_strip([ rec_diff(coeff, l+1) for coeff in g ])

    return rec_diff(f, 1)

def zzx_eval(f, x):
    """Evaluate f(x) in Z[x] using Horner scheme. """
    result = INT_ZERO

    if not x:
        return poly_TC(f)

    for a in f:
        result *= x
        result += a

    return result

def zzX_eval(f, x):
    """Evaluate f in Z[X] in variable x_1, using Horner scheme. """
    if hasattr(x, '__iter__'):
        return zzX_eval_list(f, x)

    if poly_univariate_p(f):
        return zzx_eval(f, x)

    if not x:
        return poly_TC(f)

    result = poly_LC(f)

    for coeff in f[1:]:
        result = zzX_mul_const(result, x)
        result = zzX_add(result, coeff)

    return result

def zzX_eval_for(f, k, x):
    """Evaluate f in Z[X] in variable x_k, using Horner scheme. """
    if k < 0:
        k += poly_level(f) + 1

    if k == 1:
        return zzX_eval(f, x)

    def rec_eval(g, l):
        if l == k:
            return zzX_eval(g, x)
        else:
            return zzX_strip([ rec_eval(coeff, l+1) for coeff in g ])

    return rec_eval(f, 1)

def zzX_eval_list(f, A):
    """Evaluate f(x_j,...,x_n) polynomial in Z[X]. """
    def rec_eval(g, l, L):
        if l == L:
            return zzx_eval(g, A[-1])
        else:
            h = [ rec_eval(h, l+1, L) for h in g ]

            if l <= L - len(A):
                return h
            else:
                return zzx_eval(h, A[-L+l-1])

    if not A:
        return f

    L = poly_level(f)

    if zzX_zero_p(f):
        return zzX_zero(L - len(A))

    e = rec_eval(f, 1, L)

    if L == len(A):
        return e
    else:
        return zzX_strip(e)

def zzX_diff_eval(f, k, m, x):
    """Differentiate and evaluate f in Z[X] in variable x_k. """
    if k < 0:
        k += poly_level(f) + 1

    if k == 1:
        return zzX_eval(zzX_diff(f, m), x)

    def rec_eval(g, l):
        if l == k:
            return zzX_eval(zzX_diff(g, m), x)
        else:
            return zzX_strip([ rec_eval(coeff, l+1) for coeff in g ])

    return rec_eval(f, 1)

def zzx_trunc(f, p):
    """Reduce Z[x] polynomial modulo integer p. """
    g = []

    for coeff in f:
        coeff %= p

        if coeff > p // 2:
            g.append(coeff - p)
        else:
            g.append(coeff)

    return zzx_strip(g)

def zzX_trunc(f, p):
    """Reduce Z[X] polynomial modulo polynomial p. """
    return zzX_strip([ zzX_rem(g, p) for g in f ])

def zzX_zz_trunc(f, p):
    """Reduce Z[X] polynomial modulo integer p. """
    if poly_univariate_p(f):
        return zzx_trunc(f, p)
    else:
        return zzX_strip([ zzX_zz_trunc(g, p) for g in f ])

def zzx_content(f):
    """Returns integer GCD of coefficients in Z[x]. """
    cont = INT_ZERO

    for coeff in f:
        cont = igcd(cont, coeff)

        if cont == 1:
            break

    return cont

def zzX_content(f):
    """Returns polynomial GCD of coefficients in Z[X]. """
    cont = poly_LC(f)

    for g in f[1:]:
        cont = zzX_gcd(cont, g)

        if zzX_one_p(cont):
            break

    return cont

def zzX_zz_content(f):
    """Returns integer GCD of coefficients in Z[X]. """
    if poly_univariate_p(f):
        return zzx_content(f)

    cont = INT_ZERO

    for g in f:
        gc = zzX_zz_content(g)
        cont = igcd(cont, gc)

        if cont == 1:
            break

    return cont

def zzx_primitive(f):
    """Divides all coefficients by integer content. """
    cont = zzx_content(f)

    if cont == 1:
        return 1, f
    else:
        return cont, [ coeff // cont for coeff in f ]

def zzX_primitive(f):
    """Divides all coefficients by polynomial content. """
    cont = zzX_content(f)

    if zzX_one_p(cont):
        return cont, f
    else:
        return cont, [ zzX_quo(coeff, cont) for coeff in f ]

def zzX_zz_primitive(f):
    """Divides all coefficients by integer content in Z[X]. """
    cont = zzX_zz_content(f)

    if cont == 1:
        return 1, f
    else:
        return cont, zzX_quo_const(f, cont)

def zzx_sqf_part(f):
    """Returns square-free part of a polynomial in Z[x]. """
    quo = zzx_quo(f, zzx_gcd(f, zzx_diff(f)))
    return zzx_primitive(quo)[1]

def zzX_sqf_part(f):
    """Returns square-free part of a polynomial in Z[X]. """
    quo = zzX_quo(f, zzX_gcd(f, zzX_diff(f)))
    return zzX_primitive(quo)[1]

def zzx_sqf_p(f):
    """Returns True if f is a square-free polynomial in Z[x]. """
    return zzx_one_p(zzx_gcd(zzx_primitive(f)[1], zzx_diff(f)))

def zzX_sqf_p(f):
    """Returns True if f is a square-free polynomial in Z[X]. """
    return zzX_one_p(zzX_gcd(zzX_primitive(f)[1], zzX_diff(f)))

def zzx_gcd(f, g, **flags):
    """Returns polynomial GCD in Z[x]. """
    return zzx_cofactors(f, g, **flags)[0]

def zzx_cofactors(f, g, **flags):
    """Returns polynomial GCD and its co-factors in Z[x]. """
    return zzx_heu_gcd(f, g, **flags)

def zzX_gcd(f, g, **flags):
    """Returns polynomial GCD in Z[X]. """
    return zzX_cofactors(f, g, **flags)[0]

def zzX_cofactors(f, g, **flags):
    """Returns polynomial GCD and its co-factors in Z[X]. """
    if poly_univariate_p(f):
        return zzx_heu_gcd(f, g, **flags)

    if not flags.get('reduced', True):
        return zzX_heu_gcd(f, g, **flags)
    else:
        K, (f, g) = zzX_multi_reduce(f, g)

        return [ zzX_compose_term(h, K)
            for h in zzX_heu_gcd(f, g, **flags) ]

def zzx_heu_gcd(f, g, **flags):
    """Heuristic polynomial GCD over Z[x].

       Given univariate polynomials f and g over Z[x], returns their GCD
       and cofactors, i.e. polynomials h, cff and cfg such that:

              h = gcd(f, g), cff = quo(f, h) and cfg = quo(g, h)

       The algorithm is purely heuristic which means it may fail to compute
       the GCD. This will be signaled by raising an exception. In this case
       you will need to switch to another GCD method.

       The algorithm computes the polynomial GCD by evaluating polynomials
       f and g at certain points and computing (fast) integer GCD of those
       evaluations. The polynomial GCD is recovered from the integer image
       by interpolation.  The final step is to verify if the result is the
       correct GCD. This gives cofactors of input polynomials as a side
       effect (see zzX_cofactors).

       For more details on the implemented algorithm refer to:

       [1] Hsin-Chao Liao, R. Fateman, Evaluation of the heuristic polynomial
           GCD, International Symposium on Symbolic and Algebraic Computation
           (ISSAC), ACM Press, Montreal, Quebec, Canada, 1995, pp. 240--247

    """
    def interpolate(h, x):
        f = []

        while h:
            g = h % x

            if g > x // 2:
                g -= x

            f.insert(0, g)
            h = (h-g) // x

        return f

    def finalize(h, cff, cfg, gcd):
        h = zzx_mul_const(h, gcd)
        return h, cff, cfg

    if not (f or g):
        return [], [], []
    elif not f:
        return g, [], [1]
    elif not g:
        return f, [1], []

    df = zzx_degree(f)
    dg = zzx_degree(g)

    cf = zzx_content(f)
    cg = zzx_content(g)

    gcd = igcd(cf, cg)

    f = [ c // gcd for c in f ]
    g = [ c // gcd for c in g ]

    if df == 0 or dg == 0:
        return [gcd], f, g

    f_norm = zzx_max_norm(f)
    g_norm = zzx_max_norm(g)

    B = 2*min(f_norm, g_norm) + 29

    x = max(min(B, 99*INT_TYPE(isqrt(B))),
            2*min(f_norm // abs(poly_LC(f)),
                  g_norm // abs(poly_LC(g))) + 2)

    for i in xrange(0, 6):
        ff = zzx_eval(f, x)
        gg = zzx_eval(g, x)

        if ff and gg:
            h = igcd(ff, gg)

            cff = ff // h
            cfg = gg // h

            h = interpolate(h, x)
            h = zzx_primitive(h)[1]

            cff_, r = zzx_div(f, h)

            if not r:
                cfg_, r = zzx_div(g, h)

                if not r:
                    return finalize(h, cff_, cfg_, gcd)

            cff = interpolate(cff, x)

            h, r = zzx_div(f, cff)

            if not r:
                cfg_, r = zzx_div(g, h)

                if not r:
                    return finalize(h, cff, cfg_, gcd)

            cfg = interpolate(cfg, x)

            h, r = zzx_div(g, cfg)

            if not r:
                cff_, r = zzx_div(f, h)

                if not r:
                    return finalize(h, cff_, cfg, gcd)

        x = INT_TYPE(2.7319*x*isqrt(isqrt(x)))

    raise HeuristicGCDFailed('no luck')

def zzX_heu_gcd(f, g, **flags):
    """Heuristic polynomial GCD in Z[X].

       Given univariate polynomials f and g in Z[X], returns their GCD
       and cofactors, i.e. polynomials h, cff and cfg such that:

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
       the input polynomials as a side effect (see zzX_cofactors).

       For more details on the implemented algorithm refer to:

       [1] Hsin-Chao Liao, R. Fateman, Evaluation of the heuristic polynomial
           GCD, International Symposium on Symbolic and Algebraic Computation
           (ISSAC), ACM Press, Montreal, Quebec, Canada, 1995, pp. 240--247
    """
    if poly_univariate_p(f):
        return zzx_heu_gcd(f, g, **flags)

    def interpolate(h, x):
        f = []

        while not zzX_zero_p(h):
            g = zzX_zz_trunc(h, x)
            f.insert(0, g)
            h = zzX_sub(h, g)
            h = zzX_quo_const(h, x)

        return f

    def finalize(h, cff, cfg, gcd):
        if zzX_zz_LC(h) > 0:
            h = zzX_mul_const(h, gcd)
        else:
            h = zzX_mul_const(h, -gcd)
            cff = zzX_neg(cff)
            cfg = zzX_neg(cfg)

        return h, cff, cfg

    zero_f = zzX_zero_p(f)
    zero_g = zzX_zero_p(g)

    l = poly_level(f)
    z = zzX_zero(l)

    if zero_f and zero_g:
        return z, z, z
    elif zero_f:
        return g, z, zzX_const(l, 1)
    elif zero_g:
        return f, zzX_const(l, 1), z

    df = zzX_degree(f)
    dg = zzX_degree(g)

    cf = zzX_zz_content(f)
    cg = zzX_zz_content(g)

    gcd = igcd(cf, cg)

    f = zzX_quo_const(f, gcd)
    g = zzX_quo_const(g, gcd)

    f_norm = zzX_max_norm(f)
    g_norm = zzX_max_norm(g)

    B = 2*min(f_norm, g_norm) + 29

    x = max(min(B, 99*INT_TYPE(isqrt(B))),
            2*min(f_norm // abs(zzX_zz_LC(f)),
                  g_norm // abs(zzX_zz_LC(g))) + 2)

    for i in xrange(0, 6):
        ff = zzX_eval(f, x)
        gg = zzX_eval(g, x)

        if not (zzX_zero_p(ff) or zzX_zero_p(gg)):
            h, cff, cfg = zzX_heu_gcd(ff, gg, **flags)

            h = interpolate(h, x)
            h = zzX_zz_primitive(h)[1]

            cff_, r = zzX_div(f, h)

            if zzX_zero_p(r):
                cfg_, r = zzX_div(g, h)

                if zzX_zero_p(r):
                    return finalize(h, cff_, cfg_, gcd)

            cff = interpolate(cff, x)

            h, r = zzX_div(f, cff)

            if zzX_zero_p(r):
                cfg_, r = zzX_div(g, h)

                if zzX_zero_p(r):
                    return finalize(h, cff, cfg_, gcd)

            cfg = interpolate(cfg, x)

            h, r = zzX_div(g, cfg)

            if zzX_zero_p(r):
                cff_, r = zzX_div(f, h)

                if zzX_zero_p(r):
                    return finalize(h, cff_, cfg, gcd)

        x = INT_TYPE(2.7319*x*isqrt(isqrt(x)))

    raise HeuristicGCDFailed('no luck')

# XXX: this is completely broken
def zzx_mod_gcd(f, g, **flags):
    """Modular small primes polynomial GCD over Z[x].

       Given univariate polynomials f and g over Z[x], returns their
       GCD and cofactors, i.e. polynomials h, cff and cfg such that:

              h = gcd(f, g), cff = quo(f, h) and cfg = quo(g, h)

       The algorithm uses modular small primes approach. It works by
       computing several GF(p)[x] GCDs for a set of randomly chosen
       primes and uses Chinese Remainder Theorem to recover the GCD
       over Z[x] from its images.

       The algorithm is probabilistic which means it never fails,
       however its running time depends on the number of unlucky
       primes chosen for computing GF(p)[x] images.

       For more details on the implemented algorithm refer to:

       [1] J. von zur Gathen, J. Gerhard, Modern Computer Algebra,
           First Edition, Cambridge University Press, 1999, pp. 158

    """
    if not (f or g):
        return [], [], []
    elif not f:
        return g, [], [1]
    elif not g:
        return f, [1], []

    n = zzx_degree(f)
    m = zzx_degree(g)

    cf = zzx_content(f)
    cg = zzx_content(g)

    gcd = igcd(cf, cg)

    f = [ c // gcd for c in f ]
    g = [ c // gcd for c in g ]

    if n == 0 or m == 0:
        return [gcd], f, g

    A = max(zzx_abs(f) + zzx_abs(g))
    b = igcd(poly_LC(f), poly_LC(g))

    B = int(ceil(2**n*A*b*int(sqrt(n + 1))))
    k = int(ceil(2*b*log((n + 1)**n*A**(2*n), 2)))
    l = int(ceil(log(2*B + 1, 2)))

    prime_max = max(int(ceil(2*k*log(k))), 51)

    while True:
        while True:
            primes  = set([])
            unlucky = set([])

            ff, gg, hh = {}, {}, {}

            while len(primes) < l:
                p = randprime(3, prime_max+1)

                if (p in primes) or (b % p == 0):
                    continue

                F = gf_from_int_poly(f, p)
                G = gf_from_int_poly(g, p)

                H = gf_gcd(F, G, p)

                primes.add(p)

                ff[p] = F
                gg[p] = G
                hh[p] = H

            e = min([ gf_degree(h) for h in hh.itervalues() ])

            for p in set(primes):
                if gf_degree(hh[p]) != e:
                    primes.remove(p)
                    unlucky.add(p)

                    del ff[p]
                    del gg[p]
                    del hh[p]

            if len(primes) < l // 2:
                continue

            while len(primes) < l:
                p = randprime(3, prime_max+1)

                if (p in primes) or (p in unlucky) or (b % p == 0):
                    continue

                F = gf_from_int_poly(f, p)
                G = gf_from_int_poly(g, p)

                H = gf_gcd(F, G, p)

                if gf_degree(H) != e:
                    unlucky.add(p)
                else:
                    primes.add(p)

                    ff[p] = F
                    gg[p] = G
                    hh[p] = H

            break

        fff, ggg = {}, {}

        for p in primes:
            fff[p] = gf_quo(ff[p], hh[p], p)
            ggg[p] = gf_quo(gg[p], hh[p], p)

        F, G, H = [], [], []

        crt_mm, crt_e, crt_s = crt1(primes)

        for i in xrange(0, e + 1):
            C = [ b * poly_nth(hh[p], i) for p in primes ]
            c = crt2(primes, C, crt_mm, crt_e, crt_s, True)

            H.insert(0, c)

        H = zzx_strip(H)

        for i in xrange(0, zzx_degree(f) - e + 1):
            C = [ poly_nth(fff[p], i) for p in primes ]
            c = crt2(primes, C, crt_mm, crt_e, crt_s, True)

            F.insert(0, c)

        for i in xrange(0, zzx_degree(g) - e + 1):
            C = [ poly_nth(ggg[p], i) for p in primes ]
            c = crt2(primes, C, crt_mm, crt_e, crt_s, True)

            G.insert(0, c)

        H_norm = zzx_l1_norm(H)

        F_norm = zzx_l1_norm(F)
        G_norm = zzx_l1_norm(G)

        if H_norm*F_norm <= B and H_norm*G_norm <= B:
            break

    return zzx_mul_const(H, gcd), F, G

def zzx_hensel_step(m, f, g, h, s, t):
    """One step in Hensel lifting.

       Given positive integer m and Z[x] polynomials f, g, h, s and t such that:

        [1] f == g*h (mod m)
        [2] s*g + t*h == 1 (mod m)

        [3] lc(f) not a zero divisor (mod m)
        [4] lc(h) == 1

        [5] deg(f) == deg(g) + deg(h)
        [6] deg(s) < deg(h)
        [7] deg(t) < deg(g)

       returns polynomials G, H, S and T, such that:

        [A] f == G*H (mod m**2)
        [B] S*G + T**H == 1 (mod m**2)

       For more details on the implemented algorithm refer to:

       [1] J. von zur Gathen, J. Gerhard, Modern Computer Algebra,
           First Edition, Cambridge University Press, 1999, pp. 418

    """
    M = m**2

    e = zzx_sub_mul(f, g, h)
    e = zzx_trunc(e, M)

    q, r = zzx_div(zzx_mul(s, e), h)

    q = zzx_trunc(q, M)
    r = zzx_trunc(r, M)

    u = zzx_add(zzx_mul(t, e), zzx_mul(q, g))
    G = zzx_trunc(zzx_add(g, u), M)
    H = zzx_trunc(zzx_add(h, r), M)

    u = zzx_add(zzx_mul(s, G), zzx_mul(t, H))
    b = zzx_trunc(zzx_sub(u, [1]), M)

    c, d = zzx_div(zzx_mul(s, b), H)

    c = zzx_trunc(c, M)
    d = zzx_trunc(d, M)

    u = zzx_add(zzx_mul(t, b), zzx_mul(c, G))
    S = zzx_trunc(zzx_sub(s, d), M)
    T = zzx_trunc(zzx_sub(t, u), M)

    return G, H, S, T

def zzx_hensel_lift(p, f, f_list, l):
    """Multifactor Hensel lifting.

       Given a prime p, polynomial f over Z[x] such that lc(f) is a
       unit modulo p,  monic pair-wise coprime polynomials f_i over
       Z[x] satisfying:

                    f = lc(f) f_1 ... f_r (mod p)

       and a positive integer l, returns a list of monic polynomials
       F_1, F_2, ..., F_r satisfying:

                    f = lc(f) F_1 ... F_r (mod p**l)

                    F_i = f_i (mod p), i = 1..r

       For more details on the implemented algorithm refer to:

       [1] J. von zur Gathen, J. Gerhard, Modern Computer Algebra,
           First Edition, Cambridge University Press, 1999, pp. 424

    """
    r = len(f_list)
    lc = poly_LC(f)

    if r == 1:
        F = zzx_mul_const(f, igcdex(lc, p**l)[0])
        return [ zzx_trunc(F, p**l) ]

    m = p
    k = int(r // 2)
    d = int(ceil(log(l, 2)))

    g = gf_from_int_poly([lc], p)

    for f_i in f_list[:k]:
        g = gf_mul(g, gf_from_int_poly(f_i, p), p)

    h = gf_from_int_poly(f_list[k], p)

    for f_i in f_list[k+1:]:
        h = gf_mul(h, gf_from_int_poly(f_i, p), p)

    s, t, _ = gf_gcdex(g, h, p)

    g = gf_to_int_poly(g, p)
    h = gf_to_int_poly(h, p)
    s = gf_to_int_poly(s, p)
    t = gf_to_int_poly(t, p)

    for _ in range(1, d+1):
        (g, h, s, t), m = zzx_hensel_step(m, f, g, h, s, t), m**2

    return zzx_hensel_lift(p, g, f_list[:k], l) \
         + zzx_hensel_lift(p, h, f_list[k:], l)

def zzx_zassenhaus(f):
    """Factor primitive square-free polynomials in Z[x]. """
    n = zzx_degree(f)

    if n == 1:
        return [f]

    A = zzx_max_norm(f)
    b = poly_LC(f)
    B = abs(int(sqrt(n+1))*2**n*A*b)
    C = (n+1)**(2*n)*A**(2*n-1)
    gamma = int(ceil(2*log(C, 2)))
    prime_max = int(2*gamma*log(gamma))

    for p in xrange(3, prime_max+1):
        if not isprime(p) or b % p == 0:
            continue

        F = gf_from_int_poly(f, p)

        if gf_sqf_p(F, p):
            break

    l = int(ceil(log(2*B + 1, p)))

    modular = []

    for ff in gf_factor_sqf(F, p)[1]:
        modular.append(gf_to_int_poly(ff, p))

    g = zzx_hensel_lift(p, f, modular, l)

    T = set(range(len(g)))
    factors, s = [], 1

    while 2*s <= len(T):
        for S in subsets(T, s):
            G, H = [b], [b]

            S = set(S)

            for i in S:
                G = zzx_mul(G, g[i])
            for i in T-S:
                H = zzx_mul(H, g[i])

            G = zzx_trunc(G, p**l)
            H = zzx_trunc(H, p**l)

            G_norm = zzx_l1_norm(G)
            H_norm = zzx_l1_norm(H)

            if G_norm*H_norm <= B:
                T = T - S

                G = zzx_primitive(G)[1]
                f = zzx_primitive(H)[1]

                factors.append(G)
                b = poly_LC(f)

                break
        else:
            s += 1

    return factors + [f]

def zzx_eisenstein(f):
    """Eisenstein's irreducibility criterion. """
    lc = poly_LC(f)
    tc = poly_TC(f)

    e_fc = zzx_content(f[1:])

    if not e_fc:
        return

    e_ff = factorint(e_fc)

    for p in e_ff.iterkeys():
        if (lc % p) and (tc % p**2):
            return True

def zzx_factor(f, **flags):
    """Factor (non square-free) polynomials in Z[x].

       Given a univariate polynomial f in Z[x] computes its complete
       factorization f_1, ..., f_n into irreducibles over integers:

                   f = content(f) f_1**k_1 ... f_n**k_n

       The factorization is computed by reducing the input polynomial
       into a primitive square-free polynomial and factoring it using
       Zassenhaus algorithm. Trial division is used to recover the
       multiplicities of factors.

       The result is returned as a tuple consisting of:

                 (content(f), [(f_1, k_1), ..., (f_n, k_n))

       Consider polynomial f = 2*x**4 - 2:

       >>> from sympy.polys.integerpolys import zzx_factor
       >>> zzx_factor([2, 0, 0, 0, -2])
       (2, [([1, -1], 1), ([1, 1], 1), ([1, 0, 1], 1)])

       In result we got the following factorization:

                    f = 2 (x - 1) (x + 1) (x**2 + 1)

       Note that this is a complete factorization over integers,
       however over Gaussian integers we can factor the last term.

       By default, polynomials x**n - 1 and x**n + 1 are factored
       using cyclotomic decomposition to speedup computations. To
       disable this behaviour set cyclotomic=False.

       For more details on the implemented algorithm refer to:

       [1] J. von zur Gathen, J. Gerhard, Modern Computer Algebra,
           First Edition, Cambridge University Press, 1999, pp. 427
    """
    cont, g = zzx_primitive(f)

    n = zzx_degree(g)

    if n <= 0:
        return cont, []

    if poly_LC(g) < 0:
        cont, g = -cont, zzx_neg(g)

    if n == 1 or zzx_eisenstein(g):
        return cont, [(g, 1)]

    g = zzx_sqf_part(g)
    H, factors = None, []

    if flags.get('cyclotomic', True):
        H = zzx_cyclotomic_factor(g)

    if H is None:
        H = zzx_zassenhaus(g)

    for h in H:
        k = 0

        while True:
            q, r = zzx_div(f, h)

            if not r:
                f, k = q, k+1
            else:
                break

        factors.append((h, k))

    def compare((f_a, e_a), (f_b, e_b)):
        i = len(f_a) - len(f_b)

        if not i:
            j = e_a - e_b

            if not j:
                return cmp(f_a, f_b)
            else:
                return j
        else:
            return i

    return cont, sorted(factors, compare)

def zzx_factor_sqf(f, **flags):
    """Factor square-free (non-primitive) polyomials in Z[x].  """
    cont, g = zzx_primitive(f)

    n = zzx_degree(g)

    if n <= 0:
        return cont, []

    if poly_LC(g) < 0:
        cont, g = -cont, zzx_neg(g)

    if n == 1 or zzx_eisenstein(g):
        return cont, [(g, 1)]

    factors = []

    if flags.get('cyclotomic', True):
        factors = zzx_cyclotomic_factor(g)

    if factors is None:
        factors = zzx_zassenhaus(g)

    def compare(f_a, f_b):
        i = len(f_a) - len(f_b)

        if not i:
            return cmp(f_a, f_b)
        else:
            return i

    return cont, sorted(factors, compare)

def zzx_cyclotomic_factor(f):
    """Efficiently factor polynomials x**n - 1 and x**n + 1 in Z[x].

       Given a univariate polynomial f in Z[x] returns a list of factors
       of f, provided that f is in form x**n - 1 or x**n + 1 for n >= 1.
       Otherwise returns None.

       Factorization is performed using using cyclotomic decomposition
       of f, which makes this method much faster that any other direct
       factorization approach (e.g. Zassenhaus's).

       For more details on the implemented algorithm refer to:

       [1] Eric W. Weisstein, Cyclotomic Polynomial, From MathWorld - A Wolfram
           Web Resource, http://mathworld.wolfram.com/CyclotomicPolynomial.html
    """
    lc_f, tc_f = poly_LC(f), poly_TC(f)

    if zzx_degree(f) <= 0:
        return None

    if lc_f != 1 or tc_f not in [-1, 1]:
        return None

    if any([ bool(cf) for cf in f[1:-1] ]):
        return None

    def decompose(n):
        H = [[1,-1]]

        for p, k in factorint(n).iteritems():
            Q = [ zzx_quo(zzx_compose_term(h, p), h) for h in H ]
            H.extend(Q)

            for i in xrange(1, k):
                Q = [ zzx_compose_term(q, p) for q in Q ]
                H.extend(Q)

        return H

    n = zzx_degree(f)
    F = decompose(n)

    if tc_f != 1:
        return F
    else:
        H = []

        for h in decompose(2*n):
            if h not in F:
                H.append(h)

        return H

def zzX_wang_non_divisors(E, cu, cv):
    """EEZ: Compute a set of valid divisors.  """
    result = [ cu*cv ]

    for q in E:
        q = abs(q)

        for r in reversed(result):
            while r != 1:
                r = igcd(r, q)
                q = q // r

            if q == 1:
                return None

        result.append(q)

    return result[1:]

def zzX_wang_test_points(f, V, cv, A):
    """EEZ: Test evaluation points for suitability. """
    if not zzX_eval(poly_LC(f), A):
        return None

    U = zzX_eval(f, A)

    if not zzx_sqf_p(U):
        return None

    cu, u = zzx_primitive(U)

    if poly_LC(u) < 0:
        cu, u = -cu, zzx_neg(u)

    E = [ zzX_eval(v, A) for v, _ in V ]
    D = zzX_wang_non_divisors(E, cu, cv)

    if D is not None:
        return cu, u, E, D
    else:
        return None

def zzX_wang_lead_coeffs(f, V, cu, E, H, A):
    """EEZ: Compute correct leading coefficients. """
    l = poly_level(f) - 1
    C, K = [], [0]*len(E)

    for h in H:
        c = zzX_const(l, 1)
        d = poly_LC(h)*cu

        for i in reversed(xrange(len(E))):
            k, e, (v, _) = 0, E[i], V[i]

            while not (d % e):
                d, k = d//e, k+1

            if k != 0:
                c, K[i] = zzX_mul(c, zzX_pow(v, k)), 1

        C.append(c)

    if any([ not k for k in K ]):
        raise ExtraneousFactors

    CC, HH = [], []

    for c, h in zip(C, H):
        d = zzX_eval(c, A)
        lc = poly_LC(h)

        if cu == 1:
            cc = lc//d
        else:
            g = igcd(lc, d)
            d, cc = d//g, lc//g
            h, cu = zzx_mul_const(h, d), cu//d

        c = zzX_mul_const(c, cc)

        CC.append(c)
        HH.append(h)

    if cu == 1:
        return f, HH, CC

    CCC, HHH = [], []

    for c, h in zip(CC, HH):
        CCC.append(zzX_mul_const(c, cu))
        HHH.append(zzX_mul_const(h, cu))

    f = zzX_mul_const(f, cu**(len(H)-1))

    return f, HHH, CCC

def zzX_wang_more_coeffs(f, C, H):
    pass # XXX: to be done

def zzx_diophantine(F, m, p):
    """Solve univariate Diophantine equations. """
    if len(F) == 2:
        a, b = F

        f = gf_from_int_poly(a, p)
        g = gf_from_int_poly(b, p)

        s, t, G = gf_gcdex(g, f, p)

        s = gf_lshift(s, m)
        t = gf_lshift(t, m)

        q, s = gf_div(s, f, p)

        t = gf_add_mul(t, q, g, p)

        s = gf_to_int_poly(s, p)
        t = gf_to_int_poly(t, p)

        result = [s, t]
    else:
        G = [F[-1]]

        for f in reversed(F[1:-1]):
            G.insert(0, zzx_mul(f, G[0]))

        S, T = [], [[1]]

        for f, g in zip(F, G):
            t, s = zzX_diophantine([g, f], T[-1], [], 0, p)
            T.append(t)
            S.append(s)

        result, S = [], S + [T[-1]]

        for s, f in zip(S, F):
            s = gf_from_int_poly(s, p)
            f = gf_from_int_poly(f, p)

            r = gf_rem(gf_lshift(s, m), f, p)
            s = gf_to_int_poly(r, p)

            result.append(s)

    return result

def zzX_diophantine(F, c, A, d, p):
    """Solve multivariate Diophantine equations. """
    if not A:
        S = [ [] for _ in F ]
        n = zzx_degree(c)

        for i, coeff in enumerate(c):
            if not coeff:
                continue

            T = zzx_diophantine(F, n-i, p)

            for j, (s, t) in enumerate(zip(S, T)):
                t = zzx_mul_const(t, coeff)
                S[j] = zzx_trunc(zzx_add(s, t), p)
    else:
        n = len(A) + 1
        e = zzX_expand(*F)

        a, A = A[-1], A[:-1]
        B, G = [], []

        for f in F:
            B.append(zzX_quo(e, f))
            G.append(zzX_eval_for(f, n, a))

        C = zzX_eval_for(c, n, a)

        S = zzX_diophantine(G, C, A, d, p)
        S = [ zzX_lift(1, s) for s in S ]

        for s, b in zip(S, B):
            c = zzX_sub_mul(c, s, b)

        c = zzX_zz_trunc(c, p)

        m = zzX_value(n-1, [1, -a])
        M = zzX_const(n, 1)

        for k in xrange(0, d):
            if zzX_zero_p(c):
                break

            M = zzX_mul(M, m)
            C = zzX_diff_eval(c, n, k+1, a)

            if not zzX_zero_p(C):
                C = zzX_quo_const(C, factorial(k+1))
                T = zzX_diophantine(G, C, A, d, p)

                for i, t in enumerate(T):
                    T[i] = zzX_mul(zzX_lift(1, t), M)

                for i, (s, t) in enumerate(zip(S, T)):
                    S[i] = zzX_add(s, t)

                for t, b in zip(T, B):
                    c = zzX_sub_mul(c, t, b)

                c = zzX_zz_trunc(c, p)

        S = [ zzX_zz_trunc(s, p) for s in S ]

    return S

def zzX_wang_hensel_lifting(f, H, LC, A, p):
    """EEZ: Parallel Hensel lifting algorithm. """
    U, n = [f], len(A)+1

    H = deepcopy(H)

    for i, a in enumerate(reversed(A[1:])):
        u = zzX_eval_for(U[0], n-i, a)
        U.insert(0, zzX_zz_trunc(u, p))

    d = max(zzX_degree_all(f)[1:])

    for j, u, a in zip(xrange(2, n+1), U, A):
        G = deepcopy(H)

        I, J = A[:j-2], A[j-1:]

        for i, (h, lc) in enumerate(zip(H, LC)):
            lc = zzX_zz_trunc(zzX_eval(lc, J), p)
            H[i] = [lc] + zzX_lift(1, h[1:])

        m = zzX_value(j-1, [1, -a])
        M = zzX_const(j, 1)

        c = zzX_sub(u, zzX_expand(*H))

        dj = zzX_degree_for(u, j)

        for k in xrange(0, dj):
            if zzX_zero_p(c):
                break

            M = zzX_mul(M, m)
            C = zzX_diff_eval(c, j, k+1, a)

            if not zzX_zero_p(C):
                C = zzX_quo_const(C, factorial(k+1))
                T = zzX_diophantine(G, C, I, d, p)

                for i, (h, t) in enumerate(zip(H, T)):
                    h = zzX_add_mul(h, zzX_lift(1, t), M)
                    H[i] = zzX_zz_trunc(h, p)

                h = zzX_sub(u, zzX_expand(*H))
                c = zzX_zz_trunc(h, p)

    if zzX_expand(*H) != f:
        raise ExtraneousFactors
    else:
        return H

EEZ_NUM_OK    = 3
EEZ_NUM_TRY   = 20
EEZ_MOD_STEP  = 50

def zzX_wang(f):
    """Factor primitive square-free polynomials in Z[X].

       Given a multivariate polynomial f in Z[x_1,...,x_n], which is
       primitive and square-free in x_1, computes its factorization
       of f into irreducibles over integers.

       The procedure is based on Wang's Enhanced Extended Zassenhaus
       algorithm. The algorithm works by viewing f as a univariate
       polynomial in Z[x_2,...,x_n][x_1], for which an evaluation
       mapping is computed:

                       x_2 -> a_2, ..., x_n -> a_n

       where a_i, for i=2,...,n, are carefully chosen integers. The
       mapping is used to transform f into a univariate polynomial
       in Z[x_1], which can be factored efficiently using Zassenhaus
       algorithm. The last step is to lift univariate factors to
       obtain true multivariate factors. For this purpose a parallel
       Hensel lifting procedure is used.

       For more details on the implemented algorithm refer to:

       [1] P. S. Wang, An Improved Multivariate Polynomial Factoring
           Algorithm, Math. of Computation 32, 1978, pp. 1215--1231

       [2] K. Geddes, S. R. Czapor, G. Labahn, Algorithms for
           Computer Algebra, Springer, 1992, pp. 264--272
    """
    cv, V = zzX_factor(poly_LC(f))

    b = zzX_mignotte_bound(f)
    p = nextprime(b)

    l = poly_level(f)

    bad_points = set([])
    r, mod = None, 5

    while True:
        configs = []

        while len(configs) < EEZ_NUM_OK:
            for i in xrange(EEZ_NUM_TRY):
                A = []

                for j in xrange(0, l-1):
                    A.append(randint(-mod, mod))

                if tuple(A) not in bad_points:
                    bad_points.add(tuple(A))
                else:
                    continue

                R = zzX_wang_test_points(f, V, cv, A)

                if R is not None:
                    cu, u, E, _ = R

                    _, H = zzx_factor_sqf(u)

                    rr = len(H)

                    if r is not None:
                        if rr <= r:
                            if rr < r:
                                configs, r = [], rr
                        else:
                            continue
                    else:
                        r = rr

                    configs.append((u, cu, E, H, A))

                    if len(configs) == EEZ_NUM_OK:
                        break
            else:
                mod += EEZ_MOD_STEP

        if r == 1:
            return 1, [f] # irreducible

        u_norm, u_arg, i = None, 0, 0

        for u, _, _, _, _ in configs:
            _u_norm = zzx_max_norm(u)

            if u_norm is not None:
                if _u_norm < u_norm:
                    u_norm = _u_norm
                    u_arg = i
            else:
                u_norm = _u_norm

            i += 1

        _, cu, E, H, A = configs[u_arg]

        try:
            f, H, LC = zzX_wang_lead_coeffs(f, V, cu, E, H, A)
            factors = zzX_wang_hensel_lifting(f, H, LC, A, p)
        except ExtraneousFactors:
            continue

        negative, H = 0, []

        for h in factors:
            _, h = zzX_zz_primitive(h)

            if zzX_zz_LC(h) < 0:
                h = zzX_neg(h)
                negative += 1

            H.append(h)

        if not (negative % 2):
            return  1, H
        else:
            return -1, H

def zzX_factor(f):
    """Factor (non square-free) polynomials in Z[X].

       Given a multivariate polynomial f in Z[x] computes its complete
       factorization f_1, ..., f_n into irreducibles over integers:

                    f = content(f) f_1**k_1 ... f_n**k_n

       The factorization is computed by reducing the input polynomial
       into a primitive square-free polynomial and factoring it using
       Enhanced Extended Zassenhaus (EEZ) algorithm. Trial division is
       used to recover the multiplicities of factors.

       The result is returned as a tuple consisting of:

                (content(f), [(f_1, k_1), ..., (f_n, k_n))

       Consider polynomial f = 2*(x**2 - y**2):

       >>> from sympy.polys.integerpolys import zzX_factor
       >>> f = [[2], [], [-2, 0, 0]]

       >>> zzX_factor(f)
       (2, [([[1], [-1, 0]], 1), ([[1], [1, 0]], 1)])

       In result we got the following factorization:

                       f = 2 (x - y) (x + y)

       For more details on the implemented algorithm refer to:

       [1] J. von zur Gathen, J. Gerhard, Modern Computer Algebra,
           First Edition, Cambridge University Press, 1999, pp. 427
    """
    if poly_univariate_p(f):
        return zzx_factor(f)

    cont, g = zzX_zz_primitive(f)

    if zzX_degree(g) < 0:
        return cont, []

    if zzX_zz_LC(g) < 0:
        cont, g = -cont, zzX_neg(g)

    G, g = zzX_primitive(g)

    factors = []

    if zzX_degree(g) > 0:
        g = zzX_sqf_part(g)
        s, H = zzX_wang(g)

        cont *= s

        for h in H:
            k = 0

            while True:
                q, r = zzX_div(f, h)

                if zzX_zero_p(r):
                    f, k = q, k+1
                else:
                    break

            factors.append((h, k))

    for g, k in zzX_factor(G)[1]:
        factors.insert(0, ([g], k))

    def compare((f_a, e_a), (f_b, e_b)):
        i = len(f_a) - len(f_b)

        if not i:
            j = e_a - e_b

            if not j:
                return cmp(f_a, f_b)
            else:
                return j
        else:
            return i

    return cont, sorted(factors, compare)


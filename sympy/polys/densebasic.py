"""Basic tools for dense recursive polynomials in `K[x]` or `K[X]`. """

from sympy.utilities import any, all
from sympy.core import igcd, ilcm

from sympy.polys.monomialtools import (
    monomial_key, monomial_min, monomial_div
)

from sympy.polys.groebnertools import sdp_sort

from sympy.utilities import cythonized

import random

def poly_LC(f, K):
    """Returns leading coefficient of `f`. """
    if not f:
        return K.zero
    else:
        return f[0]

def poly_TC(f, K):
    """Returns trailing coefficient of `f`. """
    if not f:
        return K.zero
    else:
        return f[-1]

dup_LC = dmp_LC = poly_LC
dup_TC = dmp_TC = poly_TC

@cythonized("u")
def dmp_ground_LC(f, u, K):
    """Returns ground leading coefficient. """
    while u:
        f = dmp_LC(f, K)
        u -= 1
    return dup_LC(f, K)

@cythonized("u")
def dmp_ground_TC(f, u, K):
    """Returns ground trailing coefficient. """
    while u:
        f = dmp_TC(f, K)
        u -= 1
    return dup_TC(f, K)

@cythonized("u")
def dmp_true_LT(f, u, K):
    """Returns leading term `c * x_1**n_1 ... x_k**n_k`. """
    monom = []

    while u:
        monom.append(len(f) - 1)
        f, u = f[0], u - 1

    if not f:
        monom.append(0)
    else:
        monom.append(len(f) - 1)

    return tuple(monom), dup_LC(f, K)

def dup_degree(f):
    """Returns leading degree of `f` in `K[x]`. """
    return len(f) - 1

@cythonized("u")
def dmp_degree(f, u):
    """Returns leading degree of `f` in `x_0` in `K[X]`. """
    if dmp_zero_p(f, u):
        return -1
    else:
        return len(f) - 1

@cythonized("v,i,j")
def _rec_degree_in(g, v, i, j):
    """XXX"""
    if i == j:
        return dmp_degree(g, v)

    v, i = v-1, i+1

    return max([ _rec_degree_in(c, v, i, j) for c in g ])

@cythonized("j,u")
def dmp_degree_in(f, j, u):
    """Returns leading degree of `f` in `x_j` in `K[X]`. """
    if not j:
        return dmp_degree(f, u)
    if j < 0 or j > u:
        raise IndexError("-%s <= j < %s expected, got %s" % (u, u, j))

    return _rec_degree_in(f, u, 0, j)

@cythonized("v,i")
def _rec_degree_list(g, v, i, degs):
    """XXX"""
    degs[i] = max(degs[i], dmp_degree(g, v))

    if v > 0:
        v, i = v-1, i+1

        for c in g:
            _rec_degree_list(c, v, i, degs)

@cythonized("u")
def dmp_degree_list(f, u):
    """Returns a list of degrees of `f` in `K[X]`. """
    degs = [-1]*(u+1)
    _rec_degree_list(f, u, 0, degs)
    return tuple(degs)

@cythonized("i")
def dup_strip(f):
    """Remove leading zeros from `f` in `K[x]`. """
    if not f or f[0]:
        return f

    i = 0

    for cf in f:
        if cf:
            break
        else:
            i += 1

    return f[i:]

@cythonized("u,v,i")
def dmp_strip(f, u):
    """Remove leading zeros from `f` in `K[X]`. """
    if not u:
        return dup_strip(f)

    if dmp_zero_p(f, u):
        return f

    i, v = 0, u-1

    for c in f:
        if not dmp_zero_p(c, v):
            break
        else:
            i += 1

    if i == len(f):
        return dmp_zero(u)
    else:
        return f[i:]

@cythonized("i,j")
def _rec_validate(f, g, i, K):
    """XXX"""
    if type(g) is not list:
        if K is not None and not K.of_type(g):
            raise TypeError("%s in %s in not of type %s" % (g, f, K.dtype))

        return set([i-1])
    elif not g:
        return set([i])
    else:
        j, levels = i+1, set([])

        for c in g:
            levels |= _rec_validate(f, c, i+1, K)

        return levels

@cythonized("v,w")
def _rec_strip(g, v):
    """XXX"""
    if not v:
        return dup_strip(g)

    w = v-1

    return dmp_strip([ _rec_strip(c, w) for c in g ], v)

@cythonized("u")
def dmp_validate(f, K=None):
    """Returns number of levels in `f` and recursively strips it. """
    levels = _rec_validate(f, f, 0, K)

    u = levels.pop()

    if not levels:
        return _rec_strip(f, u), u
    else:
        raise ValueError("invalid data structure for a multivariate polynomial")

def dup_reverse(f):
    """Compute x**n * f(1/x), i.e.: reverse `f` in `K[x]`. """
    return dup_strip(list(reversed(f)))

def dup_copy(f):
    """Create a new copy of a polynomial `f` in `K[x]`. """
    return list(f)

@cythonized("u,v")
def dmp_copy(f, u):
    """Create a new copy of a polynomial `f` in `K[X]`. """
    if not u:
        return list(f)

    v = u-1

    return [ dmp_copy(c, v) for c in f ]

def dup_normal(f, K):
    """Normalize univariate polynomial in the given domain. """
    return dup_strip([ K.normal(c) for c in f ])

@cythonized("u,v")
def dmp_normal(f, u, K):
    """Normalize multivariate polynomial in the given domain. """
    if not u:
        return dup_normal(f, K)

    v = u-1

    return dmp_strip([ dmp_normal(c, v, K) for c in f ], u)

def dup_convert(f, K0, K1):
    """Convert ground domain of `f` from `K0` to `K1`. """
    if K0 is not None and K0 == K1:
        return f
    else:
        return dup_strip([ K1.convert(c, K0) for c in f ])

@cythonized("u,v")
def dmp_convert(f, u, K0, K1):
    """Convert ground domain of `f` from `K0` to `K1`. """
    if not u:
        return dup_convert(f, K0, K1)
    if K0 is not None and K0 == K1:
        return f

    v = u-1

    return dmp_strip([ dmp_convert(c, v, K0, K1) for c in f ], u)

def dup_from_sympy(f, K):
    """Convert ground domain of `f` from SymPy to `K`. """
    return dup_strip([ K.from_sympy(c) for c in f ])

@cythonized("u,v")
def dmp_from_sympy(f, u, K):
    """Convert ground domain of `f` from SymPy to `K`. """
    if not u:
        return dup_from_sympy(f, K)

    v = u-1

    return dmp_strip([ dmp_from_sympy(c, v, K) for c in f ], u)

@cythonized("n")
def dup_nth(f, n, K):
    """Returns n-th coefficient of `f` in `K[x]`. """
    if n < 0:
        raise IndexError("`n` must be non-negative, got %i" % n)
    elif n >= len(f):
        return K.zero
    else:
        return f[dup_degree(f)-n]

@cythonized("n,u")
def dmp_nth(f, n, u, K):
    """Returns n-th coefficient of `f` in `K[x]`. """
    if n < 0:
        raise IndexError("`n` must be non-negative, got %i" % n)
    elif n >= len(f):
        return dmp_zero(u-1)
    else:
        return f[dmp_degree(f, u)-n]

@cythonized("n,u,v")
def dmp_ground_nth(f, N, u, K):
    """Returns ground n-th coefficient of `f` in `K[x]`. """
    v = u

    for n in N:
        if n < 0:
            raise IndexError("`n` must be non-negative, got %i" % n)
        elif n >= len(f):
            return K.zero
        else:
            f, v = f[dmp_degree(f, v)-n], v-1

    return f

@cythonized("u")
def dmp_zero_p(f, u):
    """Returns True if `f` is zero in `K[X]`. """
    while u:
        if len(f) != 1:
            return False
        f = f[0]
        u -= 1
    return not f

@cythonized("u")
def dmp_zero(u):
    """Returns a multivariate zero. """
    r = []
    for i in xrange(u):
        r = [r]
    return r

@cythonized("u")
def dmp_one_p(f, u, K):
    """Returns True if `f` is one in `K[X]`. """
    return dmp_ground_p(f, K.one, u)

@cythonized("u")
def dmp_one(u, K):
    """Returns a multivariate one over `K`. """
    return dmp_ground(K.one, u)

@cythonized("u")
def dmp_ground_p(f, c, u):
    """Returns True if `f` is constant in `K[X]`. """
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

@cythonized("u")
def dmp_ground(c, u):
    """Returns a multivariate constant. """
    if not c:
        return dmp_zero(u)

    for i in xrange(u + 1):
        c = [c]
    return c

@cythonized("n,u")
def dmp_zeros(n, u, K):
    """Returns a list of multivariate zeros. """
    if not n:
        return []

    if u < 0:
        return [K.zero]*n
    else:
        return [ dmp_zero(u) for i in xrange(n) ]

@cythonized("n,u")
def dmp_grounds(c, n, u):
    """Returns a list of multivariate constants. """
    if not n:
        return []

    if u < 0:
        return [c]*n
    else:
        return [ dmp_ground(c, u) for i in xrange(n) ]

@cythonized("u")
def dmp_negative_p(f, u, K):
    """Returns `True` if `LC(f)` is negative. """
    return K.is_negative(dmp_ground_LC(f, u, K))

@cythonized("u")
def dmp_positive_p(f, u, K):
    """Returns `True` if `LC(f)` is positive. """
    return K.is_positive(dmp_ground_LC(f, u, K))

@cythonized("n,k")
def dup_from_dict(f, K):
    """Create `K[x]` polynomial from a dict. """
    if not f:
        return []

    n, h = max(f.iterkeys()), []

    if type(n) is int:
        for k in xrange(n, -1, -1):
            h.append(f.get(k, K.zero))
    else:
        (n,) = n

        for k in xrange(n, -1, -1):
            h.append(f.get((k,), K.zero))

    return dup_strip(h)

@cythonized("n,k")
def dup_from_raw_dict(f, K):
    """Create `K[x]` polynomial from a raw dict. """
    if not f:
        return []

    n, h = max(f.iterkeys()), []

    for k in xrange(n, -1, -1):
        h.append(f.get(k, K.zero))

    return dup_strip(h)

@cythonized("u,v,n,k")
def dmp_from_dict(f, u, K):
    """Create `K[X]` polynomial from a dict. """
    if not u:
        return dup_from_dict(f, K)
    if not f:
        return dmp_zero(u)

    coeffs = {}

    for monom, coeff in f.iteritems():
        head, tail = monom[0], monom[1:]

        if coeffs.has_key(head):
            coeffs[head][tail] = coeff
        else:
            coeffs[head] = { tail : coeff }

    n, v, h = max(coeffs.iterkeys()), u-1, []

    for k in xrange(n, -1, -1):
        coeff = coeffs.get(k)

        if coeff is not None:
            h.append(dmp_from_dict(coeff, v, K))
        else:
            h.append(dmp_zero(v))

    return dmp_strip(h, u)

@cythonized("n,k")
def dup_to_dict(f):
    """Convert `K[x]` polynomial to a dict. """
    n, result = dup_degree(f), {}

    for k in xrange(0, n+1):
        if f[n-k]:
            result[(k,)] = f[n-k]

    return result

@cythonized("n,k")
def dup_to_raw_dict(f):
    """Convert `K[x]` polynomial to a raw dict. """
    n, result = dup_degree(f), {}

    for k in xrange(0, n+1):
        if f[n-k]:
            result[k] = f[n-k]

    return result

@cythonized("u,v,n,k")
def dmp_to_dict(f, u):
    """Convert `K[X]` polynomial to a dict. """
    if not u:
        return dup_to_dict(f)

    n, v, result = dmp_degree(f, u), u-1, {}

    for k in xrange(0, n+1):
        h = dmp_to_dict(f[n-k], v)

        for exp, coeff in h.iteritems():
            result[(k,)+exp] = coeff

    return result

@cythonized("u,i,j")
def dmp_swap(f, i, j, u, K):
    """Transform `K[..x_i..x_j..]` to `K[..x_j..x_i..]`. """
    if i < 0 or j < 0 or i > u or j > u:
        raise IndexError("0 <= i < j <= %s expected" % u)
    elif i == j:
        return f

    F, H = dmp_to_dict(f, u), {}

    for exp, coeff in F.iteritems():
        H[exp[:i]    + (exp[j],) +
          exp[i+1:j] +
          (exp[i],)  + exp[j+1:]] = coeff

    return dmp_from_dict(H, u, K)

@cythonized("u")
def dmp_permute(f, P, u, K):
    """Returns a polynomial in `K[x_{P(1)},..,x_{P(n)}]`. """
    F, H = dmp_to_dict(f, u), {}

    for exp, coeff in F.iteritems():
        new_exp = [0]*len(exp)

        for e, p in zip(exp, P):
            new_exp[p] = e

        H[tuple(new_exp)] = coeff

    return dmp_from_dict(H, u, K)

@cythonized("l")
def dmp_nest(f, l, K):
    """Returns multivariate value nested `l`-levels. """
    if not isinstance(f, list):
        return dmp_ground(f, l)
    for i in xrange(l):
        f = [f]
    return f

@cythonized("l,k,u,v")
def dmp_raise(f, l, u, K):
    """Returns multivariate polynomial raised `l`-levels. """
    if not l:
        return f

    if not u:
        if not f:
            return dmp_zero(l)

        k = l-1

        return [ dmp_ground(c, k) for c in f ]

    v = u-1

    return [ dmp_raise(c, l, v, K) for c in f ]

@cythonized("g,i")
def dup_deflate(f, K):
    """Map `x**m` to `y` in a polynomial in `K[x]`. """
    if dup_degree(f) <= 0:
        return 1, f

    g = 0

    for i in xrange(len(f)):
        if not f[-i-1]:
            continue

        g = igcd(g, i)

        if g == 1:
            return 1, f

    return g, f[::g]

@cythonized("u,i,m,a,b")
def dmp_deflate(f, u, K):
    """Map `x_i**m_i` to `y_i` in a polynomial in `K[X]`. """
    if dmp_zero_p(f, u):
        return (1,)*(u+1), f

    F = dmp_to_dict(f, u)
    B = [0]*(u+1)

    for M in F.iterkeys():
        for i, m in enumerate(M):
            B[i] = igcd(B[i], m)

    for i, b in enumerate(B):
        if not b:
            B[i] = 1

    B = tuple(B)

    if all([ b == 1 for b in B ]):
        return B, f

    H = {}

    for A, coeff in F.iteritems():
        N = [ a // b for a, b in zip(A, B) ]
        H[tuple(N)] = coeff

    return B, dmp_from_dict(H, u, K)

@cythonized("G,g,i")
def dup_multi_deflate(polys, K):
    """Map `x**m` to `y` in a set of polynomials in `K[x]`. """
    G = 0

    for p in polys:
        if dup_degree(p) <= 0:
            return 1, polys

        g = 0

        for i in xrange(len(p)):
            if not p[-i-1]:
                continue

            g = igcd(g, i)

            if g == 1:
                return 1, polys

        G = igcd(G, g)

    return G, tuple([ p[::G] for p in polys ])

@cythonized("u,G,g,m,a,b")
def dmp_multi_deflate(polys, u, K):
    """Map `x_i**m_i` to `y_i` in a set of polynomials in `K[X]`. """
    if not u:
        M, H = dup_multi_deflate(polys, K)
        return (M,), H

    F, B = [], [0]*(u+1)

    for p in polys:
        f = dmp_to_dict(p, u)

        if not dmp_zero_p(p, u):
            for M in f.iterkeys():
                for i, m in enumerate(M):
                    B[i] = igcd(B[i], m)

        F.append(f)

    for i, b in enumerate(B):
        if not b:
            B[i] = 1

    B = tuple(B)

    if all([ b == 1 for b in B ]):
        return B, polys

    H = []

    for f in F:
        h = {}

        for A, coeff in f.iteritems():
            N = [ a // b for a, b in zip(A, B) ]
            h[tuple(N)] = coeff

        H.append(dmp_from_dict(h, u, K))

    return B, tuple(H)

@cythonized("m")
def dup_inflate(f, m, K):
    """Map `y` to `x**m` in a polynomial in `K[x]`. """
    if m <= 0:
        raise IndexError("`m` must be positive, got %s" % m)
    if m == 1 or not f:
        return f

    result = [f[0]]

    for coeff in f[1:]:
        result.extend([K.zero]*(m-1))
        result.append(coeff)

    return result

@cythonized("u,v,i,j")
def _rec_inflate(g, M, v, i, K):
    """XXX"""
    if not v:
        return dup_inflate(g, M[i], K)
    if M[i] <= 0:
        raise IndexError("all `M[i]` must be positive, got %s" % M[i])

    w, j = v-1, i+1

    g = [ _rec_inflate(c, M, w, j, K) for c in g ]

    result = [g[0]]

    for coeff in g[1:]:
        for _ in xrange(1, M[i]):
            result.append(dmp_zero(w))

        result.append(coeff)

    return result

@cythonized("u,m")
def dmp_inflate(f, M, u, K):
    """Map `y_i` to `x_i**k_i` in a polynomial in `K[X]`. """
    if not u:
        return dup_inflate(f, M[0], K)

    if all([ m == 1 for m in M ]):
        return f
    else:
        return _rec_inflate(f, M, u, 0, K)

@cythonized("u,j")
def dmp_exclude(f, u, K):
    """Exclude useless levels from `f`. """
    if not u or dmp_ground_p(f, None, u):
        return [], f, u

    J, F = [], dmp_to_dict(f, u)

    for j in xrange(0, u+1):
        for monom in F.iterkeys():
            if monom[j]:
                break
        else:
            J.append(j)

    if not J:
        return [], f, u

    f = {}

    for monom, coeff in F.iteritems():
        monom = list(monom)

        for j in reversed(J):
            del monom[j]

        f[tuple(monom)] = coeff

    u -= len(J)

    return J, dmp_from_dict(f, u, K), u

@cythonized("u,j")
def dmp_include(f, J, u, K):
    """Include useless levels in `f`. """
    if not J:
        return f

    F, f = dmp_to_dict(f, u), {}

    for monom, coeff in F.iteritems():
        monom = list(monom)

        for j in J:
            monom.insert(j, 0)

        f[tuple(monom)] = coeff

    u += len(J)

    return dmp_from_dict(f, u, K)

@cythonized("u,v,w")
def dmp_inject(f, u, K, front=False):
    """Convert `f` from `K[X][Y]` to `K[X,Y]`. """
    f, h = dmp_to_dict(f, u), {}

    v = len(K.gens) - 1

    for f_monom, g in f.iteritems():
        g = g.to_dict()

        for g_monom, c in g.iteritems():
            if front:
                h[g_monom + f_monom] = c
            else:
                h[f_monom + g_monom] = c


    w = u + v + 1

    return dmp_from_dict(h, w, K.dom), w

@cythonized("u,v")
def dmp_eject(f, u, K, front=False):
    """Convert `f` from `K[X,Y]` to `K[X][Y]`. """
    f, h = dmp_to_dict(f, u), {}

    v = u - len(K.gens) + 1

    for monom, c in f.iteritems():
        if front:
            g_monom, f_monom = monom[:v], monom[v:]
        else:
            f_monom, g_monom = monom[:v], monom[v:]

        if h.has_key(f_monom):
            h[f_monom][g_monom] = c
        else:
            h[f_monom] = {g_monom: c}

    for monom, c in h.iteritems():
        h[monom] = K(c)

    return dmp_from_dict(h, v-1, K)

@cythonized("i")
def dup_terms_gcd(f, K):
    """Remove GCD of terms from `f` in `K[x]`. """
    if not f or f[-1]:
        return 0, f

    i = 0

    for c in reversed(f):
        if not c:
            i += 1
        else:
            break

    return i, f[:-i]

@cythonized("u,g")
def dmp_terms_gcd(f, u, K):
    """Remove GCD of terms from `f` in `K[X]`. """
    if dmp_ground_TC(f, u, K) or dmp_zero_p(f, u):
        return (0,)*(u+1), f

    F = dmp_to_dict(f, u)
    G = monomial_min(*F.keys())

    if all([ g == 0 for g in G ]):
        return G, f

    f = {}

    for monom, coeff in F.iteritems():
        f[monomial_div(monom, G)] = coeff

    return G, dmp_from_dict(f, u, K)

@cythonized("v,w,d,i")
def _rec_list_terms(g, v, monom):
    """XXX"""
    d, terms = dmp_degree(g, v), []

    if not v:
        for i, c in enumerate(g):
            if not c:
                continue

            terms.append((monom + (d-i,), c))
    else:
        w = v-1

        for i, c in enumerate(g):
            terms.extend(_rec_list_terms(c, v-1, monom + (d-i,)))

    return terms

@cythonized("u")
def dmp_list_terms(f, u, K, order=None):
    """List all non-zero terms from `f` in the given order order. """
    terms = _rec_list_terms(f, u, ())

    if not terms:
        return [((0,)*(u+1), K.zero)]

    if order is None:
        return terms
    else:
        return sdp_sort(terms, monomial_key(order))

@cythonized("n,m")
def dup_apply_pairs(f, g, h, args, K):
    """Apply `h` to pairs of coefficients of `f` and `g`. """
    n, m = len(f), len(g)

    if n != m:
        if n > m:
            g = [K.zero]*(n-m) + g
        else:
            f = [K.zero]*(m-n) + f

    result = []

    for a, b in zip(f, g):
        result.append(h(a, b, *args))

    return dup_strip(result)

@cythonized("u,v,n,m")
def dmp_apply_pairs(f, g, h, args, u, K):
    """Apply `h` to pairs of coefficients of `f` and `g`. """
    if not u:
        return dup_apply_pairs(f, g, h, args, K)

    n, m, v = len(f), len(g), u-1

    if n != m:
        if n > m:
            g = dmp_zeros(n-m, v, K) + g
        else:
            f = dmp_zeros(m-n, v, K) + f

    result = []

    for a, b in zip(f, g):
        result.append(dmp_apply_pairs(a, b, h, args, v, K))

    return dmp_strip(result, u)

def dup_random(n, a, b, K):
    """Return a polynomial of degree ``n`` with coefficients in ``[a, b]``. """
    f = [ K.convert(random.randint(a, b)) for _ in xrange(0, n+1) ]

    while not f[0]:
        f[0] = K.convert(random.randint(a, b))

    return f


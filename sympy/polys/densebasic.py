"""Basic tools for dense recursive polynomials in `K[x]` or `K[X]`. """

from sympy.utilities import any, all
from sympy.core import igcd, ilcm

from sympy.polys.monomialtools import (
    monomial_min, monomial_div,
)

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

def dmp_ground_LC(f, u, K):
    """Returns ground leading coefficient. """
    if not u:
        return dup_LC(f, K)
    else:
        return dmp_ground_LC(dmp_LC(f, K), u-1, K)

def dmp_ground_TC(f, u, K):
    """Returns ground trailing coefficient. """
    if not u:
        return dup_TC(f, K)
    else:
        return dmp_ground_TC(dmp_TC(f, K), u-1, K)

def dup_degree(f):
    """Returns leading degree of `f` in `K[x]`. """
    return len(f) - 1

def dmp_degree(f, u):
    """Returns leading degree of `f` in `x_0` in `K[X]`. """
    if dmp_zero_p(f, u):
        return -1
    else:
        return len(f) - 1

def dmp_degree_in(f, j, u):
    """Returns leading degree of `f` in `x_j` in `K[X]`. """
    if not j:
        return dmp_degree(f, u)
    if j < 0 or j > u:
        raise IndexError("-%s <= j < %s expected, got %s" % (u, u, j))

    def rec_degree_in(g, v, i):
        if i == j:
            return dmp_degree(g, v)
        else:
            return max([ rec_degree_in(c, v-1, i+1) for c in g ])

    return rec_degree_in(f, u, 0)

def dmp_degree_list(f, u):
    """Returns a list of degrees of `f` in `K[X]`. """
    degs = [-1]*(u+1)

    def rec_degree_list(g, v, i):
        degs[i] = max(degs[i], dmp_degree(g, v))

        if v > 0:
            for c in g:
                rec_degree_list(c, v-1, i+1)

    rec_degree_list(f, u, 0)
    return tuple(degs)

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

def dmp_validate(f, K=None):
    """Returns number of levels in `f` and recursively strips it. """
    levels = set([])

    def rec_validate(g, i):
        if type(g) is not list:
            if K is not None and not K.of_type(g):
                raise TypeError("%s in %s in not of type %s" % (g, f, K.dtype))

            levels.add(i-1)
        elif not g:
            levels.add(i)
        else:
            for c in g:
                rec_validate(c, i+1)

    def rec_strip(g, v):
        if not v:
            return dup_strip(g)
        else:
            return dmp_strip([ rec_strip(c, v-1) for c in g ], v)

    rec_validate(f, 0)
    u = levels.pop()

    if not levels:
        return rec_strip(f, u), u
    else:
        raise ValueError("invalid data structure for a multivariate polynomial")

def dup_copy(f):
    """Create a new copy of a polynomial `f` in `K[x]`. """
    return list(f)

def dmp_copy(f, u):
    """Create a new copy of a polynomial `f` in `K[X]`. """
    if not u:
        return list(f)
    else:
        return [ dmp_copy(c, u-1) for c in f ]

def dup_normal(f, K):
    """Normalize univariate polynomial in the given domain. """
    return dup_strip([ K(c) for c in f ])

def dmp_normal(f, u, K):
    """Normalize multivariate polynomial in the given domain. """
    if not u:
        return dup_normal(f, K)
    else:
        return dmp_strip([ dmp_normal(c, u-1, K) for c in f ], u)

def dup_convert(f, K0, K1):
    """Convert ground domain of `f` from `K0` to `K1`. """
    if K0 == K1:
        return f
    else:
        return dup_strip([ K1.convert(c, K0) for c in f ])

def dmp_convert(f, u, K0, K1):
    """Convert ground domain of `f` from `K0` to `K1`. """
    if K0 == K1:
        return f
    elif not u:
        return dup_convert(f, K0, K1)
    else:
        return dmp_strip([ dmp_convert(c, u-1, K0, K1) for c in f ], u)

def dup_nth(f, n, K):
    """Returns n-th coefficient of `f` in `K[x]`. """
    if n < 0:
        raise IndexError("`n` must be non-negative, got %i" % n)
    elif n >= len(f):
        return K.zero
    else:
        return f[dup_degree(f)-n]

def dmp_nth(f, n, u, K):
    """Returns n-th coefficient of `f` in `K[x]`. """
    if n < 0:
        raise IndexError("`n` must be non-negative, got %i" % n)
    elif n >= len(f):
        return dmp_zero(u-1)
    else:
        return f[dmp_degree(f, u)-n]

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

def dmp_zero_p(f, u):
    """Returns True if `f` is zero in `K[X]`. """
    if not u:
        return not f
    else:
        if len(f) == 1:
            return dmp_zero_p(f[0], u-1)
        else:
            return False

def dmp_zero(u):
    """Returns a multivariate zero. """
    if not u:
        return []
    else:
        return [dmp_zero(u-1)]

def dmp_one_p(f, u, K):
    """Returns True if `f` is one in `K[X]`. """
    return dmp_ground_p(f, K.one, u)

def dmp_one(u, K):
    """Returns a multivariate one over `K`. """
    return dmp_ground(K.one, u)

def dmp_ground_p(f, c, u):
    """Returns True if `f` is constant in `K[X]`. """
    if c is not None:
        if not c:
            return dmp_zero_p(f, u)
        if not u:
            return f == [c]
    else:
        if not u:
            return len(f) <= 1

    if len(f) == 1:
        return dmp_ground_p(f[0], c, u-1)
    else:
        return False

def dmp_ground(c, u):
    """Returns a multivariate constant. """
    if not c:
        return dmp_zero(u)
    else:
        if u < 0:
            return c
        if not u:
            return [c]
        else:
            return [dmp_ground(c, u-1)]

def dmp_zeros(n, u, K):
    """Returns a list of multivariate zeros. """
    if not n:
        return []

    if u < 0:
        return [K.zero]*n
    else:
        return [ dmp_zero(u) for i in xrange(n) ]

def dmp_grounds(c, n, u):
    """Returns a list of multivariate constants. """
    if not n:
        return []

    if u < 0:
        return [c]*n
    else:
        return [ dmp_ground(c, u) for i in xrange(n) ]

def dmp_negative_p(f, u, K):
    """Returns `True` if `LC(f)` is negative. """
    return K.is_negative(dmp_ground_LC(f, u, K))

def dmp_positive_p(f, u, K):
    """Returns `True` if `LC(f)` is positive. """
    return K.is_positive(dmp_ground_LC(f, u, K))

def dup_from_dict(f, K):
    """Create `K[x]` polynomial from a dict. """
    if not f:
        return []

    n, h = max(k for (k,) in f.iterkeys()), []

    for k in xrange(n, -1, -1):
        h.append(f.get((k,), K.zero))

    return dup_strip(h)

def dup_from_raw_dict(f, K):
    """Create `K[x]` polynomial from a raw dict. """
    if not f:
        return []

    n, h = max(f.iterkeys()), []

    for k in xrange(n, -1, -1):
        h.append(f.get(k, K.zero))

    return dup_strip(h)

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

def dup_to_dict(f):
    """Convert `K[x]` polynomial to a dict. """
    n, result = dup_degree(f), {}

    for i in xrange(0, n+1):
        if f[n-i]:
            result[(i,)] = f[n-i]

    return result

def dup_to_raw_dict(f):
    """Convert `K[x]` polynomial to a raw dict. """
    n, result = dup_degree(f), {}

    for i in xrange(0, n+1):
        if f[n-i]:
            result[i] = f[n-i]

    return result

def dmp_to_dict(f, u):
    """Convert `K[X]` polynomial to a dict. """
    if not u:
        return dup_to_dict(f)

    n, v, result = dmp_degree(f, u), u-1, {}

    for i in xrange(0, n+1):
        h = dmp_to_dict(f[n-i], v)

        for exp, coeff in h.iteritems():
            result[(i,)+exp] = coeff

    return result

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

def dmp_permute(f, P, u, K):
    """Returns a polynomial in `K[x_{P(1)},..,x_{P(n)}]`. """
    F, H = dmp_to_dict(f, u), {}

    for exp, coeff in F.iteritems():
        new_exp = [0]*len(exp)

        for e, p in zip(exp, P):
            new_exp[p] = e

        H[tuple(new_exp)] = coeff

    return dmp_from_dict(H, u, K)

def dmp_nest(f, l, K):
    """Returns multivariate value nested `l`-levels. """
    if type(f) is not list:
        return dmp_ground(f, l)
    else:
        if not l:
            return f
        else:
            return [dmp_nest(f, l-1, K)]

def dmp_raise(f, l, u, K):
    """Returns multivariate polynomial raised `l`-levels. """
    if not l:
        return f

    if not u:
        if not f:
            return dmp_zero(l)
        else:
            return [ dmp_ground(c, l-1) for c in f ]
    else:
        return [ dmp_raise(c, l, u-1, K) for c in f ]

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

def dmp_deflate(f, u, K):
    """Map `x_i**m_i` to `y_i` in a polynomial in `K[X]`. """
    if dmp_zero_p(f, u):
        return (1,)*(u+1), f

    F, H = dmp_to_dict(f, u), {}

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

    return M, dmp_from_dict(H, u, K)

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

def dmp_multi_deflate(polys, u, K):
    """Map `x_i**m_i` to `y_i` in a set of polynomials in `K[X]`. """
    def ilgcd(M):
        g = 0

        for m in M:
            g = igcd(g, m)

            if g == 1:
                break

        return g or 1

    if not u:
        M, H = dup_multi_deflate(polys, K)
        return (M,), H

    F, M, H = [], [], []

    for p in polys:
        f = dmp_to_dict(p, u)

        if dmp_zero_p(p, u):
            m = (0,)*(u+1)
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

        H.append(dmp_from_dict(h, u, K))

    return M, tuple(H)

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

def dmp_inflate(f, M, u, K):
    """Map `y_i` to `x_i**k_i` in a polynomial in `K[X]`. """
    if not u:
        return dup_inflate(f, M[0], K)

    def rec_inflate(g, v, i):
        if not v:
            return dup_inflate(g, M[i], K)

        if M[i] <= 0:
            raise IndexError("all `M[i]` must be positive, got %s" % M[i])

        g = [ rec_inflate(c, v-1, i+1) for c in g ]
        result, w = [g[0]], v-1

        for coeff in g[1:]:
            for _ in xrange(1, M[i]):
                result.append(dmp_zero(w))

            result.append(coeff)

        return result

    if all([ m == 1 for m in M ]):
        return f
    else:
        return rec_inflate(f, u, 0)

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

def dmp_inject(f, u, K, **args):
    """Convert `f` from `K[X][Y]` to `K[X,Y]`. """
    front = args.get('front', False)
    f, h = dmp_to_dict(f, u), {}

    v = len(K.gens) - 1

    for f_monom, g in f.iteritems():
        g = g.to_dict()

        for g_monom, c in g.iteritems():
            if not front:
                h[f_monom + g_monom] = c
            else:
                h[g_monom + f_monom] = c

    w = u + v + 1

    return dmp_from_dict(h, w, K.dom), w

def dmp_eject(f, u, K, **args):
    """Convert `f` from `K[X,Y]` to `K[X][Y]`. """
    front = args.get('front', False)
    f, h = dmp_to_dict(f, u), {}

    v = u - len(K.gens) + 1

    for monom, c in f.iteritems():
        if not front:
            f_monom, g_monom = monom[:v], monom[v:]
        else:
            g_monom, f_monom = monom[:v], monom[v:]

        if h.has_key(f_monom):
            h[f_monom][g_monom] = c
        else:
            h[f_monom] = {g_monom: c}

    for monom, c in h.iteritems():
        h[monom] = K(c)

    return dmp_from_dict(h, v-1, K)

def dup_terms_gcd(f, K):
    """Remove GCD of terms from `f` in `K[x]`. """
    if dup_TC(f, K) or not f:
        return 0, f

    F = dup_to_raw_dict(f)

    gcd = min(F.keys())

    if not gcd:
        return gcd, f

    f = {}

    for k, coeff in F.iteritems():
        f[k - gcd] = coeff

    return gcd, dup_from_raw_dict(f, K)

def dmp_terms_gcd(f, u, K):
    """Remove GCD of terms from `f` in `K[X]`. """
    if dmp_ground_TC(f, u, K) or dmp_zero_p(f, u):
        return (0,)*(u+1), f

    F = dmp_to_dict(f, u)

    gcd = monomial_min(*F.keys())

    if all(not n for n in gcd):
        return gcd, f

    f = {}

    for monom, coeff in F.iteritems():
        f[monomial_div(monom, gcd)] = coeff

    return gcd, dmp_from_dict(f, u, K)

def dmp_list_terms(f, u, K):
    """List all non-zero terms from `f` in lex order. """
    terms = []

    def rec_iter_terms(g, v, monom):
        d = dmp_degree(g, v)

        if not v:
            for i, c in enumerate(g):
                if not c:
                    continue

                terms.append((monom + (d-i,), c))
        else:
            for i, c in enumerate(g):
                rec_iter_terms(c, v-1, monom + (d-i,))

    rec_iter_terms(f, u, ())

    if not terms:
        return [((0,)*(u+1), K.zero)]
    else:
        return terms

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


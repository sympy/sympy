"""Different modular algorithm used in integer polynomial arithmetic."""

from sympy.polynomials.base import *

# Modular integer arithmetic:

def xgcd(a, b):
    """
    Returns g, x, y such that g = x*a + y*b = gcd(a,b).
    Input:
        a -- an integer
        b -- an integer
    Output:
        g -- an integer, the gcd of a and b
        x -- an integer
        y -- an integer
    Examples:
    >>> xgcd(2,3)
    (1, -1, 1)
    >>> xgcd(10, 12)
    (2, -1, 1)
    >>> g, x, y = xgcd(100, 2004)
    >>> print g, x, y
    4 -20 1
    >>> print x*100 + y*2004
    4
    """
    if a == 0 and b == 0: return (0, 0, 1)
    if a == 0: return (abs(b), 0, b/abs(b))
    if b == 0: return (abs(a), a/abs(a), 0)
    x_sign = 1; y_sign = 1
    if a < 0: a = -a; x_sign = -1
    if b < 0: b = -b; y_sign = -1
    x = 1; y = 0; r = 0; s = 1
    while b != 0:
        (c, q) = (a%b, a/b)
        (a, b, r, s, x, y) = (b, c, x-q*r, y-q*s, r, s)
    return (a, x*x_sign, y*y_sign)

def mod_div(a, b, p):
    g, x, y = xgcd(b, p)
    assert g == 1, "Zero division!"
    return (a * x) % p

# Simple polynomial arithmetic:

def poly_new(coeffs, p):
    result = {}
    for term in coeffs:
        c = int(term[0]) % p
        if c:
            result[int(term[1])] = c
    return result

def poly2coeffs(f, p):
    result = []
    keys = f.keys()
    keys.sort(reverse=True)
    for exponent in keys:
        if f[exponent] <= p/2:
            result.append((Integer(f[exponent]), Integer(exponent)))
        else:
            result.append((Integer(f[exponent] - p), Integer(exponent)))
    return tuple(result)

def poly_deg(f):
    if len(f) == 0:
        return -1
    return max(f.keys())

def poly_add(f, g, p):
    if len(f) < len(g):
        f, g = g, f
    result = f.copy()
    for exponent in g.keys():
        if result.has_key(exponent):
            coeff = (result[exponent] + g[exponent]) % p
            if coeff:
                result[exponent] = coeff
            else: # remove 0 coeffs
                del result[exponent]
        else:
            result[exponent] = g[exponent]
    return result

def poly_mul(f, g, p):
    if len(f) < len(g):
        f, g = g, f
    result = {}
    for exponent in g.keys():
        result = poly_add(result, poly_scale(f, g[exponent], exponent, p), p)
    return result

def poly_multi_mul(factors, p):
    if len(factors) == 1:
        return factors[0]
    elif len(factors) == 2:
        return poly_mul(factors[0], factors[1], p)
    else:
        return poly_mul(poly_multi_mul(factors[:len(factors)/2], p),
                        poly_multi_mul(factors[len(factors)/2:], p),
                        p)
    
def poly_scale(f, c, exp, p):
    if not (c % p):
        return {}
    result = {}
    for exponent in f.keys():
        result[exponent + exp] = (c * f[exponent]) % p
    return result

def poly_pow(f, n, p):
    """Repeated squaring."""
    assert isinstance(n, (int, long))
    if not n:
        return {0: 1}
    binary_n = []
    while n:
        if n % 2:
            binary_n.insert(0, 1)
            n = (n - 1) / 2
        else:
            binary_n.insert(0, 0)
            n /= 2
    result = f
    for k in binary_n[1:]:
        result = poly_mul(result, result, p)
        if k:
            result = poly_mul(result, f, p)
    return result

def poly_div(f, g, p):
    q, r = {}, f.copy()
    deg_g = poly_deg(g)
    deg_r = poly_deg(r)
    while deg_r >= deg_g:
        deg_quot = deg_r - deg_g
        c_quot = mod_div(r[deg_r], g[deg_g], p)
        q = poly_add(q, {deg_quot: c_quot}, p)
        r = poly_add(r, poly_scale(g, -c_quot, deg_quot, p), p)
        deg_r = poly_deg(r)
    return q, r

def poly_gcd(f, g, p):
    while poly_deg(g) >= 0:
        f, g = g, poly_div(f, g, p)[1]
    return f

def poly_monic(f, p):
    lc = f[poly_deg(f)]
    if lc == 1:
        return f.copy()
    factor = mod_div(1, lc, p)
    return poly_scale(f, factor, 0, p) 

def poly_rand(d_min, d_max, p):
    import random
    degree = random.randrange(d_min, d_max + 1)
    result = {degree: 1}
    for k in range(0, degree):
        coeff = random.randrange(p)
        if coeff:
            result[k] = coeff
    return result

def poly_diff(f, p):
    result = {}
    for exponent in f.keys():
        e = exponent % p
        if e:
            result[exponent-1] = (e*f[exponent]) % p
    return result

# Modular polynomial arithmetic

def mod_poly_truncate(f, n):
    result = {}
    for exponent in f.keys():
        if exponent < n:
            result[exponent] = f[exponent]
    return result

def mod_poly_pow(f, n, p, m):
    """Repeated squaring."""
    """Repeated squaring."""
    assert isinstance(n, (int, long))
    if not n:
        return {0: 1}
    binary_n = []
    while n:
        if n % 2:
            binary_n.insert(0, 1)
            n = (n - 1) / 2
        else:
            binary_n.insert(0, 0)
            n /= 2
    result = poly_div(f, m, p)[1]
    for k in binary_n[1:]:
        result = poly_mul(result, result, p)
        result = poly_div(result, m, p)[1]
        if k:
            result = poly_mul(result, f, p)
            result = poly_div(result, m, p)[1]
    return result


# Factorization

def distinct_deg_factor(f, p):
    result = []
    f = f.copy()
    s = poly_deg(f)
    h = {1: 1}
    while not (len(f) == 1 and f.has_key(0)):
        h = mod_poly_pow(h, p, p, f) # h = h**p mod (f, p)
        g = poly_gcd(poly_add(h, {1: -1}, p), f, p)
        g = poly_monic(g, p)
        f, r = poly_div(f, g, p)
        assert r == {}
        result.append(g)
        # Early abort:
        if poly_deg(f) < 2*(len(g) + 1):
            result.append(f)
            break
    return result

def equal_deg_split(f, d, p):
    a = poly_rand(1, poly_deg(f) - 1, p)
    g = poly_monic(poly_gcd(a, f, p), p)
    if g != {0: 1}: # Found a non-trivial factor
        return g
    b = mod_poly_pow(a, (p**d - 1)/2, p, f)
    g = poly_monic(poly_gcd(poly_add(b, {0: -1}, p), f, p), p)
    if g not in [{0: 1}, f]:
        return g
    return None

def equal_deg_factor(f, d, p):
    if d == poly_deg(f):
        return [f]
    g = None
    while g is None:
        g = equal_deg_split(f, d, p)
    q, r = poly_div(f, g, p)
    assert r == {}
    return equal_deg_factor(g, d, p) + equal_deg_factor(q, d, p)

def factor(f, p):
    result = []
    f = poly_monic(f, p)
    dd = distinct_deg_factor(f, p)
    for d, ff in enumerate(dd):
        if ff == {0: 1}:
            continue
        d += 1
        result += equal_deg_factor(ff, d, p)
    return result


"""Univariate polynomials with coefficients in Galois fields. """

from sympy.core.numbers import igcdex

def zp_inv(a, p):
    """Compute multiplicative inverse over GF(p). """
    s, t, g = igcdex(a, p)

    if g == 1:
        return s % p
    else:
        raise ZeroDivisionError("modular division")

def zp_pow(a, n, p):
    """Compute a**n mod p using repeated squaring. """
    if n < 0:
        a, n = zp_inv(a), -n

    b = 1

    while n:
        if n & 1:
            b *= a
            n -= 1

        a *= a
        a %= p

        n >>= 1

    return b % p

def zp_int(a, p):
    """Coerce a mod p to an integer in [-p/2, p/2] range. """
    if a <= p // 2:
        return a
    else:
        return a - p

from random import uniform
from math import ceil, sqrt, log

def gf_degree(f):
    """Returns leading degree of f. """
    return len(f)-1

def gf_LC(f):
    """Returns leading coefficient of f. """
    if not f:
        return 0
    else:
        return f[0]

def gf_TC(f):
    """Returns trailing coefficient of f. """
    if not f:
        return 0
    else:
        return f[-1]

def gf_strip(f):
    """Remove leading zeros from f. """
    if not f or f[0]:
        return f

    k = 0

    for coeff in f:
        if coeff:
            break
        else:
            k += 1

    return f[k:]

def gf_reverse(f):
    """Reverse order of coefficients in f. """
    return gf_strip(list(reversed(f)))

def gf_normal(f, p):
    """Reduce all coefficients modulo p. """
    return gf_strip([ a % p for a in f ])

def gf_from_dict(f, p):
    """Create GF(p)[x] polynomial from a dict. """
    n, h = max(f.iterkeys()), []

    for k in xrange(n, -1, -1):
        h.append(f.get(k, 0) % p)

    return gf_normal(h, p)

def gf_to_dict(f, p):
    """Convert GF(p)[x] polynomial to a dict. """
    n, result = gf_degree(f), {}

    for i in xrange(0, n+1):
        a = zp_int(f[n-i], p)
        if a: result[i] = a

    return result

def gf_from_int_poly(f, p):
    """Create GF(p)[x] polynomial from Z[x]. """
    return gf_normal(f, p)

def gf_to_int_poly(f, p):
    """Convert GF(p)[x] polynomial to Z[x]. """
    return [ zp_int(c, p) for c in f ]

def gf_neg(f, p):
    """Negate a polynomial over GF(p)[x]. """
    return [ -coeff % p for coeff in f ]

def gf_add_const(f, a, p):
    """Returns f + a where f in GF(p)[x] and a in GF(p). """
    if not f:
        a = a % p
    else:
        a = (f[-1] + a) % p

        if len(f) > 1:
            return f[:-1] + [a]

    if not a:
        return []
    else:
        return [a]

def gf_sub_const(f, a, p):
    """Returns f - a where f in GF(p)[x] and a in GF(p). """
    if not f:
        a = -a % p
    else:
        a = (f[-1] - a) % p

        if len(f) > 1:
            return f[:-1] + [a]

    if not a:
        return []
    else:
        return [a]

def gf_mul_const(f, a, p):
    """Returns f * a where f in GF(p)[x] and a in GF(p). """
    if not a:
        return []
    else:
        return [ (a*b) % p for b in f ]

def gf_div_const(f, a, p):
    """Returns f / a where f in GF(p)[x] and a in GF(p). """
    return gf_mul_const(f, zp_inv(a, p), p)

def gf_add(f, g, p):
    """Add polynomials over GF(p)[x]. """
    if not f:
        return g
    if not g:
        return f

    df = gf_degree(f)
    dg = gf_degree(g)

    if df == dg:
        return gf_strip([ (a + b) % p for a, b in zip(f, g) ])
    else:
        k = abs(df - dg)

        if df > dg:
            h, f = f[:k], f[k:]
        else:
            h, g = g[:k], g[k:]

        return h + [ (a + b) % p for a, b in zip(f, g) ]

def gf_sub(f, g, p):
    """Subtract polynomials over GF(p)[x]. """
    if not g:
        return f
    if not f:
        return gf_neg(g, p)

    df = gf_degree(f)
    dg = gf_degree(g)

    if df == dg:
        return gf_strip([ (a - b) % p for a, b in zip(f, g) ])
    else:
        k = abs(df - dg)

        if df > dg:
            h, f = f[:k], f[k:]
        else:
            h, g = gf_neg(g[:k], p), g[k:]

        return h + [ (a - b) % p for a, b in zip(f, g) ]

def gf_mul(f, g, p):
    """Multiply polynomials over GF(p)[x]. """
    df = gf_degree(f)
    dg = gf_degree(g)

    dh = df + dg
    h = [0]*(dh+1)

    for i in xrange(0, dh+1):
        coeff = 0

        for j in xrange(max(0, i-dg), min(i, df)+1):
            coeff += f[j]*g[i-j]

        h[i] = coeff % p

    return gf_strip(h)

def gf_sqr(f, p):
    """Square polynomials over GF(p)[x]. """
    df = gf_degree(f)

    dh = 2*df
    h = [0]*(dh+1)

    for i in xrange(0, dh+1):
        coeff = 0

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

        h[i] = coeff % p

    return gf_strip(h)

def gf_add_mul(f, g, h, p):
    """Returns f + g*h where f, g, h in GF(p)[x]. """
    return gf_add(f, gf_mul(g, h, p), p)

def gf_sub_mul(f, g, h, p):
    """Returns f - g*h where f, g, h in GF(p)[x]. """
    return gf_sub(f, gf_mul(g, h, p), p)

def gf_expand(F, p):
    """Expand results of factor() over GF(p)[x]. """
    if type(F) is tuple:
        LC, F = F
    else:
        LC = 1

    g = [LC]

    for f, k in F:
        f = gf_pow(f, k, p)
        g = gf_mul(g, f, p)

    return g

def gf_div(f, g, p):
    """Division with remainder over GF(p)[x].

       Given univariate polynomials f and g over a finite field with p
       elements,  returns polynomials q and r (quotient and remainder)
       such that f = q*g + r.

       Consider polynomials x**3 + x + 1 and x**2 + x over GF(2):

       >>> from sympy.polys.galoispolys import gf_div, gf_add_mul
       >>> gf_div([1, 0, 1, 1], [1, 1, 0], 2)
       ([1, 1], [1])

       As result we obtained quotient x + 1 and remainder 1, thus:

       >>> gf_add_mul([1], [1, 1], [1, 1, 0], 2)
       [1, 0, 1, 1]

       For more details on the implemented algorithms refer to:

       [1] Michael Monagan, In-place Arithmetic for Polynomials over
           Z_n, Proceedings of DISCO '92, Springer-Verlag LNCS, 721,
           1993, pp. 22-34

       [2] J. von zur Gathen, J. Gerhard, Modern Computer Algebra,
           First Edition, Cambridge University Press, 1999, pp. 247

    """
    df = gf_degree(f)
    dg = gf_degree(g)

    if not g:
        raise ZeroDivisionError("polynomial division")
    elif df < dg:
        return [], f

    inv = zp_inv(g[0], p)

    h, dq, dr = f[:], df-dg, dg-1

    for i in xrange(0, df+1):
        coeff = h[i]

        for j in xrange(max(0, dg-i), min(df-i, dr)+1):
            coeff -= h[i+j-dg] * g[dg-j]

        if i <= dq:
            coeff *= inv

        h[i] = coeff % p

    return h[:dq+1], gf_strip(h[dq+1:])

def gf_quo(f, g, p):
    """Computes polynomial quotient over GF(p)[x]. """
    df = gf_degree(f)
    dg = gf_degree(g)

    if not g:
        raise ZeroDivisionError("polynomial division")
    elif df < dg:
        return []

    inv = zp_inv(g[0], p)

    h, dq, dr = f[:], df-dg, dg-1

    for i in xrange(0, dq+1):
        coeff = h[i]

        for j in xrange(max(0, dg-i), min(df-i, dr)+1):
            coeff -= h[i+j-dg] * g[dg-j]

        h[i] = (coeff * inv) % p

    return h[:dq+1]

def gf_rem(f, g, p):
    """Returns polynomial remainder over GF(p)[x]. """
    return gf_div(f, g, p)[1]

def gf_lshift(f, n):
    """Efficiently multiply f by x**n. """
    if not f:
        return f
    else:
        return f + [0]*n

def gf_rshift(f, n):
    """Efficiently divide f by x**n. """
    if not n:
        return f, []
    else:
        return f[:-n], f[-n:]

def gf_pow(f, n, p):
    """Computes f**n over GF(p)[x] using repeated squaring. """
    assert n >= 0

    if not n:
        return [1]
    elif n == 1:
        return f
    elif n == 2:
        return gf_sqr(f, p)

    h = [1]

    while True:
        if n & 1:
            h = gf_mul(h, f, p)
            n -= 1

        n >>= 1

        if not n:
            break

        f = gf_sqr(f, p)

    return h

def gf_pow_mod(f, n, g, p):
    """Computes f**n over GF(p)[x]/(g) using repeated squaring.

       Given polynomials f and g over GF(p)[x] and  a non-negative
       integer n, efficiently computes f**n (mod g) i.e. remainder
       from division f**n by g using repeated squaring algorithm.

       For more details on the implemented algorithm refer to:

       [1] J. von zur Gathen, J. Gerhard, Modern Computer Algebra,
           First Edition, Cambridge University Press, 1999, pp. 69

    """
    assert n >= 0

    if not n:
        return [1]
    elif n == 1:
        return gf_rem(f, g, p)
    elif n == 2:
        return gf_rem(gf_sqr(f, p), g, p)

    h = [1]

    while True:
        if n & 1:
            h = gf_mul(h, f, p)
            h = gf_rem(h, g, p)
            n -= 1

        n >>= 1

        if not n:
            break

        f = gf_sqr(f, p)
        f = gf_rem(f, g, p)

    return h

def gf_gcd(f, g, p):
    """Euclidean Algorithm over GF(p)[x]. """
    while g:
        f, g = g, gf_rem(f, g, p)

    return gf_monic(f, p)[1]

def gf_gcdex(f, g, p):
    """Extended Euclidean Algorithm over GF(p)[x].

       Given polynomials f and g over GF(p)[x],  computes polynomials
       s, t and h, such that h = gcd(f, g) and s*f + t*g = h. Typical
       application of EEA is solving polynomial diophantine equations.

       Consider polynomials f = (x + 7) (x + 1), g = (x + 7) (x**2 + 1)
       over GF(11)[x]. Application of Extended Euclidean Algorithm gives:

       >>> from sympy.polys.galoispolys import gf_gcdex, gf_mul, gf_add
       >>> s, t, g = gf_gcdex([1,8,7], [1,7,1,7], 11)
       >>> s, t, g
       ([5, 6], [6], [1, 7])

       As result we obtained polynomials s = 5*x + 6 and t = 6, and
       additionally gcd(f, g) = x + 7. This is correct because:

       >>> S = gf_mul(s,   [1,8,7], 11)
       >>> T = gf_mul(t, [1,7,1,7], 11)

       >>> gf_add(S, T, 11) == [1, 7]
       True

       For more details on the implemented algorithm refer to:

       [1] J. von zur Gathen, J. Gerhard, Modern Computer Algebra,
           First Edition, Cambridge University Press, 1999, pp. 46

    """
    if not (f or g):
        return [1], [], []

    p0, r0 = gf_monic(f, p)
    p1, r1 = gf_monic(g, p)

    if not f:
        return [], [zp_inv(p1, p)], r1
    if not g:
        return [zp_inv(p0, p)], [], r0

    s0, s1 = [zp_inv(p0, p)], []
    t0, t1 = [], [zp_inv(p1, p)]

    while True:
        Q, R = gf_div(r0, r1, p)

        if not R:
            break

        (lc, r1), r0 = gf_monic(R, p), r1

        inv = zp_inv(lc, p)

        s = gf_sub_mul(s0, s1, Q, p)
        t = gf_sub_mul(t0, t1, Q, p)

        s1, s0 = gf_mul_const(s, inv, p), s1
        t1, t0 = gf_mul_const(t, inv, p), t1

    return s1, t1, r1

def gf_monic(f, p):
    """Returns LC and a monic polynomial over GF(p)[x]."""
    if not f:
        return 0, []
    else:
        LC = f[0]

        if LC == 1:
            return 1, f[:]
        else:
            return LC, gf_div_const(f, LC, p)

def gf_diff(f, p):
    """Differentiate polynomial over GF(p)[x]. """
    df = gf_degree(f)

    h, n = [0]*df, df

    for coeff in f[:-1]:
        coeff *= n
        coeff %= p

        if coeff:
            h[df-n] = coeff

        n -= 1

    return gf_strip(h)

def gf_eval(f, x, p):
    """Evaluate f(x) over GF(p) using Horner scheme. """
    result = 0

    for a in f:
        result *= x
        result += a
        result %= p

    return result

def gf_multi_eval(f, A, p):
    """Evaluate f(x) for x in { a_1, ..., a_n }. """
    return [ gf_eval(f, a, p) for a in A ]

def gf_compose(g, h, p):
    """Compute polynomial composition g(h) over GF(p)[x]. """
    if not g:
        return []

    comp = [g[0]]

    for a in g[1:]:
        comp = gf_mul(comp, h, p)
        comp = gf_add_const(comp, a, p)

    return comp

def gf_compose_mod(g, h, f, p):
    """Compute polynomial composition g(h) over GF(p)[x]/(f). """
    if not g:
        return []

    comp = [g[0]]

    for a in g[1:]:
        comp = gf_mul(comp, h, p)
        comp = gf_add_const(comp, a, p)
        comp = gf_rem(comp, f, p)

    return comp

def gf_trace_map(a, b, c, n, f, p):
    """Compute polynomial trace map over GF(p)[x]/(f).

       Given polynomial f over GF(p)[x],  polynomials a, b, c over quotient
       ring GF(p)[x]/(f) such that b = c**t (mod f) for some positive power
       t of p and a positive integer n, returns a mapping:

            a -> a**t**n, a + a**t + a**t**2 + ... + a**t**n (mod f)

       Typically, especially in factorization context, b = x**p mod f and
       c = x mod f. This way we can efficiently compute trace polynomials
       in equal degree factorization routine, much faster than with other
       methods, like iterated Frobenius method, for large degrees.

       For more details on the implemented algorithm refer to:

       [1] J. von zur Gathen, V. Shoup, Computing Frobenius Maps
           and Factoring Polynomials, ACM Symposium on Theory of
           Computing, 1992, pp. 187-224

    """
    u = gf_compose_mod(a, b, f, p)
    v = b

    if n & 1:
        U = gf_add(a, u, p)
        V = b
    else:
        U = a
        V = c

    n >>= 1

    while n:
        u = gf_add(u, gf_compose_mod(u, v, f, p), p)
        v = gf_compose_mod(v, v, f, p)

        if n & 1:
            U = gf_add(U, gf_compose_mod(u, V, f, p), p)
            V = gf_compose_mod(v, V, f, p)

        n >>= 1

    return gf_compose_mod(a, V, f, p), U

def gf_random(n, p, monic=True):
    """Generate random polynomial over GF(p)[x] of degree exactly n. """
    if not monic:
        LC = [int(uniform(1, p))]
    else:
        LC = [1]

    return LC + [ int(uniform(0, p)) for i in xrange(0, n) ]

def gf_irreducible_p(f, p):
    """Test irreducibility of f over GF(p)[x] using deterministic method. """
    raise NotImplementedError

def gf_irreducible(n, p, monic=True):
    """Generate random irreducible polynomial of degree n over GF(p)[x]. """
    while True:
        f = gf_random(n, p, monic)

        H = h = gf_pow_mod([1, 0], p, f, p)

        for i in xrange(1, n/2 + 1):
            g = gf_sub(h, [1, 0], p)

            if gf_gcd(f, g, p) == [1]:
                h = gf_compose_mod(h, H, f, p)
            else:
                break
        else:
            return f

def gf_sqf_p(f, p):
    """Returns True if f is square-free over GF(p)[x]. """
    LC, f = gf_monic(f, p)

    if not f:
        return True
    else:
        return gf_gcd(f, gf_diff(f, p), p) == [1]

def gf_sqf(f, p):
    """Returns square-free decomposition of a GF(p)[x] polynomial.

       Given a polynomial f over GF(p)[x], returns the leading coefficient
       of f and a square-free decomposition f_1**e_1 f_2**e_2 ... f_k**e_k
       such that all f_i are monic polynomials,  (f_i, f_j) for i != j are
       co-prime and e_1 ... e_k are given in increasing order. All trivial
       terms (i.e. f_i = 1) aren't included in the output.

       To check if a polynomial is square-free use gf_sqf_p() function.

       Consider polynomial f = x**11 + 1 over GF(11)[x]:

       >>> from sympy.polys.galoispolys import gf_from_dict, gf_diff
       >>> f = gf_from_dict({11: 1, 0: 1}, 11)

       Note that f'(x) = 0:

       >>> gf_diff(f, 11)
       []

       This phenomenon does not happen in characteristic zero. However
       we can still compute square-free decomposition of f using gf_sqf():

       >>> from sympy.polys.galoispolys import gf_sqf, gf_pow
       >>> gf_sqf(f, 11)
       (1, [([1, 1], 11)])

       We obtained factorization f = (x + 1)**11. This is correct because:

       >>> gf_pow([1, 1], 11, 11) == f
       True

       For more details on the implemented algorithm refer to:

       [1] K. Geddes, S. Czapor, G. Labahn, Algorithms for Computer
           Algebra, First Edition, Springer, 1992, pp. 343-347

    """
    n, sqf, factors = 1, False, []

    LC, f = gf_monic(f, p)

    if gf_degree(f) < 1:
        return LC, []

    while True:
        F = gf_diff(f, p)

        if F != []:
            g = gf_gcd(f, F, p)
            h = gf_quo(f, g, p)

            i = 1

            while h != [1]:
                G = gf_gcd(g, h, p)
                H = gf_quo(h, G, p)

                if gf_degree(H) > 0:
                    factors.append((H, i*n))

                g, h, i = gf_quo(g, G, p), G, i+1

            if g == [1]:
                sqf = True
            else:
                f = g

        if not sqf:
            d = gf_degree(f) // p

            for i in xrange(0, d+1):
                f[i] = f[i*p]

            f, n = f[:d+1], n*p
        else:
            break

    return LC, factors

def gf_ddf_zassenhaus(f, p):
    """Cantor-Zassenhaus: Deterministic Distinct Degree Factorization

       Given a monic square-free polynomial f in GF(p)[x], computes
       partial distinct degree factorization f_1 ... f_d of f where
       deg(f_i) != deg(f_j) for i != j. The result is returned as a
       list of pairs (f_i, e_i) where deg(f_i) > 0 and e_i > 0 is
       an argument to the equal degree factorization routine.

       Consider polynomial x**15 - 1 over GF(11)[x].

       >>> from sympy.polys.galoispolys import gf_from_dict, gf_ddf_zassenhaus
       >>> f = gf_from_dict({15: 1, 0: -1}, 11)

       Distinct degree factorization gives:

       >>> gf_ddf_zassenhaus(f, 11)
       [([1, 0, 0, 0, 0, 10], 1), ([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1], 2)]

       which means x**15 - 1 = (x**5 - 1) (x**10 + x**5 + 1). To obtain
       factorization into irreducibles,  use equal degree factorization
       procedure (EDF) with each of the factors.

       For more details on the implemented algorithm refer to:

       [1] J. von zur Gathen, J. Gerhard, Modern Computer Algebra,
           First Edition, Cambridge University Press, 1999, pp. 356

       [2] K. Geddes, S. Czapor, G. Labahn, Algorithms for Computer
           Algebra, First Edition, Springer, 1992, pp. 368-371

    """
    i, g, factors = 1, [1, 0], []

    while 2*i <= gf_degree(f):
        g = gf_pow_mod(g, p, f, p)
        h = gf_gcd(f, gf_sub(g, [1, 0], p), p)

        if h != [1]:
            factors.append((h, i))

            f = gf_quo(f, h, p)
            g = gf_rem(g, f, p)

        i += 1

    if f != [1]:
        return factors + [(f, gf_degree(f))]
    else:
        return factors

def gf_edf_zassenhaus(f, n, p):
    """Cantor-Zassenhaus: Probabilistic Equal Degree Factorization

       Given a monic square-free polynomial f in GF(p)[x] and integer n
       such that n divides deg(f),  returns all irreducible factors f_1
       ... f_d of f, each of degree n. This is a complete factorization
       over Galois fields.

       Consider square-free polynomial f = x**3 + x**2 + x + 1 over
       GF(5)[x]. Lets compute its irreducible factors of degree one:

       >>> from sympy.polys.galoispolys import gf_edf_zassenhaus
       >>> gf_edf_zassenhaus([1,1,1,1], 1, 5)
       [[1, 1], [1, 2], [1, 3]]

       For more details on the implemented algorithm refer to:

       [1] J. von zur Gathen, J. Gerhard, Modern Computer Algebra,
           First Edition, Cambridge University Press, 1999, pp. 358

       [2] K. Geddes, S. Czapor, G. Labahn, Algorithms for Computer
           Algebra, First Edition, Springer, 1992, pp. 371-373

    """
    factors = [f]

    if gf_degree(f) <= n:
        return factors

    N = gf_degree(f) // n

    while len(factors) < N:
        r = gf_random(2*n-1, p)

        if p == 2:
            h = r

            for i in xrange(0, 2**(n*N-1)):
                r = gf_pow_mod(r, 2, f, p)
                h = gf_add(h, r, p)

            g = gf_gcd(f, h, p)
        else:
            h = gf_pow_mod(r, (p**n-1) // 2, f, p)
            g = gf_gcd(f, gf_sub_const(h, 1, p), p)

        if g != [1] and g != f:
            factors = gf_edf_zassenhaus(g, n, p) \
                    + gf_edf_zassenhaus(gf_quo(f, g, p), n, p)

    def compare(f_a, f_b):
        i = len(f_a) - len(f_b)

        if not i:
            return cmp(f_a, f_b)
        else:
            return i

    return sorted(factors, compare)

def gf_ddf_shoup(f, p):
    """Kaltofen-Shoup: Deterministic Distinct Degree Factorization

       Given a monic square-free polynomial f in GF(p)[x], computes
       partial distinct degree factorization f_1 ... f_d of f where
       deg(f_i) != deg(f_j) for i != j. The result is returned as a
       list of pairs (f_i, e_i) where deg(f_i) > 0 and e_i > 0 is
       an argument to the equal degree factorization routine.

       This algorithm is an improved version of Zassenhaus algorithm
       for large deg(f) and modulus p (especially for deg(f) ~ lg(p)).

       For more details on the implemented algorithm refer to:

       [1] E. Kaltofen, V. Shoup, Subquadratic-time Factoring of
           Polynomials over Finite Fields, Mathematics of Compu-
           tation, Volume 67, Issue 223, 1998, pp. 1179-1197

       [2] V. Shoup, A New Polynomial Factorization Algorithm and
           its Implementation, Journal of Symbolic Computation,
           Volume 20, Issue 4, 1995, pp. 363-397

       [3] J. von zur Gathen, V. Shoup, Computing Frobenius Maps
           and Factoring Polynomials, ACM Symposium on Theory of
           Computing, 1992, pp. 187-224

    """
    n = gf_degree(f)
    k = int(ceil(sqrt(n//2)))

    h = gf_pow_mod([1, 0], p, f, p)

    U = [[1,0], h] + [0]*(k-1)

    for i in xrange(2, k+1):
        U[i] = gf_compose_mod(U[i-1], h, f, p)

    h, U = U[k], U[:k]
    V = [h] + [0]*(k-1)

    for i in xrange(1, k):
        V[i] = gf_compose_mod(V[i-1], h, f, p)

    factors = []

    for i, v in enumerate(V):
        h, j = [1], k-1

        for u in U:
            g = gf_sub(v, u, p)
            h = gf_mul(h, g, p)
            h = gf_rem(h, f, p)

        g = gf_gcd(f, h, p)
        f = gf_quo(f, g, p)

        for u in reversed(U):
            h = gf_sub(v, u, p)
            F = gf_gcd(g, h, p)

            if F != [1]:
                factors.append((F, k*(i+1)-j))

            g, j = gf_quo(g, F, p), j-1

    if f != [1]:
        factors.append((f, gf_degree(f)))

    return factors

def gf_edf_shoup(f, n, p):
    """Gathen-Shoup: Probabilistic Equal Degree Factorization

       Given a monic square-free polynomial f in GF(p)[x] and integer n
       such that n divides deg(f),  returns all irreducible factors f_1
       ... f_d of f, each of degree n. This is a complete factorization
       over Galois fields.

       This algorithm is an improved version of Zassenhaus algorithm
       for large deg(f) and modulus p (especially for deg(f) ~ lg(p)).

       For more details on the implemented algorithm refer to:

       [1] V. Shoup, A Fast Deterministic Algorithm for Factoring
           Polynomials over Finite Fields of Small Characteristic,
           In Proceedings of International Symposium on Symbolic
           and Algebraic Computation, 1991, pp. 14-21

       [2] J. von zur Gathen, V. Shoup, Computing Frobenius Maps
           and Factoring Polynomials, ACM Symposium on Theory of
           Computing, 1992, pp. 187-224

    """
    N = gf_degree(f)

    if not N:
        return []
    if N <= n:
        return [f]

    factors, x = [f], [1, 0]

    r = gf_random(N-1, p)

    h = gf_pow_mod(x, p, f, p)
    H = gf_trace_map(r, h, x, n-1, f, p)[1]

    if p == 2:
        h1 = gf_gcd(f, H, p)
        h2 = gf_quo(f, h1, p)

        factors = gf_edf_shoup(h1, n, p) \
                + gf_edf_shoup(h2, n, p)
    else:
        h = gf_pow_mod(H, (p-1)//2, f, p)

        h1 = gf_gcd(f, h, p)
        h2 = gf_gcd(f, gf_sub_const(h, 1, p), p)
        h3 = gf_quo(f, gf_mul(h1, h2, p), p)

        factors = gf_edf_shoup(h1, n, p) \
                + gf_edf_shoup(h2, n, p) \
                + gf_edf_shoup(h3, n, p)

    def compare(f_a, f_b):
        i = len(f_a) - len(f_b)

        if not i:
            return cmp(f_a, f_b)
        else:
            return i

    return sorted(factors, compare)

_ddf_methods = {
    'zassenhaus' : gf_ddf_zassenhaus,
    'shoup'      : gf_ddf_shoup,
}

def gf_ddf(f, p, **flags):
    """Distinct Degree Factorization over GF(p)[x]. """
    method = flags.get('ddf')

    if method is None:
        method = flags.get('method')

    try:
        if method is not None:
            ddf_method = _ddf_methods[method.lower()]
        else:
            # TODO: use cross-over
            ddf_method = gf_ddf_zassenhaus
    except KeyError:
        raise ValueError("'%s' is not a valid DDF method" % method)

    try:
        return ddf_method(f, p)
    except NotImplementedError:
        return gf_ddf_zassenhaus(f, p)

_edf_methods = {
    'zassenhaus' : gf_edf_zassenhaus,
    'shoup'      : gf_edf_shoup,
}

def gf_edf(f, n, p, **flags):
    """Equal Degree Factorization over GF(p)[x]. """
    method = flags.get('edf')

    if method is None:
        method = flags.get('method')

    try:
        if method is not None:
            edf_method = _edf_methods[method.lower()]
        else:
            # TODO: use cross-over
            edf_method = gf_edf_zassenhaus
    except KeyError:
        raise ValueError("'%s' is not a valid EDF method" % method)

    try:
        return edf_method(f, n, p)
    except NotImplementedError:
        return gf_edf_zassenhaus(f, n, p)

def gf_factor(f, p, **flags):
    """Factor (non square-free) polynomials over GF(p)[x].

       Given a possibly non square-free polynomial f over GF(p)[x], returns
       complete factorization into irreducibles f_1(x)**e_1 f_2(x)**e_2 ...
       f_d(x)**e_d of f(x), where each f_i is monic and gcd(f_i, f_j) == 1,
       for i != j. The result is given as a tuple consisting of the leading
       coefficient of f and a list of factors with their multiplicities.

       The algorithm proceeds by first computing square-free decomposition
       of f and then iteratively factoring each of the square-free factors.

       Consider a non square-free polynomial f = (7*x + 1) (x + 2)**2 over
       GF(11)[x]. We obtain its factorization into irreducibles as follows:

       >>> from sympy.polys.galoispolys import gf_factor
       >>> gf_factor([5, 2, 7, 2], 11)
       (5, [([1, 2], 1), ([1, 8], 2)])

       We arrived with factorization f = 5 (x + 2) (x + 8)**2.  We didn't
       recover exact form of the input polynomial because we requested to
       get monic factors of f and its leading coefficient separately.

       For more details on the implemented algorithm refer to:

       [1] J. von zur Gathen, J. Gerhard, Modern Computer Algebra,
           First Edition, Cambridge University Press, 1999, pp. 365

    """
    LC, f = gf_monic(f, p)

    if gf_degree(f) < 1:
        return LC, []

    factors = []

    for g, k in gf_sqf(f, p)[1]:
        n = gf_degree(g)

        if n == 0:
            continue

        for h in gf_factor_sqf(g, p)[1]:
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

    return LC, sorted(factors, compare)

def gf_factor_sqf(f, p, **flags):
    """Factor square-free polynomials over GF(p)[x]. """
    LC, f = gf_monic(f, p)

    if gf_degree(f) < 1:
        return LC, []

    factors = []

    for factor, n in gf_ddf(f, p, **flags):
        factors += gf_edf(factor, n, p, **flags)

    return LC, factors


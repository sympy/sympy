"""Univariate polynomials with coefficients in the ring of integers. """

from sympy.polys.galoispolys import gf_from_int_poly, gf_to_int_poly, gf_degree, \
    gf_from_dict, gf_mul, gf_quo, gf_gcd, gf_gcdex, gf_sqf_p, gf_factor_sqf

from sympy.ntheory import randprime, isprime

from sympy.utilities.iterables import subsets
from sympy.ntheory.modular import crt1, crt2
from sympy.core.numbers import igcd, igcdex

from math import floor, ceil, sqrt, log

class HeuristicGCDFailed(Exception):
    pass

def zzx_degree(f):
    """Returns leading degree of f. """
    return len(f) - 1

def zzx_LC(f):
    """Returns leading coefficient of f. """
    if not f:
        return 0
    else:
        return f[0]

def zzx_TC(f):
    """Returns trailing coefficient of f. """
    if not f:
        return 0
    else:
        return f[-1]

def zzx_nth(f, n):
    """Returns n-th coefficient of f. """
    return f[zzx_degree(f)-n]

def zzx_strip(f):
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

def zzx_from_dict(f):
    """Create Z[x] polynomial from a dict. """
    n, h = max(f.iterkeys()), []

    for k in xrange(n, -1, -1):
        h.append(f.get(k, 0))

    return zzx_strip(h)

def zzx_to_dict(f):
    """Convert Z[x] polynomial to a dict. """
    n, result = zzx_degree(f), {}

    for i in xrange(0, n+1):
        if f[n-i]:
            result[i] = f[n-i]

    return result

def zzx_add_term(f, c, k):
    """Add c*x**k to f over Z[x]. """
    if not c:
        return f

    n = len(f)
    m = n-k-1

    if k == n-1:
        return zzx_strip([f[0]+c] + f[1:])
    else:
        if k >= n:
            return [c] + [0] * (k-n) + f
        else:
            return f[:m] + [f[m]+c] + f[m+1:]

def zzx_sub_term(f, c, k):
    """Subtract c*x**k from f over Z[x]. """
    if not c:
        return f

    n = len(f)
    m = n-k-1

    if k == n-1:
        return zzx_strip([f[0]-c] + f[1:])
    else:
        if k >= n:
            return [-c] + [0] * (k-n) + f
        else:
            return f[:m] + [f[m]-c] + f[m+1:]

def zzx_mul_term(f, c, k):
    """Multiply f with c*x**k over Z[x]. """
    if not c or not f:
        return []
    else:
        return [ c * cc for cc in f ] + [0] * k

def zzx_abs(f):
    """Make all coefficients positive. """
    return [ abs(coeff) for coeff in f ]

def zzx_neg(f):
    """Negate a polynomial over Z[x]. """
    return [ -coeff for coeff in f ]

def zzx_add(f, g):
    """Add polynomials over Z[x]. """
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

def zzx_sub(f, g):
    """Subtract polynomials over Z[x]. """
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

def zzx_add_mul(f, g, h):
    """Returns f + g*h where f, g, h in Z[x]. """
    return zzx_add(f, zzx_mul(g, h))

def zzx_sub_mul(f, g, h):
    """Returns f - g*h where f, g, h in Z[x]. """
    return zzx_sub(f, zzx_mul(g, h))

def zzx_mul(f, g):
    """Multiply polynomials over Z[x]. """
    df = zzx_degree(f)
    dg = zzx_degree(g)

    dh = df + dg
    h = [0]*(dh+1)

    for i in xrange(0, dh+1):
        coeff = 0

        for j in xrange(max(0, i-dg), min(i, df)+1):
            coeff += f[j]*g[i-j]

        h[i] = coeff

    return zzx_strip(h)

def zzx_sqr(f):
    """Square polynomials over Z[x]. """
    df = zzx_degree(f)

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

        h[i] = coeff

    return zzx_strip(h)

def zzx_div(f, g):
    """Returns quotient and remainder over Z[x]. """
    df = zzx_degree(f)
    dg = zzx_degree(g)

    if not g:
        raise ZeroDivisionError("polynomial division")
    elif df < dg:
        return [], f

    q, r = [], f

    while True:
        deg_r = zzx_degree(r)
        deg_g = zzx_degree(g)

        if deg_r < deg_g:
            break

        lc_r = zzx_LC(r)
        lc_g = zzx_LC(g)

        if lc_r % lc_g != 0:
            break

        c = lc_r // lc_g
        k = deg_r - deg_g

        q = zzx_add_term(q, c, k)
        h = zzx_mul_term(g, c, k)
        r = zzx_sub(r, h)

    return q, r

def zzx_quo(f, g):
    """Returns polynomial remainder over Z[x]. """
    return zzx_div(f, g)[0]

def zzx_rem(f, g):
    """Returns polynomial remainder over Z[x]. """
    return zzx_div(f, g)[1]

def zzx_max_norm(f):
    """Returns maximum norm of a polynomial. """
    if not f:
        return 0
    else:
        return max(zzx_abs(f))

def zzx_l1_norm(f):
    """Returns l1 norm of a polynomial. """
    if not f:
        return 0
    else:
        return sum(zzx_abs(f))

def zzx_diff(f):
    """Differentiate polynomial over Z[x]. """
    n, deriv = zzx_degree(f), []

    for coeff in f[:-1]:
        deriv.append(n*coeff)
        n -= 1

    return deriv

def zzx_eval(f, x):
    """Evaluate f(x) using Horner scheme. """
    result = 0

    for a in f:
        result *= x
        result += a

    return result

def zzx_trunc(f, m):
    """Reduce Z[x] polynomial modulo m. """
    g = []

    for coeff in f:
        coeff %= m

        if coeff > m // 2:
            g.append(coeff - m)
        else:
            g.append(coeff)

    return zzx_strip(g)

def zzx_content(f):
    """Returns integer GCD of coefficients. """
    cont = 0

    for coeff in f:
        cont = igcd(cont, coeff)

        if cont == 1:
            break

    return cont

def zzx_primitive(f):
    """Divides all coefficients by content. """
    cont = zzx_content(f)

    if cont == 1:
        return 1, f
    else:
        return cont, [ coeff // cont for coeff in f ]

def zzx_sqf_part(f):
    """Returns square-free part of a polynomial over Z[x]. """
    quo = zzx_quo(f, zzx_gcd(f, zzx_diff(f)))
    return zzx_primitive(quo)[1]

def zzx_gcd(f, g, **flags):
    """Polynomial GCD over Z[x]. """
    if flags.get("heu", True):
        try:
            return zzx_heu_gcd(f, g)[0]
        except HeuristicGCDFailed:
            pass

    return zzx_mod_gcd(f, g)[0]

def zzx_heu_gcd(f, g, **flags):
    """Heuristic polynomial GCD over Z[x].

       Given univariate polynomials f and g over Z[x], returns their GCD
       and cofactors, i.e. polynomials h, cff and cfg such that:

              h = gcd(f, g), cff = quo(f, h) and cfg = quo(g, h)

       The algorithm is purely heuristic which means it may fail to compute
       the GCD. In this case an exception is raised. However it seems that
       this algorithm fails very rarely. In case of failure one can switch
       to modular GCD algorithm (zzx_mod_gcd). In high-level code, instead
       of direct application of heuristic GCD method, use zzx_gcd function
       which will do the switch automatically when necessary.

       The algorithm computes the polynomial GCD by evaluating polynomials
       f and g at certain points and computing (fast) integer GCD of those
       evaluations. The polynomial GCD is recovered from the integer image
       by interpolation. The final step is to verify correctness of the
       result. By default the algorithm will perform at most six tries
       to compute the GCD.

       Note this method can be easily adapted to multivariate case.  The
       only things to change is to modify interpolation function and add
       recursive invocation of the heuristic gcd algorithm.

       For more details on the implemented algorithm refer to:

       [1] Hsin-Chao Liao, R. Fateman, Evaluation of the heuristic polynomial
           GCD, International Symposium on Symbolic and Algebraic Computation
           (ISSAC), ACM Press, Montreal, Quebec, Canada, 1995, pp. 240--247

    """
    def sp_interpolate(h, c):
        f = []

        while h:
            rem = h % c

            if rem > c // 2:
                rem -= c

            f.insert(0, rem)
            h = (h-rem) // c

        return zzx_strip(f)

    df = zzx_degree(f)
    dg = zzx_degree(g)

    cf = zzx_content(f)
    cg = zzx_content(g)

    h = igcd(cf, cg)

    f = [ c // h for c in f ]
    g = [ c // h for c in g ]

    if df <= 0 or dg <= 0:
        return zzx_strip([h]), f, g
    else:
        gcd = h

    f_norm = zzx_max_norm(f)
    g_norm = zzx_max_norm(g)

    B = 2*min(f_norm, g_norm) + 29

    x = max(min(B, 99*sqrt(B)),
            2*min(f_norm // abs(zzx_LC(f)),
                  g_norm // abs(zzx_LC(g))) + 2)

    for i in xrange(0, 6):
        length_x = int(floor(log(x, 2)) + 1)

        if length_x * max(df, dg) > 4000:
            raise HeuristicGCDFailed

        ff = zzx_eval(f, x)
        gg = zzx_eval(g, x)

        h = igcd(ff, gg)

        cff = ff // h
        cfg = gg // h

        h = sp_interpolate(h, x)
        h = zzx_primitive(h)[1]

        cff_, r = zzx_div(f, h)

        if not r:
            cfg_, r = zzx_div(g, h)

            if not r:
                h = zzx_mul_term(h, gcd, 0)
                return h, cff_, cfg_

        cff = sp_interpolate(cff, x)

        h, r = zzx_div(f, cff)

        if not r:
            cfg_, r = zzx_div(g, h)

            if not r:
                h = zzx_mul_term(h, gcd, 0)
                return h, cff, cfg_

        cfg = sp_interpolate(cfg, x)

        h, r = zzx_div(g, cfg)

        if not r:
            cff_, r = zzx_div(f, h)

            if not r:
                h = zzx_mul_term(h, gcd, 0)
                return h, cff_, cfg

        x = int(2.7319*x*sqrt(sqrt(x)))

    raise HeuristicGCDFailed

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
    n = zzx_degree(f)
    m = zzx_degree(g)

    cf = zzx_content(f)
    cg = zzx_content(g)

    h = igcd(cf, cg)

    f = [ c // h for c in f ]
    g = [ c // h for c in g ]

    if n <= 0 or m <= 0:
        return zzx_strip([h]), f, g
    else:
        gcd = h

    A = max(zzx_abs(f) + zzx_abs(g))
    b = igcd(zzx_LC(f), zzx_LC(g))

    B = int(ceil(2**n*A*b*sqrt(n + 1)))
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

                if (p in primes) and (b % p == 0):
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
            C = [ b * zzx_nth(hh[p], i) for p in primes ]
            c = crt2(primes, C, crt_mm, crt_e, crt_s, True)

            H.insert(0, c)

        H = zzx_strip(H)

        for i in xrange(0, zzx_degree(f) - e + 1):
            C = [ zzx_nth(fff[p], i) for p in primes ]
            c = crt2(primes, C, crt_mm, crt_e, crt_s, True)

            F.insert(0, c)

        for i in xrange(0, zzx_degree(g) - e + 1):
            C = [ zzx_nth(ggg[p], i) for p in primes ]
            c = crt2(primes, C, crt_mm, crt_e, crt_s, True)

            G.insert(0, c)

        H_norm = zzx_l1_norm(H)

        F_norm = zzx_l1_norm(F)
        G_norm = zzx_l1_norm(G)

        if H_norm*F_norm <= B and H_norm*G_norm <= B:
            break

    return zzx_mul_term(H, gcd, 0), F, G

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
    lc = zzx_LC(f)

    if r == 1:
        F = zzx_mul_term(f, igcdex(lc, p**l)[0], 0)
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
    """Factor square-free polynomials over Z[x]. """
    n = zzx_degree(f)

    if n == 1:
        return [f]

    A = zzx_max_norm(f)
    b = zzx_LC(f)
    B = abs(int(sqrt(n+1)*2**n*A*b))
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
                b = zzx_LC(f)

                break
        else:
            s += 1

    return factors + [f]

def zzx_factor(f):
    """Factor (non square-free) polynomials over Z[x].

       Given a univariate polynomial f over Z[x] computes its complete
       factorization f_1, ..., f_n into irreducibles over integers:

                    f = content(f) f_1**k_1 ... f_n**k_n

       The factorization is computed by reducing the input polynomial
       into a primitive square-free polynomial and factoring it using
       Zassenhaus algorithm. Trial division is used to recover the
       multiplicities of factors.

       The result is returned as a tuple consisting of:

                 (content(f), [(f_1, k_1), ..., (f_n, k_n))

       Consider polynomial f = x**4 - 1:

       >>> zzx_factor([2, 0, 0, 0, -2])
       (2, [([1, -1], 1), ([1, 1], 1), ([1, 0, 1], 1)])

       In result we got the following factorization:

                    f = 2 (x - 1) (x + 1) (x**2 + 1)

       Note that this is a complete factorization over integers,
       however over Gaussian integers we can factor the last term.

       For more details on the implemented algorithm refer to:

       [1] J. von zur Gathen, J. Gerhard, Modern Computer Algebra,
           First Edition, Cambridge University Press, 1999, pp. 427
    """
    cont, g = zzx_primitive(f)

    if zzx_degree(g) < 1:
        return cont, []

    if zzx_LC(g) < 0:
        g = zzx_neg(g)
        cont = -cont

    g = zzx_sqf_part(g)

    factors = []

    for h in zzx_zassenhaus(g):
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


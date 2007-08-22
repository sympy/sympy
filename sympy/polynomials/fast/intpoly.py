"""Univariate polynomials with integer coefficients."""

import math

from sympy import ntheory
from sympy.polynomials.fast import modint, sparse_poly, gfpoly

class IntPoly(sparse_poly.SparsePolynomial):
    def primitive(self):
        content = reduce(modint.gcd, self.coeffs.itervalues(), 0)
        result_dict = {}
        for e, c in self.coeffs.iteritems():
            result_dict[e] = c/content
        return content, IntPoly(result_dict)

    def mod_int(self, m, symmetric=False):
        result_dict = {}
        for e, c in self.coeffs.iteritems():
            cc = c % m
            if cc:
                if symmetric and cc > m/2:
                    result_dict[e] = cc - m
                else:
                    result_dict[e] = cc
        return IntPoly(result_dict)

# Division algorithms:

def div(f, g):
    """Division with remainder over the integers."""
    q = IntPoly()
    r = f
    if not g:
        return q, r
    while r.degree >= g.degree and not (r[r.degree] % g[g.degree]):
        quot = IntPoly({r.degree - g.degree: r[r.degree]/g[g.degree]})
        q += quot
        r -= quot*g
    return q, r

def gcd_small_primes(f, g):
    """Modular small primes version for primitive polynomials."""
    if f.degree < g.degree:
        f, g = g, f
    if not g:
        return f
    if g.degree == 0:
        c = reduce(modint.gcd, f.coeffs.itervalues(), g[0])
        return IntPoly({0: c})

    n = f.degree
    A = max([abs(c) for c in f.coeffs.itervalues()]
            + [abs(c) for c in g.coeffs.itervalues()])
    b = modint.gcd(f[f.degree], g[g.degree])
    B = int(math.ceil(2**n*A*b*math.sqrt(n+1)))
    k = int(math.ceil(2*math.log((n+1)**n*b*A**(2*n), 2)))
    l = int(math.ceil(math.log(2*B + 1, 2)))
    # TODO: the minimum is needed for very small polynomials?
    prime_border = max(int(math.ceil(2*k*math.log(k))), 51)

    while True:
        while True:
            # Choose primes.
            S = []
            while len(S) < l:
                p = ntheory.generate.randprime(2, prime_border + 1)
                if (p not in S) and (b % p): # p doesn't divide b.
                    S.append(p)

            # Call the modular gcd.
            v, ff, gg = {}, {}, {} 
            for p in S:
                poly_type = gfpoly.GFPolyFactory(p)
                ff[p] = poly_type.from_int_dict(f.coeffs)
                gg[p] = poly_type.from_int_dict(g.coeffs)
                v[p] = gfpoly.gcd(ff[p], gg[p])

            e = min([v[p].degree for p in S])
            unlucky = []
            for p in S:
                if v[p].degree != e:
                    unlucky.append(p)
                    S.remove(p)
                    del v[p]
                    del ff[p]
                    del gg[p]

            if len(S) < l/2: # Forget all primes.
                continue
                    
            # Replace the unlucky primes.
            while len(S) < l:
                p = ntheory.generate.randprime(2, prime_border)
                if (p in unlucky) or (p in S) or (b % p == 0): 
                    continue
                poly_type = gfpoly.GFPolyFactory(p)
                ff[p] = poly_type.from_int_dict(f.coeffs)
                gg[p] = poly_type.from_int_dict(g.coeffs)
                v[p] = gfpoly.gcd(ff[p], gg[p])
                if v[p].degree == e:
                    S.append(p)
                else:
                    unlucky.append(p)
                    del v[p]
                    del ff[p]
                    del gg[p]
            break # The primes are good.
        
        fff, ggg = {}, {}
        for p in S:
            fff[p], r = gfpoly.div(ff[p], v[p])
            assert not r
            ggg[p], r = gfpoly.div(gg[p], v[p])
            assert not r
        w_dict, fff_dict, ggg_dict = {}, {}, {}
        crt_mm, crt_e, crt_s = modint.crt1(S)
        for i in xrange(0, e+1):
            C = [int(v[p][i]*v[p].__class__.coeff_type(b)) for p in S]
            c = modint.crt2(S, C, crt_mm, crt_e, crt_s, True)
            if c:
                w_dict[i] = c
        for i in xrange(0, f.degree - e + 1):
            c = modint.crt2(S, [int(fff[p][i]) for p in S], crt_mm,
                            crt_e, crt_s, True)
            if c:
                fff_dict[i] = c
        for i in xrange(0, g.degree - e + 1):
            c = modint.crt2(S, [int(ggg[p][i]) for p in S], crt_mm,
                            crt_e, crt_s, True)
            if c:
                ggg_dict[i] = c
        w_norm = sum([abs(c) for c in w_dict.itervalues()])
        fff_norm = sum([abs(c) for c in fff_dict.itervalues()])
        ggg_norm = sum([abs(c) for c in ggg_dict.itervalues()])
        if w_norm*fff_norm <= B and w_norm*ggg_norm <= B:
            break

    content, result =  IntPoly(w_dict).primitive()
    return result

def gcd_heuristic(f, g):
    """Heuristic gcd for primitive univariate polynomials."""
    def reconstruct_poly(u, c):
        result_dict = {}
        i = 0
        while c:
            rem = c % u
            if rem:
                if rem > u/2:
                    rem -= u
                result_dict[i] = rem
                c -= rem
            i += 1
            c /= u
        return IntPoly(result_dict)

    A = max([abs(c) for c in f.coeffs.itervalues()]
            + [abs(c) for c in g.coeffs.itervalues()])
    u = 4*A + 1

    while True:
        ff = f.evaluate(u)
        gg = g.evaluate(u)
        hh = modint.gcd(ff, gg)
        h = reconstruct_poly(u, hh)
        q, r = div(f, h)
        if not r:
            q, r = div(g, h)
            if not r:
                return h
        u *= 2
    
gcd = gcd_small_primes

def hensel_step(m, f, g, h, s, t):
    """One step in Hensel lifting.

    Takes an integer m and integer polynomials f, g, h, s and t as
    input, such that:
        f == g*h mod m
        s*g + t*h == 1 mod m
        lc(f) not a zero divisor mod m, h is monic
        deg(f) == deg(g) + deg(h)
        deg(s) < deg(h) and deg(t) < deg(g)

    Outputs integer polynomials gg, hh, ss and tt, such that:
        f == gg*hh mod m**2
        ss*gg + tt**hh == 1 mod m**2
    
    """
    
    mm = m**2

    e = (f - g*h).mod_int(mm, symmetric=True)
    q, r = div(s*e, h)
    q, r = q.mod_int(mm, symmetric=True), r.mod_int(mm, symmetric=True)
    gg = (g + t*e + q*g).mod_int(mm, symmetric=True)
    hh = (h + r).mod_int(mm, symmetric=True)

    b = (s*gg + t*hh - IntPoly({0: 1})).mod_int(mm, symmetric=True)
    c, d = div(s*b, hh)
    c, d = c.mod_int(mm, symmetric=True), d.mod_int(mm, symmetric=True)
    ss = (s - d).mod_int(mm, symmetric=True)
    tt = (t - t*b - c*gg).mod_int(mm, symmetric=True)
    
    return gg, hh, ss, tt

def multi_hensel_lift(p, f, f_list, l):
    """Multifactor Hensel lifting.

    Input: an integer p, an univariate integer polynomial f such that
    f's leading coefficient lc(f) is a unit mod p. Monic polynomials
    f_i that are pair-wise coprime mod p satisfying
        f = lc(f)*f_1*...*f_r mod p
    and an integer l.

    Output: monic polynomials ff_1, ..., ff_r satisfying
        f = lc(f)*ff_1*...*ff_r mod p**l
        ff_i = f_i mod p

    """
    r = len(f_list)
    lc = f[f.degree]

    if r == 1:
        lc_g, lc_s, lc_t = modint.xgcd(lc, p**l)
        return [f.scale(lc_s).mod_int(p**l, symmetric=True)]
    k = int(r/2)
    d = int(math.ceil(math.log(l, 2)))

    # Divide and conquer the factors. 
    IntModpPoly = gfpoly.GFPolyFactory(p)
    g = IntModpPoly.from_int_dict({0:lc})
    for f_i in f_list[0:k]:
        g *= IntModpPoly.from_int_dict(f_i.coeffs)
    h = IntModpPoly.from_int_dict(f_list[k].coeffs)
    for f_i in f_list[k+1:]:
        h *= IntModpPoly.from_int_dict(f_i.coeffs)
    x, s, t = gfpoly.xgcd(g, h)
    g = IntPoly(g.to_sym_int_dict())
    h = IntPoly(h.to_sym_int_dict())
    s = IntPoly(s.to_sym_int_dict())
    t = IntPoly(t.to_sym_int_dict())

    # Lift the two coprime parts.
    m = p
    for j in range(1, d+1):
        g, h, s, t = hensel_step(m, f, g, h, s, t)
        m *= m

    # Call recursively.
    return multi_hensel_lift(p, g, f_list[0:k], l) \
           + multi_hensel_lift(p, h, f_list[k:], l)

def zassenhaus(f):
    """Factors a square-free primitive polynomial.

    Returns a list of the unique factors.
    """

    def subsets(M, k):
        """Generates all k-subsets of M."""
        def recursion(result, M, k):
            if k == 0:
                 yield result
            else:
                for i, result2 in enumerate(M[0 : len(M) + 1 - k]):
                    for el in recursion(result + [result2], M[i + 1:], k - 1):
                        yield el

        for i, result in enumerate(M[0 : len(M) + 1 - k]):
            for el in recursion([result], M[i + 1:], k - 1):
                yield el

    n = f.degree
    if n == 1:
        return [f]
    A = max([abs(c) for c in f.coeffs.itervalues()])
    b = f[n]
    B = int(math.sqrt(n+1)*2**n*A*b)
    C = (n+1)**(2*n)*A**(2*n-1)
    gamma = int(math.ceil(2*math.log(C, 2)))
    prime_border = int(2*gamma*math.log(gamma)) + 1

    # Choose a prime.
    while True:
        p = ntheory.generate.randprime(2, prime_border)
        if b % p:
            poly_type = gfpoly.GFPolyFactory(p)
            ff = poly_type.from_int_dict(f.coeffs)
            gg = gfpoly.gcd(ff, ff.diff())
            if gg.degree == 0:
                break
    l = int(math.ceil(math.log(2*B + 1, p)))

    # Modular factorization.
    mod_factors = gfpoly.factor_sqf(ff)
    bb = mod_factors[0] # Leading coefficient mod p.
    h = [IntPoly(hh.to_sym_int_dict()) for hh in mod_factors[1:]]

    # Hensel lifting.
    g = multi_hensel_lift(p, f, h, l)

    # Factor combination and trial division.
    G = []
    T = range(len(g))
    s = 1
    while 2*s <= len(T):
        for S in subsets(T, s):
            gg = IntPoly({0:b})
            for i in S:
                gg *= g[i]
            gg = gg.mod_int(p**l, symmetric=True)
            hh = IntPoly({0:b})
            for i in [i for i in T if i not in S]: # T \ S
                hh *= g[i]
            hh = hh.mod_int(p**l, symmetric=True)

            gg_norm = sum([abs(c) for c in gg.coeffs.itervalues()])
            hh_norm = sum([abs(c) for c in hh.coeffs.itervalues()])
            if gg_norm*hh_norm <= B: # Found divisor
                T = [i for i in T if i not in S] # T \ S
                G.append(gg.primitive()[1])
                f = hh.primitive()[1]
                b = f[f.degree]
                break
        else: # No factors of degree s
            s += 1
    G.append(f)
    return G

def squarefree_part(f):
    """Computes the primitive squarefree part of a polynomial."""
    g = gcd(f, f.diff())
    q, r = div(f, g)
    assert not r
    return q.primitive()[1]

def factor(f):
    """Factorization of univariate integer polynomials.

    Outputs a list of factors with their multiplicities, the first
    being constant. 
    """
    
    content, pp = f.primitive()
    sqf_part = squarefree_part(pp)
    factors = zassenhaus(sqf_part)

    result = [(IntPoly({0:content}), 1)]
    # Determine multiplicities of factors.
    for ff in factors:
        mult = 0
        while True:
            q, r = div(f, ff)
            if r: # Not divisible.
                break
            else:
                mult += 1
                f = q
        result.append((ff, mult))
    
    return result

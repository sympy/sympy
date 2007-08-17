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
    print content
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
    

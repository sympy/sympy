"""Algorithms for the computation of Groebner bases"""

from sympy.polynomials.base import *
from sympy.polynomials import div_

def groebner(f, reduced=True):
    """Computes a (reduced) Groebner base for a given list of polynomials.

    Using an improved version of Buchberger's algorithm, following
    Cox, Little, O'Shea: Ideals, Varieties and Algorithms.
    """
    
    # Filter out the zero elements
    f = filter(lambda p: p.sympy_expr is not S.Zero, f)

    # empty ideal
    if len(f) == 0:
        return [Polynomial(S.Zero)]

    b = [] # Stores the unchecked combinations for s-poly's.
    s = len(f)
    for i in range(0, s-1):
        for j in range(i+1, s):
            b.append((i, j))

    while b:
        # TODO: Choose better pair: sugar?
        i, j = b[0]
        crit = False
        lcm = term_lcm(f[i].coeffs[0], f[j].coeffs[0])
        # Check if leading terms are relativly prime.
        if  lcm[1:] != term_mult(f[i].coeffs[0],
                                 f[j].coeffs[0])[1:]:
            # TODO: Don't operate on the whole lists, do nested ifs instead?
            kk = filter(lambda k: k!=i and k!=j,range(0, s))
            kk = filter(lambda k: not (min(i,k),max(i,k)) in b, kk)
            kk = filter(lambda k: not (min(j,k),max(j,k)) in b, kk)
            # Check if the lcm is divisible by another base element.
            kk = filter(lambda k: term_is_mult(lcm,f[k].coeffs[0]), kk)
            crit = not bool(kk)
        if crit:
            factor_i = Polynomial(coeffs=(term_div(lcm, f[i].coeffs[0]),),
                                  var=f[0].var, order=f[0].order)
            factor_j = Polynomial(coeffs=(term_div(lcm, f[j].coeffs[0]),),
                                  var=f[0].var, order=f[0].order)
            s_poly = f[i]*factor_i - f[j]*factor_j
            s_poly = div_.div(s_poly, f)[-1] # reduce
            if s_poly.sympy_expr is not S.Zero:
                # we still have to add it to the base.
                s += 1
                f.append(s_poly)
                for t in range(0, s-1): # With a new element come
                    b.append((t, s-1))  # new combinationas to test.
        b = b[1:] # Checked one more.

    # We now have one possible Groebner base, probably too big.
    if not reduced:
        return f

    # We can get rid of all elements, where the leading term can be
    # reduced in the ideal of the remaining leading terms, that is,
    # can be divided by one of the other leading terms.
    blacklist = []
    for p in f:
        if filter(lambda x: term_is_mult(p.coeffs[0], x.coeffs[0]),
               filter(lambda x: not x in blacklist and x != p, f)):
            blacklist.append(p)
    for p in blacklist:
        f.remove(p)

    # We can now sort the basis elements according to their leading
    # term.
    f.sort(cmp=lambda a,b: term_cmp(a.coeffs[0], b.coeffs[0],
                                    a.order), reverse=True)

    # Divide all basis elements by their leading coefficient, to get a
    # leading 1.
    for i, p in enumerate(f):
        c = p.coeffs[0][0]
        f[i] = Polynomial(coeffs=tuple(map(lambda t:(t[0]/c,) + t[1:],
                                           p.coeffs)),
                          var=p.var, order=p.order)

    # We now have a minimal Groebner basis, which is still not unique.
    # The next step is to reduce all basis elements in respect to the
    # rest of the base (without touching the leading terms).
    # As the basis is already sorted, the rest gets smaller each time.
    for i,p in enumerate(f[0:-1]):
        pp = div_.div(p, f[i+1:])[-1]
        f[i] = pp

    return f

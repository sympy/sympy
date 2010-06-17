"""
Algorithms for solving the Risch differential equation.

Given a differential field K of characteristic 0 that is a simple monomial
extension of a base field k and f, g in K, the Risch Differential Equation
problem is to decide if there exist y in K such that Dy + f*y == g and to find
one if there are some.  If t is a monomial over k and the coefficients of f and
g are in k(t), then y is in k(t), and the out line of the algorithm here is
given as:

1. Compute the normal part n of the denominator of y.  The problem is then
reduced to finding y' in k<t>, where y == y'/n.
2. Compute the special part s of the denominator of y.   The problem is then
reduced to finding y'' in k[t], where y == y''/(n*s)
3. Bound the degree of y''.
4. Reduce the equation Dy + f*y == g to a similar equation with f, g in k[t].
5. Find the solutions in k[t] of bounded degree of the reduced equation.

See Chapter 6 of "Symbolic Integration I: Transcendental Functions" by Manuel
Bronstein.
"""
from sympy.core.symbol import Symbol

from sympy.polys import Poly, gcd, ZZ

from sympy.integrals.risch import (gcdex_diophantine, derivation, splitfactor)

from operator import mul
#    from pudb import set_trace; set_trace() # Debugging

def weak_normalizer(a, d, D, x, t, z=None):
    """
    Weak normalization.

    Given a derivation D on k[t] and f == a/d in k(t), return q in k[t] such
    that f - Dq/q is weakly normalized with respect to t.

    f in k(t) is said to be "weakly normalized" with respect to t if
    residue_p(f) is not a positive integer for any normal irreducible p in k[t]
    such that f is in R_p (Definition 6.1.1).  If f has an elementary integral,
    this is equivalent to no logarithm of integral(f) whose argument depends on
    t has a positive integer coefficient, where the arguments of the logarithms
    not in k(t) are in k[t].

    Returns (q, f - Dq/q)
    """
    z = z or Symbol('z', dummy=True)
    dn, ds = splitfactor(d, D, x, t)

    # Compute d1, where dn == d1*d2**2*...*dn**n is a square-free
    # factorization of d.
    g = gcd(dn, dn.diff(t))
    d_sqf_part = dn.quo(g)
    d1 = d_sqf_part.quo(gcd(d_sqf_part, g))

    a1, b = gcdex_diophantine(d.quo(d1), d1, a)
    r = (a - Poly(z, t)*derivation(d1, D, x, t)).as_poly(t).resultant(d1.as_poly(t))
    r = Poly(r, z)

    if not r.has(z):
        return (Poly(1, t), (a, d))

    N = [i for i in r.real_roots() if i in ZZ and i > 0]

    q = reduce(mul, [gcd(a - Poly(n, t)*derivation(d1, D, x, t), d1) for n in N],
        Poly(1, t))

    dq = derivation(q, D, x, t)
    sn = q*a - d*dq
    sd = q*d
    sn, sd = sn.cancel(sd, include=True)

    return (q, (sn, sd))

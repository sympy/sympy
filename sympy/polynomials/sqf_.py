"""Algorithms for square-free decomposition of polynomials"""

from sympy.polynomials.base import *
from sympy.polynomials import div_

def uv(f):
    """Returns a decomposition of f in a1 * a2**2 * ... * an**n.

    Here, the ai are pairwise prime and square-free polynomials, returned
    in a list. f is assumed to be a univariate instance of Polynomial.

    """
    f = [f]
    while f[-1].coeffs[0][1] is not S.Zero:
        f.append(div_.gcd(f[-1], f[-1].diff(f[-1].var[0])))
    g = []
    for i in range(1, len(f)):
        g.append(div_.div(f[i-1], f[i])[0])
    a = []
    for i in range(0, len(g)-1):
        a.append(div_.div(g[i], g[i+1])[0])
    a.append(g[-1])
    return a

def uv_int(f):
    """Returns a decomposition of f in a1 * a2**2 * ... * an**n.

    Here, the ai are pairwise prime and square-free polynomials, returned
    in a list. f is assumed to be a univariate instance of Polynomial.
    The numeric factor is re-distributed to keep integer coeffs.
    """
    a = uv(f)

    ca = int(a[0].content())
    c = ca
    for i,p in enumerate(a[1:]):
        # Compute lcm of denominators in coeffs:
        l = 1
        for term in p.coeffs:
            l = l*term[0].q / numbers.gcd(l ,term[0].q)
        assert c % (l**(i+2)) == 0
        c /= l**(i+2)
        p = Polynomial(coeffs=tuple(map(lambda t:(t[0]*Rational(l),) + t[1:],
                                        p.coeffs)),
                       var=p.var, order=p.order)
    ca /= c
    a[0] = Polynomial(coeffs=tuple(map(lambda t:(t[0]*Rational(1,ca),) + t[1:],
                                       a[0].coeffs)),
                      var=a[0].var, order=a[0].order)
    return a

def uv_part(f):
    """Returns the square-free part of f.

    f is assumed to be a univariate instance of Polynomial.

    """
    ff = div_.gcd(f, f.diff(f.var[0]))
    return div_.div(f, ff)[0]

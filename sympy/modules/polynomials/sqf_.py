"""Algorithms for square-free decomposition of polynomials"""

from sympy.modules.polynomials.base import *
from sympy.modules.polynomials import gcd_
from sympy.modules.polynomials import div_

def uv(f):
    """Returns a decomposition of f in a1 * a2**2 * ... * an**n.

    Here, the ai are pairwise prime and square-free polynomials, returned
    in a list. f is assumed to be a univariate instance of Polynomial.

    """
    f = [f]
    while f[-1].cl[0][1] != 0:
        f.append(gcd_.uv(f[-1], f[-1].diff(f[-1].var[0])))
    g = []
    for i in range(1, len(f)):
        g.append(div_.mv(f[i-1], f[i])[0][0])
    a = []
    for i in range(0, len(g)-1):
        a.append(div_.mv(g[i], g[i+1])[0][0])
    a.append(g[-1])
    return a

def uv_int(f):
    """Returns a decomposition of f in a1 * a2**2 * ... * an**n.

    Here, the ai are pairwise prime and square-free polynomials, returned
    in a list. f is assumed to be a univariate instance of Polynomial.
    The numeric factor is re-distributed to keep integer coefficients.
    """
    assert f.coeff == 'int'
    a = uv(f)

    ca = int(a[0].content())
    c = ca
    for i,p in enumerate(a[1:]):
        # Compute lcm of denominators in coefficients:
        l = 1
        for term in p.cl:
            l = l*term[0].q / numbers.gcd(l ,term[0].q)
        assert c % (l**(i+2)) == 0
        c /= l**(i+2)
        p.cl = map(lambda t:[t[0]*Rational(l)] + t[1:], p.cl)
    ca /= c
    a[0].cl = map(lambda t:[t[0]*Rational(1,ca)] + t[1:], a[0].cl)
    return a

def uv_part(f):
    """Returns the square-free part of f.

    f is assumed to be a univariate instance of Polynomial.

    """
    ff = gcd_.uv(f, f.diff(f.var[0]))
    return div_.mv(f, ff)[0][0]

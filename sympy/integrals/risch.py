from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.symbol import Symbol, Wild
from sympy.core.basic import S, C, sympify
from sympy.core.numbers import Rational

from sympy.functions import exp, sin , cos , tan , cot , asin
from sympy.functions import log, sinh, cosh, tanh, coth, asinh

from sympy.functions import sqrt, erf

from sympy.solvers import solve

from sympy.polys import quo, gcd, lcm, \
    monomials, factor, cancel, PolynomialError, Poly
from sympy.polys.polyroots import root_factors

from sympy.utilities.iterables import make_list


#Zero = Poly(0, *V)

#One = Poly(1, *V)

def derivation(p, D, x, t):
    """
    Computes the Dp, given the derivation D with D = d/dx and p is a polynomial
    in t over K(x)
    """
    px = p.as_poly(t, x)
    if px is None:
        px = p.as_basic()
    return p.diff(t)*D + px.diff(x).as_poly(t)

def splitfactor(p, D, x, t):
    """
    Splitting factorization.

    Given a derivation D on k[t] and p in k[t], return (p_n, p_s) in k[t] x k[t]
    such that p = p_n*p_s, p_s is special, and each square factor of p_n is
    normal.

    Page. 100
    """
    Zero = Poly(0, t)
    One = Poly(1, t)
    if not p.has_any_symbols(t):
        return (p, One)

    Dp = derivation(p, D, x, t)
    if Dp != Zero:
        h = p.gcd(Dp)
        g = p.gcd(p.diff(t))
        s = h.quo(g)

        if s.degree(t) == 0:
            return (p, One)

        q_split = splitfactor(p.quo(s), D, x, t)

        return (q_split[0], q_split[1]*s)
    else:
        return (p, One)

def canonical_representation(a, d, D, x, t):
    """
    Canonical Representation.

    Given a derivation D on k[t] and f = a/d in k(t), return (f_p, f_s, f_n) in
    k[t] x k(t) x k(t) such that f = f_p + f_s + f_n is the canonical
    representation of f (f_p is a polynomial, f_s is reduced (has a special
    denominator), and f_n is simple (has a normal denominator).
    """
    Zero = Poly(0, t)
    # Make d monic
    l = Poly(1/d.LC(), t)
    a, d = a.mul(l), d.mul(l)

    q, r = a.div(d)
    dn, ds = splitfactor(d, D, x, t)

    # Extended Euclidean Algorithm (Diophantine Version) pg. 13
    # XXX: Should this go in densetools.py ?
    b, c, g = dn.gcdex(ds)
    q1 = r.quo(g)
    b, c = q1*b, q1*c
    if b != Zero and b.degree() >= b.degree():
        q1, r1 = b.div(ds)
        b, r1 = r1, c + q1*dn

    return (q, (b,ds), (c,dn))

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

    if derivation(p, D, x, t) != Zero:
        c, q = p.as_poly(t).primitive()
        c, q = Poly(c, t), Poly(q, t)

        h = q.gcd(derivation(q, D, x, t))
        s = h.quo(q.gcd(q.diff(t)))

        c_split = splitfactor(c, D, x, t)

        if s.degree(t) == 0:
            return (c_split[0], q * c_split[1])

        q_split = splitfactor(q.quo(s), D, x, t)

        return (c_split[1]*q_split[1], c_split[0]*q_split[0]*s)
    else:
        return (p, One)


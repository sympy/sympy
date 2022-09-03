"""Efficient functions for generating Appell sequences."""


from sympy.core.symbol import Dummy
from sympy.polys.densearith import dup_mul_ground, dup_add_ground
from sympy.polys.densetools import dup_eval, dup_integrate
from sympy.polys.domains import ZZ, QQ
from sympy.polys.polyclasses import DMP
from sympy.polys.polytools import Poly, PurePoly
from sympy.utilities import public

def dup_appell(n, seq, v, f, K):
    """Common low-level implementation of Appell sequences."""
    for i in range(2, n+1):
        p = dup_integrate(dup_mul_ground(seq[-1], i, K), 1, K)
        seq.append(dup_add_ground(p, f(dup_eval(p, v, K), i), K))
    return seq[n]

def appell_poly(n, seq, v, f, K, x=None, polys=False):
    """Generates the nth polynomial in `x` of the Appell sequence with
    parameters `seq`, `v` and `f`.

    Parameters
    ==========

    n : int
        Order of the polynomial.
    seq : iterable
    v : Expr
    f : callable
    K : Domain
        Domain in which to perform computations and in which the coefficients
        of the specified sequence's polynomials lie in.
    x : optional
    polys : bool, optional
        If True, return a Poly, otherwise (default) return an expression.
    """
    if n < 0:
        raise ValueError(
            "Cannot generate Appell sequence polynomial of order %s" % n)
    poly = DMP(dup_appell(int(n), seq, v, f, K), K)
    if x is not None:
        poly = Poly.new(poly, x)
    else:
        poly = PurePoly.new(poly, Dummy('x'))
    return poly if polys else poly.as_expr()


@public
def bernoulli_poly(n, x=None, polys=False):
    """Generates the Bernoulli polynomial of degree `n` in `x`.

    Parameters
    ==========

    n : int
        Degree of the polynomial.
    x : optional
    polys : bool, optional
        If True, return a Poly, otherwise (default) return an expression.
    """
    return appell_poly(n, [[1], [1, QQ(-1,2)]], QQ(1,2),
            lambda p, i: p * QQ(1<<(i-1), 1-(1<<i)), QQ, x, polys)


@public
def genocchi_poly(n, x=None, polys=False):
    """Generates the Genocchi polynomial of degree `n-1` in `x`.

    Parameters
    ==========

    n : int
        Degree of the polynomial plus one.
    x : optional
    polys : bool, optional
        If True, return a Poly, otherwise (default) return an expression.
    """
    return appell_poly(n, [[0], [-1]], 1, lambda p, i: -(p>>1), ZZ, x, polys)


@public
def euler_poly(n, x=None, polys=False):
    """Generates the Euler polynomial of degree `n` in `x`.

    Parameters
    ==========

    n : int
        Degree of the polynomial.
    x : optional
    polys : bool, optional
        If True, return a Poly, otherwise (default) return an expression.
    """
    return appell_poly(n, [[1], [1, QQ(-1,2)]], 1, lambda p, i: -p / 2, QQ, x, polys)


@public
def andre_poly(n, x=None, polys=False):
    """Generates the André polynomial of degree `n` in `x`.

    Parameters
    ==========

    n : int
        Degree of the polynomial.
    x : optional
    polys : bool, optional
        If True, return a Poly, otherwise (default) return an expression.
    """
    return appell_poly(n, [[1], [1, 0]], 1, lambda p, i: p * ((i&1)-1), ZZ, x, polys)

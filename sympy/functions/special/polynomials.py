"""
This module mainly implements special orthogonal polynomials.

See also functions.combinatorial.numbers which contains some
combinatorial polynomials.

"""

from sympy.core.basic import C
from sympy.core.singleton import S
from sympy.core import Rational
from sympy.core.function import Function
from sympy.functions.combinatorial.factorials import factorial

from sympy.polys.orthopolys import (
    chebyshevt_poly,
    chebyshevu_poly,
    laguerre_poly,
    hermite_poly,
    legendre_poly,
)

_x = C.Dummy('x')

class PolynomialSequence(Function):
    """Polynomial sequence with one index and n >= 0. """

    nargs = 2

    @classmethod
    def eval(cls, n, x):
        if n.is_integer and n >= 0:
            return cls._ortho_poly(int(n), _x).subs(_x, x)
        if n.is_negative:
            raise ValueError("%s index must be nonnegative integer (got %r)" % (cls, n))

#----------------------------------------------------------------------------
# Chebyshev polynomials of first and second kind
#

class chebyshevt(PolynomialSequence):
    """
    chebyshevt(n, x) gives the nth Chebyshev polynomial (of the first
    kind) of x, T_n(x)

    The Chebyshev polynomials of the first kind are orthogonal on
    [-1, 1] with respect to the weight 1/sqrt(1-x**2).

    Examples
    ========
        >>> from sympy import chebyshevt
        >>> from sympy.abc import x
        >>> chebyshevt(0, x)
        1
        >>> chebyshevt(1, x)
        x
        >>> chebyshevt(2, x)
        -1 + 2*x**2

    References
    ==========
    * http://en.wikipedia.org/wiki/Chebyshev_polynomial
    """

    """
    Chebyshev polynomial of the first kind, T_n(x)
    """

    _ortho_poly = staticmethod(chebyshevt_poly)

class chebyshevu(PolynomialSequence):
    """
    chebyshevu(n, x) gives the nth Chebyshev polynomial of the second
    kind of x, U_n(x)

    The Chebyshev polynomials of the second kind are orthogonal on
    [-1, 1] with respect to the weight sqrt(1-x**2).

    Examples
    ========
        >>> from sympy import chebyshevu
        >>> from sympy.abc import x
        >>> chebyshevu(0, x)
        1
        >>> chebyshevu(1, x)
        2*x
        >>> chebyshevu(2, x)
        -1 + 4*x**2

    """

    _ortho_poly = staticmethod(chebyshevu_poly)

class chebyshevt_root(Function):
    """
    chebyshev_root(n, k) returns the kth root (indexed from zero) of
    the nth Chebyshev polynomial of the first kind; that is, if
    0 <= k < n, chebyshevt(n, chebyshevt_root(n, k)) == 0.

    Examples
    ========

    >>> from sympy import chebyshevt, chebyshevt_root
    >>> chebyshevt_root(3, 2)
    -3**(1/2)/2
    >>> chebyshevt(3, chebyshevt_root(3, 2))
    0

    """

    nargs = 2

    @classmethod
    def eval(cls, n, k):
        if not 0 <= k < n:
            raise ValueError("must have 0 <= k < n")
        return C.cos(S.Pi*(2*k+1)/(2*n))

class chebyshevu_root(Function):
    """
    chebyshevu_root(n, k) returns the kth root (indexed from zero) of the
    nth Chebyshev polynomial of the second kind; that is, if 0 <= k < n,
    chebyshevu(n, chebyshevu_root(n, k)) == 0.

    Examples
    ========

        >>> from sympy import chebyshevu, chebyshevu_root
        >>> chebyshevu_root(3, 2)
        -2**(1/2)/2
        >>> chebyshevu(3, chebyshevu_root(3, 2))
        0

    """

    nargs = 2

    @classmethod
    def eval(cls, n, k):
        if not 0 <= k < n:
            raise ValueError("must have 0 <= k < n")
        return C.cos(S.Pi*(k+1)/(n+1))

#----------------------------------------------------------------------------
# Legendre polynomials  and  Associated Legendre polynomials
#

class legendre(PolynomialSequence):
    """
    legendre(n, x) gives the nth Legendre polynomial of x, P_n(x)

    The Legendre polynomials are orthogonal on [-1, 1] with respect to
    the constant weight 1. They satisfy P_n(1) = 1 for all n; further,
    P_n is odd for odd n and even for even n

    Examples
    ========
        >>> from sympy import legendre
        >>> from sympy.abc import x
        >>> legendre(0, x)
        1
        >>> legendre(1, x)
        x
        >>> legendre(2, x)
        -1/2 + 3*x**2/2

    References
    ==========
    * http://en.wikipedia.org/wiki/Legendre_polynomial
    """

    _ortho_poly = staticmethod(legendre_poly)

class assoc_legendre(Function):
    """
    assoc_legendre(n,m, x) gives P_nm(x), where n and m are the degree
    and order or an expression which is related to the nth order
    Legendre polynomial, P_n(x) in the following manner:

        P_nm(x) = (-1)**m * (1 - x**2)**(m/2) * diff(P_n(x), x, m)

    Associated Legendre polynomial are orthogonal on [-1, 1] with:

    - weight = 1            for the same m, and different n.
    - weight = 1/(1-x**2)   for the same n, and different m.

    Examples
    ========
        >>> from sympy import assoc_legendre
        >>> from sympy.abc import x
        >>> assoc_legendre(0,0, x)
        1
        >>> assoc_legendre(1,0, x)
        x
        >>> assoc_legendre(1,1, x)
        -(1 - x**2)**(1/2)

    References
    ==========
    * http://en.wikipedia.org/wiki/Associated_Legendre_polynomials
    """

    nargs = 3

    @classmethod
    def calc(cls, n, m):
        P = legendre_poly(n, _x, polys=True).diff((_x, m))
        return (-1)**m * (1 - _x**2)**Rational(m, 2) * P.as_expr()

    @classmethod
    def eval(cls, n, m, x):
        if n.is_integer and n >= 0 and m.is_integer and abs(m) <= n:
            assoc = cls.calc(int(n), abs(int(m)))

            if m < 0:
                assoc *= (-1)**(-m) * (C.Factorial(n + m)/C.Factorial(n - m))

            return assoc.subs(_x, x)

        if n.is_negative:
            raise ValueError("%s : 1st index must be nonnegative integer (got %r)" % (cls, n))

        if abs(m) > n:
            raise ValueError("%s : abs('2nd index') must be <= '1st index' (got %r, %r)" % (cls, n, m))

#----------------------------------------------------------------------------
# Hermite polynomials
#

class hermite(PolynomialSequence):
    """
    hermite(n, x) gives the nth Hermite polynomial in x, H_n(x)

    The Hermite polynomials are orthogonal on (-oo, oo) with respect to
    the weight exp(-x**2/2).

    Examples
    ========
        >>> from sympy import hermite
        >>> from sympy.abc import x
        >>> hermite(0, x)
        1
        >>> hermite(1, x)
        2*x
        >>> hermite(2, x)
        -2 + 4*x**2

    References
    ==========
    * http://mathworld.wolfram.com/HermitePolynomial.html
    """

    _ortho_poly = staticmethod(hermite_poly)

def laguerre_l(n, alpha, x):
    """
    Returns the generalized Laguerre polynomial.

    n     ... 0, 1, 2, 3, ...
    alpha ... any symbol (alpha=0 gives regular Laguerre polynomials)

    Examples::

    >>> from sympy import laguerre_l, var
    >>> var("alpha, x")
    (alpha, x)
    >>> laguerre_l(0, alpha, x)
    1
    >>> laguerre_l(1, alpha, x)
    1 + alpha - x
    >>> laguerre_l(2, alpha, x)
    (1 + alpha)*(2 + alpha)/2 - x*(2 + alpha) + x**2/2

    If you set alpha=0, you get regular Laguerre polynomials::

    >>> laguerre_l(1, 0, x)
    1 - x
    >>> laguerre_l(2, 0, x)
    1 - 2*x + x**2/2
    >>> laguerre_l(3, 0, x)
    1 - 3*x + 3*x**2/2 - x**3/6
    >>> laguerre_l(4, 0, x)
    1 - 4*x + 3*x**2 - 2*x**3/3 + x**4/24

    """
    n, alpha, x = S(n), S(alpha), S(x)
    r = 0
    for m in range(n+1):
        c = 1
        for i in range(m+1, n+1):
            c *= alpha+i
        r += (-1)**m * c * x**m/(factorial(m)*factorial(n-m))
    return r

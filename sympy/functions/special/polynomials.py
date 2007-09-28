"""
This module mainly implements special orthogonal polynomials.

See also functions.combinatorial.numbers which contains some
combinatorial polynomials.

"""

from sympy.core.basic import Basic, S
from sympy.core import Rational, Symbol
from sympy.core.function import SingleValuedFunction
from sympy.utilities.memoization import recurrence_memo


_x = Basic.Symbol('x', dummy=True)

class PolynomialSequence(SingleValuedFunction):

    nofargs = 2
    precedence = Basic.Apply_precedence

    @classmethod
    def _eval_apply(cls, n, x):
        if n.is_integer and n >= 0:
            return cls.calc(int(n)).subs(_x, x)
        if n.is_negative:
            raise ValueError("%s index must be nonnegative integer (got %r)" % (self, i))



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
        >>> x = Symbol('x')
        >>> chebyshevt(0, x)
        1
        >>> chebyshevt(1, x)
        x
        >>> chebyshevt(2, x)
        (-1) + 2*x**2

    References
    ==========
    * http://en.wikipedia.org/wiki/Chebyshev_polynomial
    """

    """
    Chebyshev polynomial of the first kind, T_n(x)
    """
    @staticmethod
    @recurrence_memo([Basic.One(), _x])
    def calc(n, prev):
        return (2*_x*prev[n-1] - prev[n-2]).expand()


class chebyshevu(PolynomialSequence):
    """
    chebyshevu(n, x) gives the nth Chebyshev polynomial of the second
    kind of x, U_n(x)

    The Chebyshev polynomials of the second kind are orthogonal on
    [-1, 1] with respect to the weight sqrt(1-x**2).

    Examples
    ========
        >>> x = Symbol('x')
        >>> chebyshevu(0, x)
        1
        >>> chebyshevu(1, x)
        2*x
        >>> chebyshevu(2, x)
        (-1) + 4*x**2
    """
    @staticmethod
    @recurrence_memo([Basic.One(), 2*_x])
    def calc(n, prev):
        return (2*_x*prev[n-1] - prev[n-2]).expand()


class chebyshevt_root(SingleValuedFunction):
    """
    chebyshev_root(n, k) returns the kth root (indexed from zero) of
    the nth Chebyshev polynomial of the first kind; that is, if
    0 <= k < n, chebyshevt(n, chebyshevt_root(n, k)) == 0.

    Examples
    ========

    >>> chebyshevt_root(3, 2)
    -1/2*3**(1/2)
    >>> chebyshevt(3, chebyshevt_root(3, 2))
    0
    """
    nofargs = 2

    @classmethod
    def _eval_apply(cls, n, k):
        if not 0 <= k < n:
            raise ValueError, "must have 0 <= k < n"
        return Basic.cos(S.Pi*(2*k+1)/(2*n))

class chebyshevu_root(SingleValuedFunction):
    """
    chebyshevu_root(n, k) returns the kth root (indexed from zero) of the
    nth Chebyshev polynomial of the second kind; that is, if 0 <= k < n,
    chebyshevu(n, chebyshevu_root(n, k)) == 0.

    Examples
    ========

        >>> chebyshevu_root(3, 2)
        -1/2*2**(1/2)
        >>> chebyshevu(3, chebyshevu_root(3, 2))
        0
    """
    nofargs = 2

    @classmethod
    def _eval_apply(cls, n, k):
        if not 0 <= k < n:
            raise ValueError, "must have 0 <= k < n"
        return Basic.cos(S.Pi*(k+1)/(n+1))


#----------------------------------------------------------------------------
# Legendre polynomials
#

class legendre(PolynomialSequence):
    """
    legendre(n, x) gives the nth Legendre polynomial of x, P_n(x)

    The Legendre polynomials are orthogonal on [-1, 1] with respect to
    the constant weight 1. They satisfy P_n(1) = 1 for all n; further,
    P_n is odd for odd n and even for even n

    Examples
    ========
        >>> x = Symbol('x')
        >>> legendre(0, x)
        1
        >>> legendre(1, x)
        x
        >>> legendre(2, x)
        (-1/2) + (3/2)*x**2

    References
    ========
    * http://en.wikipedia.org/wiki/Legendre_polynomial
    """
    @staticmethod
    @recurrence_memo([Basic.One(), _x])
    def calc(n, prev):
        return (((2*n-1)*_x*prev[n-1] - (n-1)*prev[n-2])/n).expand()


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
        >>> x = Symbol('x')
        >>> hermite(0, x)
        1
        >>> hermite(1, x)
        2*x
        >>> hermite(2, x)
        (-2) + 4*x**2

    References
    ==========
    * http://mathworld.wolfram.com/HermitePolynomial.html    
    """
    @staticmethod
    @recurrence_memo([Basic.One(), 2*_x])
    def calc(n, prev):
        return (2*_x*prev[n-1]-2*(n-1)*prev[n-2]).expand()

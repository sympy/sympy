"""
This module mainly implements special orthogonal polynomials.

See also functions.combinatorial.numbers which contains some
combinatorial polynomials.

"""

from sympy.core.basic import S, C
from sympy.core import Rational
from sympy.core.function import Function
from sympy.utilities.memoization import recurrence_memo, assoc_recurrence_memo

_x = C.Symbol('x', dummy=True)

class PolynomialSequence(Function):
    """Polynomial sequence with 1 index

       n >= 0
    """

    nargs = 2

    @classmethod
    def eval(cls, n, x):
        if n.is_integer and n >= 0:
            return cls.calc(int(n)).subs(_x, x)
        if n.is_negative:
            raise ValueError("%s index must be nonnegative integer (got %r)" % (cls, n))


class PolynomialSequence2(Function):
    """Polynomial sequence with 2 indexes

       n >= 0
       abs(m) <= n
    """

    nargs = 3

    @classmethod
    def eval(cls, n, m, x):
        if n.is_integer and n >= 0 and m.is_integer and abs(m) <= n:
            return cls.calc2(int(n), int(m)).subs(_x, x)

        if n.is_negative:
            raise ValueError("%s : 1st index must be nonnegative integer (got %r)" % (cls, n))

        if abs(m) > n:
            raise ValueError("%s : abs('2nd index') must be <= '1st index' (got %r, %r)" % (cls, n, m))



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
    @staticmethod
    @recurrence_memo([S.One, _x])
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
        >>> from sympy import chebyshevu
        >>> from sympy.abc import x
        >>> chebyshevu(0, x)
        1
        >>> chebyshevu(1, x)
        2*x
        >>> chebyshevu(2, x)
        -1 + 4*x**2

    """
    @staticmethod
    @recurrence_memo([S.One, 2*_x])
    def calc(n, prev):
        return (2*_x*prev[n-1] - prev[n-2]).expand()


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
    @staticmethod
    @recurrence_memo([S.One, _x])
    def calc(n, prev):
        return (((2*n-1)*_x*prev[n-1] - (n-1)*prev[n-2])/n).expand()


class assoc_legendre(PolynomialSequence2):
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

    @staticmethod
    @assoc_recurrence_memo(legendre.calc)
    def _calc2(n, m, prev):
        P=prev

        # this is explicit (not recurrence) formula
        #
        # the result is pretty and we still benefit from memoization
        Pnm = (-1)**m * (1-_x**2)**Rational(m,2) * P[n][0].diff(_x, m)
        return Pnm


        # this is reccurence formula, but the cost to keep it pretty (simplify)
        # is too high
        Pnm = (n-(m-1))*_x*P[n][m-1] - (n+(m-1))*P[n-1][m-1]

        # hack to simplify the expression
        # FIXME something more lightweight?
        from sympy.simplify import simplify

        Pnm = Pnm / (1-_x**2)**Rational(m+1,2)
        Pnm = simplify(Pnm)

        Pnm*= (1-_x**2)**Rational(m+1-1,2)

        return Pnm

    @staticmethod
    def calc2(n,m):
        if m >= 0:
            return assoc_legendre._calc2(n,m)
        else:
            factorial = C.Factorial
            m = -m
            return (-1)**m *factorial(n-m)/factorial(n+m) * assoc_legendre._calc2(n, m)



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
    @staticmethod
    @recurrence_memo([S.One, 2*_x])
    def calc(n, prev):
        return (2*_x*prev[n-1]-2*(n-1)*prev[n-2]).expand()

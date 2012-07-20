"""
This module mainly implements special orthogonal polynomials.

See also functions.combinatorial.numbers which contains some
combinatorial polynomials.

"""

from sympy.core.basic import C
from sympy.core.singleton import S
from sympy.core import Rational
from sympy.core.function import Function, ArgumentIndexError
from sympy.functions.elementary.miscellaneous import sqrt

from sympy.polys.orthopolys import (
    chebyshevt_poly,
    chebyshevu_poly,
    laguerre_poly,
    hermite_poly,
    legendre_poly
)

_x = C.Dummy('x')


class OrthogonalPolynomial(Function):
    """Base class for orthogonal polynomials.
    """
    nargs = 2

    @classmethod
    def _eval_at_order(cls, n, x):
        if n.is_integer and n >= 0:
            return cls._ortho_poly(int(n), _x).subs(_x, x)

    def _eval_conjugate(self):
        return self.func(self.args[0], self.args[1].conjugate())


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
    2*x**2 - 1

    See Also
    ========

    chebysevt_root, chebyshevu, chebyshevu_root
    legendre, assoc_legendre

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
    4*x**2 - 1

    See Also
    ========

    chebyshevu_root, chebyshevt, chebyshevt_root
    legendre, assoc_legendre

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
    -sqrt(3)/2
    >>> chebyshevt(3, chebyshevt_root(3, 2))
    0

    See Also
    ========

    chebyshevt, chebyshevu, chebyshevu_root
    legendre, assoc_legendre

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
    -sqrt(2)/2
    >>> chebyshevu(3, chebyshevu_root(3, 2))
    0

    See Also
    ========

    chebyshevu, chebyshevt, chebyshevt_root
    legendre, assoc_legendre

    """

    nargs = 2

    @classmethod
    def eval(cls, n, k):
        if not 0 <= k < n:
            raise ValueError("must have 0 <= k < n")
        return C.cos(S.Pi*(k+1)/(n+1))

#----------------------------------------------------------------------------
# Legendre polynomials and Associated Legendre polynomials
#

class legendre(OrthogonalPolynomial):
    r"""
    legendre(n, x) gives the nth Legendre polynomial of x, :math:`P_n(x)`

    The Legendre polynomials are orthogonal on [-1, 1] with respect to
    the constant weight 1. They satisfy :math:`P_n(1) = 1` for all n; further,
    :math:`P_n` is odd for odd n and even for even n.

    Examples
    ========

    >>> from sympy import legendre, diff
    >>> from sympy.abc import x, n
    >>> legendre(0, x)
    1
    >>> legendre(1, x)
    x
    >>> legendre(2, x)
    3*x**2/2 - 1/2
    >>> legendre(n, x)
    legendre(n, x)
    >>> diff(legendre(n,x), x)
    n*(x*legendre(n, x) - legendre(n - 1, x))/(x**2 - 1)

    See Also
    ========

    assoc_legendre
    chebyshevu, chebyshevt, chebyshevu_root, chebyshevt_root

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Legendre_polynomial
    """

    _ortho_poly = staticmethod(legendre_poly)

    @classmethod
    def eval(cls, n, x):
        if not n.is_Number:
            # Symbolic result L_n(x)
            # L_n(-x)  --->  (-1)**n * L_n(x)
            if x.could_extract_minus_sign():
                return S.NegativeOne**n * legendre(n,-x)
            # L_{-n}(x)  --->  L_{n-1}(x)
            if n.could_extract_minus_sign():
                return legendre(-n - S.One,x)
            # We can evaluate for some special values of x
            if x == S.Zero:
                return sqrt(S.Pi)/(C.gamma(S.Half - n/2)*C.gamma(S.One + n/2))
            elif x == S.One:
                return S.One
            elif x == S.Infinity:
                return S.Infinity
        else:
            # n is a given fixed integer, evaluate into polynomial
            return cls._eval_at_order(n, x)

    def fdiff(self, argindex=2):
        if argindex == 1:
            # Diff wrt n
            raise ArgumentIndexError(self, argindex)
        elif argindex == 2:
            # Diff wrt x
            # Find better formula, this is unsuitable for x = 1
            n, x = self.args
            return n/(x**2 - 1)*(x*legendre(n, x) - legendre(n - 1, x))
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_polynomial(self, n, x):
        k = C.Dummy("k")
        kern = (-1)**k*C.binomial(n, k)**2*((1 + x)/2)**(n - k)*((1 - x)/2)**k
        return C.Sum(kern, (k, 0, n))

class assoc_legendre(Function):
    r"""
    assoc_legendre(n,m, x) gives :math:`P_n^m(x)`, where n and m are
    the degree and order or an expression which is related to the nth
    order Legendre polynomial, :math:`P_n(x)` in the following manner:

    .. math::
        P_n^m(x) = (-1)^m (1 - x^2)^{\frac{m}{2}}
                   \frac{\mathrm{d}^m P_n(x)}{\mathrm{d} x^m}

    Associated Legendre polynomial are orthogonal on [-1, 1] with:

    - weight = 1            for the same m, and different n.
    - weight = 1/(1-x**2)   for the same n, and different m.

    Examples
    ========

    >>> from sympy import assoc_legendre
    >>> from sympy.abc import x, m, n
    >>> assoc_legendre(0,0, x)
    1
    >>> assoc_legendre(1,0, x)
    x
    >>> assoc_legendre(1,1, x)
    -sqrt(-x**2 + 1)
    >>> assoc_legendre(n,m,x)
    assoc_legendre(n, m, x)

    See Also
    ========

    legendre
    chebyshevu, chebyshevt, chebyshevu_root, chebyshevt_root

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Associated_Legendre_polynomials
    """

    nargs = 3

    @classmethod
    def _eval_at_order(cls, n, m):
        P = legendre_poly(n, _x, polys=True).diff((_x, m))
        return (-1)**m * (1 - _x**2)**Rational(m, 2) * P.as_expr()

    @classmethod
    def eval(cls, n, m, x):
        if m.could_extract_minus_sign():
            # P^{-m}_n  --->  F * P^m_n
            return S.NegativeOne**(-m) * (C.factorial(m + n)/C.factorial(n - m)) * assoc_legendre(n, -m, x)
        if m == 0:
            # P^0_n  --->  L_n
            return legendre(n, x)
        if x == 0:
            return 2**m*sqrt(S.Pi) / (C.gamma((1 - m - n)/2)*C.gamma(1 - (m - n)/2))
        if n.is_Number and m.is_Number and n.is_integer and m.is_integer:
            if n.is_negative:
                raise ValueError("%s : 1st index must be nonnegative integer (got %r)" % (cls, n))
            if abs(m) > n:
                raise ValueError("%s : abs('2nd index') must be <= '1st index' (got %r, %r)" % (cls, n, m))
            return cls._eval_at_order(int(n), abs(int(m))).subs(_x, x)

    def fdiff(self, argindex=3):
        if argindex == 1:
            # Diff wrt n
            raise ArgumentIndexError(self, argindex)
        elif argindex == 2:
            # Diff wrt m
            raise ArgumentIndexError(self, argindex)
        elif argindex == 3:
            # Diff wrt x
            # Find better formula, this is unsuitable for x = 1
            n, m, x = self.args
            return 1/(x**2 - 1)*(x*n*assoc_legendre(n,m,x) - (m + n)*assoc_legendre(n - 1,m,x))
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_polynomial(self, n, m, x):
        k = C.Dummy("k")
        kern = C.factorial(2*n - 2*k)/(2**n*C.factorial(n - k)*C.factorial(k)*C.factorial(n - 2*k - m))*(-1)**k*x**(n - m - 2*k)
        return (1 - x**2)**(m/2) * C.Sum(kern, (k, 0, C.floor((n - m)*S.Half)))

#----------------------------------------------------------------------------
# Hermite polynomials
#

class hermite(PolynomialSequence):
    """
    hermite(n, x) gives the nth Hermite polynomial in x, H_n(x)

    The Hermite polynomials are orthogonal on (-oo, oo) with respect to
    the weight `exp(-x**2/2)`.

    Examples
    ========

    >>> from sympy import hermite
    >>> from sympy.abc import x
    >>> hermite(0, x)
    1
    >>> hermite(1, x)
    2*x
    >>> hermite(2, x)
    4*x**2 - 2

    See Also
    ========

    sympy.polys.orthopolys.hermite_poly

    References
    ==========

    * http://mathworld.wolfram.com/HermitePolynomial.html

    """

    _ortho_poly = staticmethod(hermite_poly)

def laguerre_l(n, alpha, x):
    """
    Returns the generalized Laguerre polynomial.

    Parameters
    ==========

    n : int
        Degree of Laguerre polynomial. Must be ``n >= 0``.

    alpha : Expr
        Arbitrary expression. For ``alpha=0`` regular Laguerre
        polynomials will be generated.

    Examples
    ========

    To construct generalized Laguerre polynomials issue::

        >>> from sympy import laguerre_l, var
        >>> var("alpha, x")
        (alpha, x)

        >>> laguerre_l(0, alpha, x)
        1
        >>> laguerre_l(1, alpha, x)
        alpha - x + 1
        >>> laguerre_l(2, alpha, x)
        alpha**2/2 + 3*alpha/2 + x**2/2 + x*(-alpha - 2) + 1

    If you set ``alpha=0``, you get regular Laguerre polynomials::

        >>> laguerre_l(1, 0, x)
        -x + 1
        >>> laguerre_l(2, 0, x)
        x**2/2 - 2*x + 1
        >>> laguerre_l(3, 0, x)
        -x**3/6 + 3*x**2/2 - 3*x + 1
        >>> laguerre_l(4, 0, x)
        x**4/24 - 2*x**3/3 + 3*x**2 - 4*x + 1

    See Also
    ========

    sympy.polys.orthopolys.laguerre_poly

    """
    return laguerre_poly(n, x, alpha)

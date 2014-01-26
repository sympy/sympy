"""
This module mainly implements special orthogonal polynomials.

See also functions.combinatorial.numbers which contains some
combinatorial polynomials.

"""

from __future__ import print_function, division

from sympy.core.basic import C
from sympy.core.singleton import S
from sympy.core import Rational
from sympy.core.numbers import I
from sympy.core.function import Function, ArgumentIndexError
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.exponential import exp
from sympy.functions.special.gamma_functions import gamma
from sympy.functions.special.hyper import hyper
from sympy.functions.combinatorial.factorials import factorial, RisingFactorial



from sympy.polys.orthopolys import (
    jacobi_poly,
    gegenbauer_poly,
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

    @classmethod
    def _eval_at_order(cls, n, x):
        if n.is_integer and n >= 0:
            return cls._ortho_poly(int(n), _x).subs(_x, x)

    def _eval_conjugate(self):
        return self.func(self.args[0], self.args[1].conjugate())

#----------------------------------------------------------------------------
# Jacobi polynomials
#


class jacobi(OrthogonalPolynomial):
    r"""
    Jacobi polynomial :math:`P_n^{\left(\alpha, \beta\right)}(x)`

    jacobi(n, alpha, beta, x) gives the nth Jacobi polynomial
    in x, :math:`P_n^{\left(\alpha, \beta\right)}(x)`.

    The Jacobi polynomials are orthogonal on :math:`[-1, 1]` with respect
    to the weight :math:`\left(1-x\right)^\alpha \left(1+x\right)^\beta`.

    Examples
    ========

    >>> from sympy import jacobi, S, conjugate, diff
    >>> from sympy.abc import n,a,b,x

    >>> jacobi(0, a, b, x)
    1
    >>> jacobi(1, a, b, x)
    a/2 - b/2 + x*(a/2 + b/2 + 1)
    >>> jacobi(2, a, b, x)   # doctest:+SKIP
    (a**2/8 - a*b/4 - a/8 + b**2/8 - b/8 + x**2*(a**2/8 + a*b/4 + 7*a/8 +
    b**2/8 + 7*b/8 + 3/2) + x*(a**2/4 + 3*a/4 - b**2/4 - 3*b/4) - 1/2)

    >>> jacobi(n, a, b, x)
    jacobi(n, a, b, x)

    >>> jacobi(n, a, a, x)
    RisingFactorial(a + 1, n)*gegenbauer(n,
        a + 1/2, x)/RisingFactorial(2*a + 1, n)

    >>> jacobi(n, 0, 0, x)
    legendre(n, x)

    >>> jacobi(n, S(1)/2, S(1)/2, x)
    RisingFactorial(3/2, n)*chebyshevu(n, x)/factorial(n + 1)

    >>> jacobi(n, -S(1)/2, -S(1)/2, x)
    RisingFactorial(1/2, n)*chebyshevt(n, x)/factorial(n)

    >>> jacobi(n, a, b, -x)
    (-1)**n*jacobi(n, b, a, x)

    >>> jacobi(n, a, b, 0)
    2**(-n)*gamma(a + n + 1)*hyper((-b - n, -n), (a + 1,), -1)/(factorial(n)*gamma(a + 1))
    >>> jacobi(n, a, b, 1)
    RisingFactorial(a + 1, n)/factorial(n)

    >>> conjugate(jacobi(n, a, b, x))
    jacobi(n, conjugate(a), conjugate(b), conjugate(x))

    >>> diff(jacobi(n,a,b,x), x)
    (a/2 + b/2 + n/2 + 1/2)*jacobi(n - 1, a + 1, b + 1, x)

    See Also
    ========

    gegenbauer,
    chebyshevt_root, chebyshevu, chebyshevu_root,
    legendre, assoc_legendre,
    hermite,
    laguerre, assoc_laguerre,
    sympy.polys.orthopolys.jacobi_poly,
    sympy.polys.orthopolys.gegenbauer_poly
    sympy.polys.orthopolys.chebyshevt_poly
    sympy.polys.orthopolys.chebyshevu_poly
    sympy.polys.orthopolys.hermite_poly
    sympy.polys.orthopolys.legendre_poly
    sympy.polys.orthopolys.laguerre_poly

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Jacobi_polynomials
    .. [2] http://mathworld.wolfram.com/JacobiPolynomial.html
    .. [3] http://functions.wolfram.com/Polynomials/JacobiP/
    """


    @classmethod
    def eval(cls, n, a, b, x):
        # Simplify to other polynomials
        # P^{a, a}_n(x)
        if a == b:
            if a == -S.Half:
                return C.RisingFactorial(S.Half, n) / C.factorial(n) * chebyshevt(n, x)
            elif a == S.Zero:
                return legendre(n, x)
            elif a == S.Half:
                return C.RisingFactorial(3*S.Half, n) / C.factorial(n + 1) * chebyshevu(n, x)
            else:
                return C.RisingFactorial(a + 1, n) / C.RisingFactorial(2*a + 1, n) * gegenbauer(n, a + S.Half, x)
        elif b == -a:
            # P^{a, -a}_n(x)
            return C.gamma(n + a + 1) / C.gamma(n + 1) * (1 + x)**(a/2) / (1 - x)**(a/2) * assoc_legendre(n, -a, x)
        elif a == -b:
            # P^{-b, b}_n(x)
            return C.gamma(n - b + 1) / C.gamma(n + 1) * (1 - x)**(b/2) / (1 + x)**(b/2) * assoc_legendre(n, b, x)

        if not n.is_Number:
            # Symbolic result P^{a,b}_n(x)
            # P^{a,b}_n(-x)  --->  (-1)**n * P^{b,a}_n(-x)
            if x.could_extract_minus_sign():
                return S.NegativeOne**n * jacobi(n, b, a, -x)
            # We can evaluate for some special values of x
            if x == S.Zero:
                return (2**(-n) * C.gamma(a + n + 1) / (C.gamma(a + 1) * C.factorial(n)) *
                        C.hyper([-b - n, -n], [a + 1], -1))
            if x == S.One:
                return C.RisingFactorial(a + 1, n) / C.factorial(n)
            elif x == S.Infinity:
                if n.is_positive:
                    # TODO: Make sure a+b+2*n \notin Z
                    return C.RisingFactorial(a + b + n + 1, n) * S.Infinity
        else:
            # n is a given fixed integer, evaluate into polynomial
            return jacobi_poly(n, a, b, x)

    def fdiff(self, argindex=4):
        if argindex == 1:
            # Diff wrt n
            raise ArgumentIndexError(self, argindex)
        elif argindex == 2:
            # Diff wrt a
            n, a, b, x = self.args
            k = C.Dummy("k")
            f1 = 1 / (a + b + n + k + 1)
            f2 = ((a + b + 2*k + 1) * C.RisingFactorial(b + k + 1, n - k) /
                  ((n - k) * C.RisingFactorial(a + b + k + 1, n - k)))
            return C.Sum(f1 * (jacobi(n, a, b, x) + f2*jacobi(k, a, b, x)), (k, 0, n - 1))
        elif argindex == 3:
            # Diff wrt b
            n, a, b, x = self.args
            k = C.Dummy("k")
            f1 = 1 / (a + b + n + k + 1)
            f2 = (-1)**(n - k) * ((a + b + 2*k + 1) * C.RisingFactorial(a + k + 1, n - k) /
                  ((n - k) * C.RisingFactorial(a + b + k + 1, n - k)))
            return C.Sum(f1 * (jacobi(n, a, b, x) + f2*jacobi(k, a, b, x)), (k, 0, n - 1))
        elif argindex == 4:
            # Diff wrt x
            n, a, b, x = self.args
            return S.Half * (a + b + n + 1) * jacobi(n - 1, a + 1, b + 1, x)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_polynomial(self, n, a, b, x):
        # TODO: Make sure n \in N
        k = C.Dummy("k")
        kern = (C.RisingFactorial(-n, k) * C.RisingFactorial(a + b + n + 1, k) * C.RisingFactorial(a + k + 1, n - k) /
                C.factorial(k) * ((1 - x)/2)**k)
        return 1 / C.factorial(n) * C.Sum(kern, (k, 0, n))

    def _eval_conjugate(self):
        n, a, b, x = self.args
        return self.func(n, a.conjugate(), b.conjugate(), x.conjugate())


def jacobi_normalized(n, a, b, x):
    r"""
    Jacobi polynomial :math:`P_n^{\left(\alpha, \beta\right)}(x)`

    jacobi_normalized(n, alpha, beta, x) gives the nth Jacobi polynomial
    in x, :math:`P_n^{\left(\alpha, \beta\right)}(x)`.

    The Jacobi polynomials are orthogonal on :math:`[-1, 1]` with respect
    to the weight :math:`\left(1-x\right)^\alpha \left(1+x\right)^\beta`.

    This functions returns the polynomials normilzed:

    .. math::

        \int_{-1}^{1}
          P_m^{\left(\alpha, \beta\right)}(x)
          P_n^{\left(\alpha, \beta\right)}(x)
          (1-x)^{\alpha} (1+x)^{\beta} \mathrm{d}x
        = \delta_{m,n}

    Examples
    ========

    >>> from sympy import jacobi_normalized
    >>> from sympy.abc import n,a,b,x

    >>> jacobi_normalized(n, a, b, x)
    jacobi(n, a, b, x)/sqrt(2**(a + b + 1)*gamma(a + n + 1)*gamma(b + n + 1)/((a + b + 2*n + 1)*factorial(n)*gamma(a + b + n + 1)))

    See Also
    ========

    gegenbauer,
    chebyshevt_root, chebyshevu, chebyshevu_root,
    legendre, assoc_legendre,
    hermite,
    laguerre, assoc_laguerre,
    sympy.polys.orthopolys.jacobi_poly,
    sympy.polys.orthopolys.gegenbauer_poly
    sympy.polys.orthopolys.chebyshevt_poly
    sympy.polys.orthopolys.chebyshevu_poly
    sympy.polys.orthopolys.hermite_poly
    sympy.polys.orthopolys.legendre_poly
    sympy.polys.orthopolys.laguerre_poly

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Jacobi_polynomials
    .. [2] http://mathworld.wolfram.com/JacobiPolynomial.html
    .. [3] http://functions.wolfram.com/Polynomials/JacobiP/
    """
    nfactor = (S(2)**(a + b + 1) * (gamma(n + a + 1) * gamma(n + b + 1))
               / (2*n + a + b + 1) / (factorial(n) * gamma(n + a + b + 1)))

    return jacobi(n, a, b, x) / sqrt(nfactor)


#----------------------------------------------------------------------------
# Gegenbauer polynomials
#


class gegenbauer(OrthogonalPolynomial):
    r"""
    Gegenbauer polynomial :math:`C_n^{\left(\alpha\right)}(x)`

    gegenbauer(n, alpha, x) gives the nth Gegenbauer polynomial
    in x, :math:`C_n^{\left(\alpha\right)}(x)`.

    The Gegenbauer polynomials are orthogonal on :math:`[-1, 1]` with
    respect to the weight :math:`\left(1-x^2\right)^{\alpha-\frac{1}{2}}`.

    Examples
    ========

    >>> from sympy import gegenbauer, conjugate, diff
    >>> from sympy.abc import n,a,x
    >>> gegenbauer(0, a, x)
    1
    >>> gegenbauer(1, a, x)
    2*a*x
    >>> gegenbauer(2, a, x)
    -a + x**2*(2*a**2 + 2*a)
    >>> gegenbauer(3, a, x)
    x**3*(4*a**3/3 + 4*a**2 + 8*a/3) + x*(-2*a**2 - 2*a)

    >>> gegenbauer(n, a, x)
    gegenbauer(n, a, x)
    >>> gegenbauer(n, a, -x)
    (-1)**n*gegenbauer(n, a, x)

    >>> gegenbauer(n, a, 0)
    2**n*sqrt(pi)*gamma(a + n/2)/(gamma(a)*gamma(-n/2 + 1/2)*gamma(n + 1))
    >>> gegenbauer(n, a, 1)
    gamma(2*a + n)/(gamma(2*a)*gamma(n + 1))

    >>> conjugate(gegenbauer(n, a, x))
    gegenbauer(n, conjugate(a), conjugate(x))

    >>> diff(gegenbauer(n, a, x), x)
    2*a*gegenbauer(n - 1, a + 1, x)

    See Also
    ========

    jacobi,
    chebyshevt_root, chebyshevu, chebyshevu_root,
    legendre, assoc_legendre,
    hermite,
    laguerre, assoc_laguerre,
    sympy.polys.orthopolys.jacobi_poly
    sympy.polys.orthopolys.gegenbauer_poly
    sympy.polys.orthopolys.chebyshevt_poly
    sympy.polys.orthopolys.chebyshevu_poly
    sympy.polys.orthopolys.hermite_poly
    sympy.polys.orthopolys.legendre_poly
    sympy.polys.orthopolys.laguerre_poly

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Gegenbauer_polynomials
    .. [2] http://mathworld.wolfram.com/GegenbauerPolynomial.html
    .. [3] http://functions.wolfram.com/Polynomials/GegenbauerC3/
    """


    @classmethod
    def eval(cls, n, a, x):
        # For negative n the polynomials vanish
        # See http://functions.wolfram.com/Polynomials/GegenbauerC3/03/01/03/0012/
        if n.is_negative:
            return S.Zero

        # Some special values for fixed a
        if a == S.Half:
            return legendre(n, x)
        elif a == S.One:
            return chebyshevu(n, x)
        elif a == S.NegativeOne:
            return S.Zero

        if not n.is_Number:
            # Handle this before the general sign extraction rule
            if x == S.NegativeOne:
                if (C.re(a) > S.Half) is True:
                    return S.ComplexInfinity
                else:
                    # No sec function available yet
                    #return (C.cos(S.Pi*(a+n)) * C.sec(S.Pi*a) * C.gamma(2*a+n) /
                    #            (C.gamma(2*a) * C.gamma(n+1)))
                    return None

            # Symbolic result C^a_n(x)
            # C^a_n(-x)  --->  (-1)**n * C^a_n(x)
            if x.could_extract_minus_sign():
                return S.NegativeOne**n * gegenbauer(n, a, -x)
            # We can evaluate for some special values of x
            if x == S.Zero:
                return (2**n * sqrt(S.Pi) * C.gamma(a + S.Half*n) /
                        (C.gamma((1 - n)/2) * C.gamma(n + 1) * C.gamma(a)) )
            if x == S.One:
                return C.gamma(2*a + n) / (C.gamma(2*a) * C.gamma(n + 1))
            elif x == S.Infinity:
                if n.is_positive:
                    return C.RisingFactorial(a, n) * S.Infinity
        else:
            # n is a given fixed integer, evaluate into polynomial
            return gegenbauer_poly(n, a, x)

    def fdiff(self, argindex=3):
        if argindex == 1:
            # Diff wrt n
            raise ArgumentIndexError(self, argindex)
        elif argindex == 2:
            # Diff wrt a
            n, a, x = self.args
            k = C.Dummy("k")
            factor1 = 2 * (1 + (-1)**(n - k)) * (k + a) / ((k +
                           n + 2*a) * (n - k))
            factor2 = 2*(k + 1) / ((k + 2*a) * (2*k + 2*a + 1)) + \
                2 / (k + n + 2*a)
            kern = factor1*gegenbauer(k, a, x) + factor2*gegenbauer(n, a, x)
            return C.Sum(kern, (k, 0, n - 1))
        elif argindex == 3:
            # Diff wrt x
            n, a, x = self.args
            return 2*a*gegenbauer(n - 1, a + 1, x)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_polynomial(self, n, a, x):
        k = C.Dummy("k")
        kern = ((-1)**k * C.RisingFactorial(a, n - k) * (2*x)**(n - 2*k) /
                (C.factorial(k) * C.factorial(n - 2*k)))
        return C.Sum(kern, (k, 0, C.floor(n/2)))

    def _eval_conjugate(self):
        n, a, x = self.args
        return self.func(n, a.conjugate(), x.conjugate())

#----------------------------------------------------------------------------
# Chebyshev polynomials of first and second kind
#


class chebyshevt(OrthogonalPolynomial):
    r"""
    Chebyshev polynomial of the first kind, :math:`T_n(x)`

    chebyshevt(n, x) gives the nth Chebyshev polynomial (of the first
    kind) in x, :math:`T_n(x)`.

    The Chebyshev polynomials of the first kind are orthogonal on
    :math:`[-1, 1]` with respect to the weight :math:`\frac{1}{\sqrt{1-x^2}}`.

    Examples
    ========

    >>> from sympy import chebyshevt, diff
    >>> from sympy.abc import n,x
    >>> chebyshevt(0, x)
    1
    >>> chebyshevt(1, x)
    x
    >>> chebyshevt(2, x)
    2*x**2 - 1

    >>> chebyshevt(n, x)
    chebyshevt(n, x)
    >>> chebyshevt(n, -x)
    (-1)**n*chebyshevt(n, x)
    >>> chebyshevt(-n, x)
    chebyshevt(n, x)

    >>> chebyshevt(n, 0)
    cos(pi*n/2)
    >>> chebyshevt(n, -1)
    (-1)**n

    >>> diff(chebyshevt(n, x), x)
    n*chebyshevu(n - 1, x)

    See Also
    ========

    jacobi, gegenbauer,
    chebyshevt_root, chebyshevu, chebyshevu_root,
    legendre, assoc_legendre,
    hermite,
    laguerre, assoc_laguerre,
    sympy.polys.orthopolys.jacobi_poly
    sympy.polys.orthopolys.gegenbauer_poly
    sympy.polys.orthopolys.chebyshevt_poly
    sympy.polys.orthopolys.chebyshevu_poly
    sympy.polys.orthopolys.hermite_poly
    sympy.polys.orthopolys.legendre_poly
    sympy.polys.orthopolys.laguerre_poly

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Chebyshev_polynomial
    .. [2] http://mathworld.wolfram.com/ChebyshevPolynomialoftheFirstKind.html
    .. [3] http://mathworld.wolfram.com/ChebyshevPolynomialoftheSecondKind.html
    .. [4] http://functions.wolfram.com/Polynomials/ChebyshevT/
    .. [5] http://functions.wolfram.com/Polynomials/ChebyshevU/
    """

    _ortho_poly = staticmethod(chebyshevt_poly)

    @classmethod
    def eval(cls, n, x):
        if not n.is_Number:
            # Symbolic result T_n(x)
            # T_n(-x)  --->  (-1)**n * T_n(x)
            if x.could_extract_minus_sign():
                return S.NegativeOne**n * chebyshevt(n, -x)
            # T_{-n}(x)  --->  T_n(x)
            if n.could_extract_minus_sign():
                return chebyshevt(-n, x)
            # We can evaluate for some special values of x
            if x == S.Zero:
                return C.cos(S.Half * S.Pi * n)
            if x == S.One:
                return S.One
            elif x == S.Infinity:
                return S.Infinity
        else:
            # n is a given fixed integer, evaluate into polynomial
            if n.is_negative:
                # T_{-n}(x) == T_n(x)
                return cls._eval_at_order(-n, x)
            else:
                return cls._eval_at_order(n, x)

    def fdiff(self, argindex=2):
        if argindex == 1:
            # Diff wrt n
            raise ArgumentIndexError(self, argindex)
        elif argindex == 2:
            # Diff wrt x
            n, x = self.args
            return n * chebyshevu(n - 1, x)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_polynomial(self, n, x):
        k = C.Dummy("k")
        kern = C.binomial(n, 2*k) * (x**2 - 1)**k * x**(n - 2*k)
        return C.Sum(kern, (k, 0, C.floor(n/2)))


class chebyshevu(OrthogonalPolynomial):
    r"""
    Chebyshev polynomial of the second kind, :math:`U_n(x)`

    chebyshevu(n, x) gives the nth Chebyshev polynomial of the second
    kind in x, :math:`U_n(x)`.

    The Chebyshev polynomials of the second kind are orthogonal on
    :math:`[-1, 1]` with respect to the weight :math:`\sqrt{1-x^2}`.

    Examples
    ========

    >>> from sympy import chebyshevu, diff
    >>> from sympy.abc import n,x
    >>> chebyshevu(0, x)
    1
    >>> chebyshevu(1, x)
    2*x
    >>> chebyshevu(2, x)
    4*x**2 - 1

    >>> chebyshevu(n, x)
    chebyshevu(n, x)
    >>> chebyshevu(n, -x)
    (-1)**n*chebyshevu(n, x)
    >>> chebyshevu(-n, x)
    -chebyshevu(n - 2, x)

    >>> chebyshevu(n, 0)
    cos(pi*n/2)
    >>> chebyshevu(n, 1)
    n + 1

    >>> diff(chebyshevu(n, x), x)
    (-x*chebyshevu(n, x) + (n + 1)*chebyshevt(n + 1, x))/(x**2 - 1)

    See Also
    ========

    jacobi, gegenbauer,
    chebyshevt, chebyshevt_root, chebyshevu_root,
    legendre, assoc_legendre,
    hermite,
    laguerre, assoc_laguerre,
    sympy.polys.orthopolys.jacobi_poly
    sympy.polys.orthopolys.gegenbauer_poly
    sympy.polys.orthopolys.chebyshevt_poly
    sympy.polys.orthopolys.chebyshevu_poly
    sympy.polys.orthopolys.hermite_poly
    sympy.polys.orthopolys.legendre_poly
    sympy.polys.orthopolys.laguerre_poly

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Chebyshev_polynomial
    .. [2] http://mathworld.wolfram.com/ChebyshevPolynomialoftheFirstKind.html
    .. [3] http://mathworld.wolfram.com/ChebyshevPolynomialoftheSecondKind.html
    .. [4] http://functions.wolfram.com/Polynomials/ChebyshevT/
    .. [5] http://functions.wolfram.com/Polynomials/ChebyshevU/
    """

    _ortho_poly = staticmethod(chebyshevu_poly)

    @classmethod
    def eval(cls, n, x):
        if not n.is_Number:
            # Symbolic result U_n(x)
            # U_n(-x)  --->  (-1)**n * U_n(x)
            if x.could_extract_minus_sign():
                return S.NegativeOne**n * chebyshevu(n, -x)
            # U_{-n}(x)  --->  -U_{n-2}(x)
            if n.could_extract_minus_sign():
                if n == S.NegativeOne:
                    return S.Zero
                else:
                    return -chebyshevu(-n - 2, x)
            # We can evaluate for some special values of x
            if x == S.Zero:
                return C.cos(S.Half * S.Pi * n)
            if x == S.One:
                return S.One + n
            elif x == S.Infinity:
                return S.Infinity
        else:
            # n is a given fixed integer, evaluate into polynomial
            if n.is_negative:
                # U_{-n}(x)  --->  -U_{n-2}(x)
                if n == S.NegativeOne:
                    return S.Zero
                else:
                    return -cls._eval_at_order(-n - 2, x)
            else:
                return cls._eval_at_order(n, x)

    def fdiff(self, argindex=2):
        if argindex == 1:
            # Diff wrt n
            raise ArgumentIndexError(self, argindex)
        elif argindex == 2:
            # Diff wrt x
            n, x = self.args
            return ((n + 1) * chebyshevt(n + 1, x) - x * chebyshevu(n, x)) / (x**2 - 1)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_polynomial(self, n, x):
        k = C.Dummy("k")
        kern = S.NegativeOne**k * C.factorial(
            n - k) * (2*x)**(n - 2*k) / (C.factorial(k) * C.factorial(n - 2*k))
        return C.Sum(kern, (k, 0, C.floor(n/2)))


class chebyshevt_root(Function):
    r"""
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

    jacobi, gegenbauer,
    chebyshevt, chebyshevu, chebyshevu_root,
    legendre, assoc_legendre,
    hermite,
    laguerre, assoc_laguerre,
    sympy.polys.orthopolys.jacobi_poly
    sympy.polys.orthopolys.gegenbauer_poly
    sympy.polys.orthopolys.chebyshevt_poly
    sympy.polys.orthopolys.chebyshevu_poly
    sympy.polys.orthopolys.hermite_poly
    sympy.polys.orthopolys.legendre_poly
    sympy.polys.orthopolys.laguerre_poly
    """


    @classmethod
    def eval(cls, n, k):
        if not ((0 <= k) is (k < n) is True):
            raise ValueError("must have 0 <= k < n, "
                "got k = %s and n = %s" % (k, n))
        return C.cos(S.Pi*(2*k + 1)/(2*n))


class chebyshevu_root(Function):
    r"""
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

    chebyshevt, chebyshevt_root, chebyshevu,
    legendre, assoc_legendre,
    hermite,
    laguerre, assoc_laguerre,
    sympy.polys.orthopolys.jacobi_poly
    sympy.polys.orthopolys.gegenbauer_poly
    sympy.polys.orthopolys.chebyshevt_poly
    sympy.polys.orthopolys.chebyshevu_poly
    sympy.polys.orthopolys.hermite_poly
    sympy.polys.orthopolys.legendre_poly
    sympy.polys.orthopolys.laguerre_poly
    """


    @classmethod
    def eval(cls, n, k):
        if not ((0 <= k) is (k < n) is True):
            raise ValueError("must have 0 <= k < n, "
                "got k = %s and n = %s" % (k, n))
        return C.cos(S.Pi*(k + 1)/(n + 1))

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

    jacobi, gegenbauer,
    chebyshevt, chebyshevt_root, chebyshevu, chebyshevu_root,
    assoc_legendre,
    hermite,
    laguerre, assoc_laguerre,
    sympy.polys.orthopolys.jacobi_poly
    sympy.polys.orthopolys.gegenbauer_poly
    sympy.polys.orthopolys.chebyshevt_poly
    sympy.polys.orthopolys.chebyshevu_poly
    sympy.polys.orthopolys.hermite_poly
    sympy.polys.orthopolys.legendre_poly
    sympy.polys.orthopolys.laguerre_poly

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Legendre_polynomial
    .. [2] http://mathworld.wolfram.com/LegendrePolynomial.html
    .. [3] http://functions.wolfram.com/Polynomials/LegendreP/
    .. [4] http://functions.wolfram.com/Polynomials/LegendreP2/
    """

    _ortho_poly = staticmethod(legendre_poly)

    @classmethod
    def eval(cls, n, x):
        if not n.is_Number:
            # Symbolic result L_n(x)
            # L_n(-x)  --->  (-1)**n * L_n(x)
            if x.could_extract_minus_sign():
                return S.NegativeOne**n * legendre(n, -x)
            # L_{-n}(x)  --->  L_{n-1}(x)
            if n.could_extract_minus_sign():
                return legendre(-n - S.One, x)
            # We can evaluate for some special values of x
            if x == S.Zero:
                return sqrt(S.Pi)/(C.gamma(S.Half - n/2)*C.gamma(S.One + n/2))
            elif x == S.One:
                return S.One
            elif x == S.Infinity:
                return S.Infinity
        else:
            # n is a given fixed integer, evaluate into polynomial
            if n.is_negative:
                raise ValueError(
                    "The index n must be nonnegative integer (got %r)" % n)
            else:
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

    jacobi, gegenbauer,
    chebyshevt, chebyshevt_root, chebyshevu, chebyshevu_root,
    legendre,
    hermite,
    laguerre, assoc_laguerre,
    sympy.polys.orthopolys.jacobi_poly
    sympy.polys.orthopolys.gegenbauer_poly
    sympy.polys.orthopolys.chebyshevt_poly
    sympy.polys.orthopolys.chebyshevu_poly
    sympy.polys.orthopolys.hermite_poly
    sympy.polys.orthopolys.legendre_poly
    sympy.polys.orthopolys.laguerre_poly

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Associated_Legendre_polynomials
    .. [2] http://mathworld.wolfram.com/LegendrePolynomial.html
    .. [3] http://functions.wolfram.com/Polynomials/LegendreP/
    .. [4] http://functions.wolfram.com/Polynomials/LegendreP2/
    """


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
            return 1/(x**2 - 1)*(x*n*assoc_legendre(n, m, x) - (m + n)*assoc_legendre(n - 1, m, x))
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_polynomial(self, n, m, x):
        k = C.Dummy("k")
        kern = C.factorial(2*n - 2*k)/(2**n*C.factorial(n - k)*C.factorial(
            k)*C.factorial(n - 2*k - m))*(-1)**k*x**(n - m - 2*k)
        return (1 - x**2)**(m/2) * C.Sum(kern, (k, 0, C.floor((n - m)*S.Half)))

#----------------------------------------------------------------------------
# Hermite polynomials
#


class hermite(OrthogonalPolynomial):
    r"""
    hermite(n, x) gives the nth Hermite polynomial in x, :math:`H_n(x)`

    The Hermite polynomials are orthogonal on :math:`(-\infty, \infty)`
    with respect to the weight :math:`\exp\left(-\frac{x^2}{2}\right)`.

    Examples
    ========

    >>> from sympy import hermite, diff
    >>> from sympy.abc import x, n
    >>> hermite(0, x)
    1
    >>> hermite(1, x)
    2*x
    >>> hermite(2, x)
    4*x**2 - 2
    >>> hermite(n, x)
    hermite(n, x)
    >>> diff(hermite(n,x), x)
    2*n*hermite(n - 1, x)
    >>> diff(hermite(n,x), x)
    2*n*hermite(n - 1, x)
    >>> hermite(n, -x)
    (-1)**n*hermite(n, x)

    See Also
    ========

    jacobi, gegenbauer,
    chebyshevt, chebyshevt_root, chebyshevu, chebyshevu_root,
    legendre, assoc_legendre,
    laguerre, assoc_laguerre,
    sympy.polys.orthopolys.jacobi_poly
    sympy.polys.orthopolys.gegenbauer_poly
    sympy.polys.orthopolys.chebyshevt_poly
    sympy.polys.orthopolys.chebyshevu_poly
    sympy.polys.orthopolys.hermite_poly
    sympy.polys.orthopolys.legendre_poly
    sympy.polys.orthopolys.laguerre_poly

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Hermite_polynomial
    .. [2] http://mathworld.wolfram.com/HermitePolynomial.html
    .. [3] http://functions.wolfram.com/Polynomials/HermiteH/
    """

    _ortho_poly = staticmethod(hermite_poly)

    @classmethod
    def eval(cls, n, x):
        if not n.is_Number:
            # Symbolic result H_n(x)
            # H_n(-x)  --->  (-1)**n * H_n(x)
            if x.could_extract_minus_sign():
                return S.NegativeOne**n * hermite(n, -x)
            # We can evaluate for some special values of x
            if x == S.Zero:
                return 2**n * sqrt(S.Pi) / C.gamma((S.One - n)/2)
            elif x == S.Infinity:
                return S.Infinity
        else:
            # n is a given fixed integer, evaluate into polynomial
            if n.is_negative:
                raise ValueError(
                    "The index n must be nonnegative integer (got %r)" % n)
            else:
                return cls._eval_at_order(n, x)

    def fdiff(self, argindex=2):
        if argindex == 1:
            # Diff wrt n
            raise ArgumentIndexError(self, argindex)
        elif argindex == 2:
            # Diff wrt x
            n, x = self.args
            return 2*n*hermite(n - 1, x)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_polynomial(self, n, x):
        k = C.Dummy("k")
        kern = (-1)**k / (C.factorial(k)*C.factorial(n - 2*k)) * (2*x)**(n - 2*k)
        return C.factorial(n)*C.Sum(kern, (k, 0, C.floor(n/2)))

#----------------------------------------------------------------------------
# Laguerre polynomials
#


class laguerre(OrthogonalPolynomial):
    r"""
    Returns the nth Laguerre polynomial in x, :math:`L_n(x)`.

    Parameters
    ==========

    n : int
        Degree of Laguerre polynomial. Must be ``n >= 0``.

    Examples
    ========

    >>> from sympy import laguerre, diff
    >>> from sympy.abc import x, n
    >>> laguerre(0, x)
    1
    >>> laguerre(1, x)
    -x + 1
    >>> laguerre(2, x)
    x**2/2 - 2*x + 1
    >>> laguerre(3, x)
    -x**3/6 + 3*x**2/2 - 3*x + 1

    >>> laguerre(n, x)
    laguerre(n, x)

    >>> diff(laguerre(n, x), x)
    -assoc_laguerre(n - 1, 1, x)

    See Also
    ========

    jacobi, gegenbauer,
    chebyshevt, chebyshevt_root, chebyshevu, chebyshevu_root,
    legendre, assoc_legendre,
    hermite,
    assoc_laguerre,
    sympy.polys.orthopolys.jacobi_poly
    sympy.polys.orthopolys.gegenbauer_poly
    sympy.polys.orthopolys.chebyshevt_poly
    sympy.polys.orthopolys.chebyshevu_poly
    sympy.polys.orthopolys.hermite_poly
    sympy.polys.orthopolys.legendre_poly
    sympy.polys.orthopolys.laguerre_poly

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Laguerre_polynomial
    .. [2] http://mathworld.wolfram.com/LaguerrePolynomial.html
    .. [3] http://functions.wolfram.com/Polynomials/LaguerreL/
    .. [4] http://functions.wolfram.com/Polynomials/LaguerreL3/
    """

    @classmethod
    def eval(cls, n, x):
        if not n.is_Number:
            # Symbolic result L_n(x)
            # L_{n}(-x)  --->  exp(-x) * L_{-n-1}(x)
            # L_{-n}(x)  --->  exp(x) * L_{n-1}(-x)
            if n.could_extract_minus_sign():
                return C.exp(x) * laguerre(n - 1, -x)
            # We can evaluate for some special values of x
            if x == S.Zero:
                return S.One
            elif x == S.NegativeInfinity:
                return S.Infinity
            elif x == S.Infinity:
                return S.NegativeOne**n * S.Infinity
        else:
            # n is a given fixed integer, evaluate into polynomial
            if n.is_negative:
                raise ValueError(
                    "The index n must be nonnegative integer (got %r)" % n)
            else:
                return laguerre_poly(n, x, 0)

    def fdiff(self, argindex=2):
        if argindex == 1:
            # Diff wrt n
            raise ArgumentIndexError(self, argindex)
        elif argindex == 2:
            # Diff wrt x
            n, x = self.args
            return -assoc_laguerre(n - 1, 1, x)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_polynomial(self, n, x):
        # TODO: Should make sure n is in N_0
        k = C.Dummy("k")
        kern = C.RisingFactorial(-n, k) / C.factorial(k)**2 * x**k
        return C.Sum(kern, (k, 0, n))


class assoc_laguerre(OrthogonalPolynomial):
    r"""
    Returns the nth generalized Laguerre polynomial in x, :math:`L_n(x)`.

    Parameters
    ==========

    n : int
        Degree of Laguerre polynomial. Must be ``n >= 0``.

    alpha : Expr
        Arbitrary expression. For ``alpha=0`` regular Laguerre
        polynomials will be generated.

    Examples
    ========

    >>> from sympy import assoc_laguerre, diff
    >>> from sympy.abc import x, n, a
    >>> assoc_laguerre(0, a, x)
    1
    >>> assoc_laguerre(1, a, x)
    a - x + 1
    >>> assoc_laguerre(2, a, x)
    a**2/2 + 3*a/2 + x**2/2 + x*(-a - 2) + 1
    >>> assoc_laguerre(3, a, x)
    a**3/6 + a**2 + 11*a/6 - x**3/6 + x**2*(a/2 + 3/2) +
        x*(-a**2/2 - 5*a/2 - 3) + 1

    >>> assoc_laguerre(n, a, 0)
    binomial(a + n, a)

    >>> assoc_laguerre(n, a, x)
    assoc_laguerre(n, a, x)

    >>> assoc_laguerre(n, 0, x)
    laguerre(n, x)

    >>> diff(assoc_laguerre(n, a, x), x)
    -assoc_laguerre(n - 1, a + 1, x)

    >>> diff(assoc_laguerre(n, a, x), a)
    Sum(assoc_laguerre(_k, a, x)/(-a + n), (_k, 0, n - 1))

    See Also
    ========

    jacobi, gegenbauer,
    chebyshevt, chebyshevt_root, chebyshevu, chebyshevu_root,
    legendre, assoc_legendre,
    hermite,
    laguerre,
    sympy.polys.orthopolys.jacobi_poly
    sympy.polys.orthopolys.gegenbauer_poly
    sympy.polys.orthopolys.chebyshevt_poly
    sympy.polys.orthopolys.chebyshevu_poly
    sympy.polys.orthopolys.hermite_poly
    sympy.polys.orthopolys.legendre_poly
    sympy.polys.orthopolys.laguerre_poly

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Laguerre_polynomial#Assoc_laguerre_polynomials
    .. [2] http://mathworld.wolfram.com/AssociatedLaguerrePolynomial.html
    .. [3] http://functions.wolfram.com/Polynomials/LaguerreL/
    .. [4] http://functions.wolfram.com/Polynomials/LaguerreL3/
    """


    @classmethod
    def eval(cls, n, alpha, x):
        # L_{n}^{0}(x)  --->  L_{n}(x)
        if alpha == S.Zero:
            return laguerre(n, x)

        if not n.is_Number:
            # We can evaluate for some special values of x
            if x == S.Zero:
                return C.binomial(n + alpha, alpha)
            elif x == S.Infinity and n > S.Zero:
                return S.NegativeOne**n * S.Infinity
            elif x == S.NegativeInfinity and n > S.Zero:
                return S.Infinity
        else:
            # n is a given fixed integer, evaluate into polynomial
            if n.is_negative:
                raise ValueError(
                    "The index n must be nonnegative integer (got %r)" % n)
            else:
                return laguerre_poly(n, x, alpha)

    def fdiff(self, argindex=3):
        if argindex == 1:
            # Diff wrt n
            raise ArgumentIndexError(self, argindex)
        elif argindex == 2:
            # Diff wrt alpha
            n, alpha, x = self.args
            k = C.Dummy("k")
            return C.Sum(assoc_laguerre(k, alpha, x) / (n - alpha), (k, 0, n - 1))
        elif argindex == 3:
            # Diff wrt x
            n, alpha, x = self.args
            return -assoc_laguerre(n - 1, alpha + 1, x)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_polynomial(self, n, alpha, x):
        # TODO: Should make sure n is in N_0
        k = C.Dummy("k")
        kern = C.RisingFactorial(
            -n, k) / (C.gamma(k + alpha + 1) * C.factorial(k)) * x**k
        return C.gamma(n + alpha + 1) / C.factorial(n) * C.Sum(kern, (k, 0, n))


#----------------------------------------------------------------------------
# Charlier polynomials
#

class charlier(OrthogonalPolynomial):
    r"""
    The Charlier polynomial in x, :math:`C_n(a, x)`.

    .. math::
        C_n(a, x) := {}_2F_0\left(
                     \begin{matrix} -n, -x \\ - \end{matrix}
                     \middle| -\frac{1}{a}\right)

    The orthogonality condition is

    .. math::
        \sum_{x=0}^\infty \frac{a^x}{x!} C_n(a, x) C_m(a, x)
        = a^{-n} e^a n! \delta_{n,m}

    with :math:`a > 0`.

    Examples
    ========

    >>> from sympy import Symbol, charlier
    >>> x = Symbol('x')
    >>> n = Symbol("n")
    >>> a = Symbol("a")

    >>> C = charlier(n, a, x)
    >>> C
    charlier(n, a, x)

    >>> from sympy import hyper
    >>> C.rewrite(hyper)
    hyper((-n, -x), (), -1/a)

    See Also
    ========

    jacobi, gegenbauer,
    chebyshevt, chebyshevt_root, chebyshevu, chebyshevu_root,
    legendre, assoc_legendre,
    hermite,
    assoc_laguerre,
    sympy.polys.orthopolys.jacobi_poly
    sympy.polys.orthopolys.gegenbauer_poly
    sympy.polys.orthopolys.chebyshevt_poly
    sympy.polys.orthopolys.chebyshevu_poly
    sympy.polys.orthopolys.hermite_poly
    sympy.polys.orthopolys.legendre_poly
    sympy.polys.orthopolys.laguerre_poly

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Charlier_polynomials
    .. [2] http://dlmf.nist.gov/18.19
    """

    @classmethod
    def eval(cls, n, a, x):
        if n.is_Number:
            if n.is_negative:
                raise ValueError("The index n must be nonnegative integer (got %r)" % n)
        if a.is_Number:
            if a.is_nonnegative:
                raise ValueError("The index a must be a positive number (got %r)" % n)


    def _eval_rewrite_as_hyper(self, n, a, x):
        return hyper([-n, -x], [], -1/a)


#----------------------------------------------------------------------------
# Meixner polynomials
#

class meixner(OrthogonalPolynomial):
    r"""
    The Meixner polynomial in x, :math:`M_n(\beta, c, x)`.

    .. math::
        M_n(\beta, c, x) := {}_2F_1\left(\begin{matrix} -n, -x \\ \beta \end{matrix} \middle| 1 - \frac{1}{c}\right)

    The orthogonality condition is

    .. math::
        \sum_{x=0}^\infty \frac{(\beta)_x}{x!} c^x M_n(\beta, c, x) M_m(\beta, c, x)
        = \frac{c^{-n} n!}{(\beta)_n (1-c)^\beta} \delta_{n,m}

    with :math:`\beta > 0` and :math:`0 < c < 1`.

    Examples
    ========

    >>> from sympy import Symbol, meixner
    >>> x = Symbol('x')
    >>> n = Symbol("n")
    >>> b = Symbol("beta")
    >>> c = Symbol("c")

    >>> M = meixner(n, b, c, x)
    >>> M
    meixner(n, beta, c, x)

    >>> from sympy import hyper
    >>> M.rewrite(hyper)
    hyper((-n, -x), (beta,), 1 - 1/c)

    See Also
    ========

    jacobi, gegenbauer,
    chebyshevt, chebyshevt_root, chebyshevu, chebyshevu_root,
    legendre, assoc_legendre,
    hermite,
    assoc_laguerre,
    sympy.polys.orthopolys.jacobi_poly
    sympy.polys.orthopolys.gegenbauer_poly
    sympy.polys.orthopolys.chebyshevt_poly
    sympy.polys.orthopolys.chebyshevu_poly
    sympy.polys.orthopolys.hermite_poly
    sympy.polys.orthopolys.legendre_poly
    sympy.polys.orthopolys.laguerre_poly

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Meixner_polynomials
    .. [2] http://dlmf.nist.gov/18.19
    """

    @classmethod
    def eval(cls, n, b, c, x):
        if n.is_Number:
            if n.is_negative:
                raise ValueError("The index n must be nonnegative integer (got %r)" % n)
        if b.is_Number:
            if b.is_nonnegative:
                raise ValueError("The index beta must be a positive number (got %r)" % n)


    def _eval_rewrite_as_hyper(self, n, b, c, x):
        return hyper([-n, -x], [b], 1 - 1/c)

    def _eval_rewrite_as_jacobi(self, n, b, c, x):
        return jacobi(n, b - 1, -n - b - x, (2 - c)/c) * factorial(n) / RisingFactorial(b, n)


#----------------------------------------------------------------------------
# Krawtchouk polynomials
#

class krawtchouk(OrthogonalPolynomial):
    r"""
    The Krawtchouk polynomial in x, :math:`K_n(p, N, x)`.

    .. math::
        K_n(p, N, x) := {}_2F_1\left(\begin{matrix} -n, -x \\ -N \end{matrix} \middle| \frac{1}{p} \right)

    and :math:`n = 0, 1, 2, \ldots, N`.

    The orthogonality condition is

    .. math::
        \sum_{x=0}^N \binom{N}{x} p^x (1-p)^{N-x} K_n(p, N, x) K_m(p, N, x)
        = \frac{(-1)^n n!}{(-N)_n} \left(\frac{1-p}{p}\right)^n \delta_{m,n}

    with :math:`0 < p < 1`.

    Examples
    ========

    >>> from sympy import Symbol, krawtchouk
    >>> x = Symbol('x')
    >>> n = Symbol("n")
    >>> p = Symbol("p")
    >>> N = Symbol("N")

    >>> K = krawtchouk(n, p, N, x)
    >>> K
    krawtchouk(n, p, N, x)

    >>> from sympy import hyper
    >>> K.rewrite(hyper)
    hyper((-n, -x), (-N,), 1/p)

    See Also
    ========

    jacobi, gegenbauer,
    chebyshevt, chebyshevt_root, chebyshevu, chebyshevu_root,
    legendre, assoc_legendre,
    hermite,
    assoc_laguerre,
    sympy.polys.orthopolys.jacobi_poly
    sympy.polys.orthopolys.gegenbauer_poly
    sympy.polys.orthopolys.chebyshevt_poly
    sympy.polys.orthopolys.chebyshevu_poly
    sympy.polys.orthopolys.hermite_poly
    sympy.polys.orthopolys.legendre_poly
    sympy.polys.orthopolys.laguerre_poly

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Krawtchouk_polynomials
    .. [2] http://dlmf.nist.gov/18.19
    .. [3] http://mathworld.wolfram.com/KrawtchoukPolynomial.html
    """

    @classmethod
    def eval(cls, n, p, N, x):
        if n.is_Number:
            if n.is_negative:
                raise ValueError("The index n must be nonnegative integer (got %r)" % n)

    def _eval_rewrite_as_hyper(self, n, p, N, x):
        return hyper([-n, -x], [-N], 1/p)

    def _eval_rewrite_as_meixner(self, n, p, N, x):
        return meixner(n, -N, p/(p - 1), x)

#----------------------------------------------------------------------------
# Meixner-Pollaczek polynomials
#

class meixner_pollaczek(OrthogonalPolynomial):
    r"""
    The Meixner-Pollaczek polynomial in x, :math:`P_n(\lambda, \phi, x)`.

    .. math::
        P_n(\lambda, \phi, x) := \frac{(2\lambda)_n}{n!} \exp(\imath n \phi)
        {}_2F_1\left(\begin{matrix} -n, \lambda + \imath x \\ 2 \lambda \end{matrix}
        \middle| 1 - \exp(-2\imath \phi) \right)

    The orthogonality condition is


    with :math:`\lambda > 0` and :math:`0 < \phi < \pi`

    Examples
    ========

    >>> from sympy import Symbol, meixner_pollaczek
    >>> x = Symbol('x')
    >>> n = Symbol("n")
    >>> l = Symbol("l")
    >>> p = Symbol("p")

    >>> M = meixner_pollaczek(n, l, p, x)
    >>> M
    meixner_pollaczek(n, l, p, x)

    >>> from sympy import hyper
    >>> M.rewrite(hyper)
    exp(I*n*p)*RisingFactorial(2*l, n)*hyper((-n, l + I*x), (2*l,), 1 - exp(-2*I*p))/factorial(n)

    See Also
    ========

    jacobi, gegenbauer,
    chebyshevt, chebyshevt_root, chebyshevu, chebyshevu_root,
    legendre, assoc_legendre,
    hermite,
    assoc_laguerre,
    sympy.polys.orthopolys.jacobi_poly
    sympy.polys.orthopolys.gegenbauer_poly
    sympy.polys.orthopolys.chebyshevt_poly
    sympy.polys.orthopolys.chebyshevu_poly
    sympy.polys.orthopolys.hermite_poly
    sympy.polys.orthopolys.legendre_poly
    sympy.polys.orthopolys.laguerre_poly

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Meixner%E2%80%93Pollaczek_polynomials
    .. [2] http://dlmf.nist.gov/18.19
    """

    @classmethod
    def eval(cls, n, l, p, x):
        if n.is_Number:
            if n.is_negative:
                raise ValueError("The index n must be nonnegative integer (got %r)" % n)

    def _eval_rewrite_as_hyper(self, n, l, p, x):
        return (RisingFactorial(2*l, n) / factorial(n) * exp(I*n*p) *
                hyper([-n, l + I*x], [2*l], 1 - exp(-2*I*p)))

#----------------------------------------------------------------------------
# Hahn polynomials
#

class hahn(OrthogonalPolynomial):
    r"""
    The Hahn polynomial in x, :math:`Q_n(\alpha, \beta, N, x)`

    .. math::
        Q_n(\alpha, \beta, N, x) := {}_3F_2\left(
        \begin{matrix}
        -n, n + \alpha + \beta + 1, -x \\
        \alpha + 1, N
        \end{matrix}
        \middle| 1 \right)

    for :math:`n = 0, 1, 2, \ldots, N`.

    Examples
    ========

    >>> from sympy import Symbol, hahn
    >>> x = Symbol('x')
    >>> n = Symbol("n")
    >>> a = Symbol("a")
    >>> b = Symbol("b")
    >>> N = Symbol("N")

    >>> H = hahn(n, a, b, N, x)
    >>> H
    hahn(n, a, b, N, x)

    >>> from sympy import hyper
    >>> H.rewrite(hyper)
    hyper((-n, a + b + n + 1, -x), (a + 1, -N), 1)

    See Also
    ========

    jacobi, gegenbauer,
    chebyshevt, chebyshevt_root, chebyshevu, chebyshevu_root,
    legendre, assoc_legendre,
    hermite,
    assoc_laguerre,
    sympy.polys.orthopolys.jacobi_poly
    sympy.polys.orthopolys.gegenbauer_poly
    sympy.polys.orthopolys.chebyshevt_poly
    sympy.polys.orthopolys.chebyshevu_poly
    sympy.polys.orthopolys.hermite_poly
    sympy.polys.orthopolys.legendre_poly
    sympy.polys.orthopolys.laguerre_poly

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Hahn_polynomials
    .. [2] http://dlmf.nist.gov/18.19
    """

    @classmethod
    def eval(cls, n, a, b, N, x):
        if n.is_Number:
            if n.is_negative:
                raise ValueError("The index n must be nonnegative integer (got %r)" % n)

    def _eval_rewrite_as_hyper(self, n, a, b, N, x):
        return hyper([-n, n + a + b + 1, -x], [a + 1, -N], 1)


#----------------------------------------------------------------------------
# Dual Hahn polynomials
#

class hahn_dual(OrthogonalPolynomial):
    r"""
    The dual Hahn polynomial in x, :math:`R_n(\lambda(x), \gamma, \delta, N, x)`

    .. math::
        R_n(\lambda(x), \gamma, \delta, N, x) := {}_3F_2\left(
        \begin{matrix}
        -n, -x, x + \gamma + \delta + 1 \\
        \gamma + 1, -N
        \end{matrix}
        \middle| 1 \right)

    for :math:`n = 0, 1, 2, \ldots, N` and with
    :math:`\lambda(x) = x (x + \gamma + \delta + 1)`.

    Examples
    ========

    >>> from sympy import Symbol, hahn_dual
    >>> x = Symbol('x')
    >>> n = Symbol("n")
    >>> c = Symbol("c")
    >>> d = Symbol("d")
    >>> N = Symbol("N")

    >>> H = hahn_dual(n, c, d, N, x)
    >>> H
    hahn_dual(n, c, d, N, x)

    >>> from sympy import hyper
    >>> H.rewrite(hyper)    # doctest:+SKIP
    hyper((-n, c/2 + d/2 + sqrt(c**2 + 2*c*d + 2*c + d**2 + 2*d + 4*x + 1)/2 + 1/2,
    c/2 + d/2 - sqrt(c**2 + 2*c*d + 2*c + d**2 + 2*d + 4*x + 1)/2 + 1/2), (c + 1, -N), 1)

    See Also
    ========

    jacobi, gegenbauer,
    chebyshevt, chebyshevt_root, chebyshevu, chebyshevu_root,
    legendre, assoc_legendre,
    hermite,
    assoc_laguerre,
    sympy.polys.orthopolys.jacobi_poly
    sympy.polys.orthopolys.gegenbauer_poly
    sympy.polys.orthopolys.chebyshevt_poly
    sympy.polys.orthopolys.chebyshevu_poly
    sympy.polys.orthopolys.hermite_poly
    sympy.polys.orthopolys.legendre_poly
    sympy.polys.orthopolys.laguerre_poly

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Dual_Hahn_polynomials
    .. [2] http://dlmf.nist.gov/18.19
    """

    @classmethod
    def eval(cls, n, c, d, N, lx):
        if n.is_Number:
            if n.is_negative:
                raise ValueError("The index n must be nonnegative integer (got %r)" % n)

    def _eval_rewrite_as_hyper(self, n, c, d, N, lx):
        # It does not matter which root of l(x) = x * (x + c + d + 1) we take
        x = -d/2 - c/2 - sqrt(d**2 + 2*d*c + 2*d + c**2 + 2*c + 4*lx + 1)/2 - S.Half
        return hyper([-n, -x, x + c + d + 1], [c + 1, -N], 1)

#----------------------------------------------------------------------------
# Continuous Hahn polynomials
#

class hahn_continuous(OrthogonalPolynomial):
    r"""
    The continuous Hahn polynomial in x, :math:`p_n(a, b, c, d, x)`

    .. math::
        p_n(a, b, c, d, x) := \imath^n \frac{(a + c)_n (a + d)_n}{n!}
        {}_3F_2\left(
        \begin{matrix}
        -n, n + a + b + c + d - 1, a + \imath x \\
        a + c, a + d
        \end{matrix}
        \middle| 1 \right)

    Examples
    ========

    >>> from sympy import Symbol, hahn_continuous
    >>> x = Symbol('x')
    >>> n = Symbol("n")
    >>> a = Symbol("a")
    >>> b = Symbol("b")
    >>> c = Symbol("c")
    >>> d = Symbol("d")
    >>> N = Symbol("N")

    >>> H = hahn_continuous(n, a, b, c, d, x)
    >>> H
    hahn_continuous(n, a, b, c, d, x)

    >>> from sympy import hyper
    >>> H.rewrite(hyper)    # doctest:+SKIP
    I**n*RisingFactorial(a + c, n)*RisingFactorial(a + d, n)*
    hyper((-n, a + b + c + d + n - 1, a + I*x), (a + c, a + d), 1)/factorial(n)

    See Also
    ========

    jacobi, gegenbauer,
    chebyshevt, chebyshevt_root, chebyshevu, chebyshevu_root,
    legendre, assoc_legendre,
    hermite,
    assoc_laguerre,
    sympy.polys.orthopolys.jacobi_poly
    sympy.polys.orthopolys.gegenbauer_poly
    sympy.polys.orthopolys.chebyshevt_poly
    sympy.polys.orthopolys.chebyshevu_poly
    sympy.polys.orthopolys.hermite_poly
    sympy.polys.orthopolys.legendre_poly
    sympy.polys.orthopolys.laguerre_poly

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Continuous_Hahn_polynomials
    .. [2] http://dlmf.nist.gov/18.19
    """

    @classmethod
    def eval(cls, n, a, b, c, d, x):
        if n.is_Number:
            if n.is_negative:
                raise ValueError("The index n must be nonnegative integer (got %r)" % n)

    def _eval_rewrite_as_hyper(self, n, a, b, c, d, x):
        return (I**n * RisingFactorial(a + c,n) * RisingFactorial(a + d,n) / factorial(n)
                * hyper([-n, n + a + b + c + d - 1, a + I*x], [a + c, a + d], 1))

#----------------------------------------------------------------------------
# Continuous dual Hahn polynomials
#

class hahn_dual_continuous(OrthogonalPolynomial):
    r"""
    The continuous dual Hahn polynomial in x, :math:`S_n(a, b, c, x^2)`

    .. math::
        S_n(a, b, c, x^2) := (a + b)_n (a + c)_n {}_3F_2\left(
        \begin{matrix}
        -n, a + \imath x, a - \imath x \\
        a + b, a + c
        \end{matrix}
        \middle| 1 \right)

    Examples
    ========

    >>> from sympy import Symbol, hahn_dual_continuous
    >>> x = Symbol('x')
    >>> n = Symbol("n")
    >>> a = Symbol("a")
    >>> b = Symbol("b")
    >>> c = Symbol("c")
    >>> N = Symbol("N")

    >>> H = hahn_dual_continuous(n, a, b, c, x)
    >>> H
    hahn_dual_continuous(n, a, b, c, x)

    >>> from sympy import hyper
    >>> H.rewrite(hyper)    # doctest:+SKIP
    RisingFactorial(a + b, n)*RisingFactorial(a + c, n)*
    hyper((-n, a + I*sqrt(x), a - I*sqrt(x)), (a + b, a + c), 1)

    See Also
    ========

    jacobi, gegenbauer,
    chebyshevt, chebyshevt_root, chebyshevu, chebyshevu_root,
    legendre, assoc_legendre,
    hermite,
    assoc_laguerre,
    sympy.polys.orthopolys.jacobi_poly
    sympy.polys.orthopolys.gegenbauer_poly
    sympy.polys.orthopolys.chebyshevt_poly
    sympy.polys.orthopolys.chebyshevu_poly
    sympy.polys.orthopolys.hermite_poly
    sympy.polys.orthopolys.legendre_poly
    sympy.polys.orthopolys.laguerre_poly

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Continuous_dual_Hahn_polynomials
    .. [2] http://dlmf.nist.gov/18.19
    """

    @classmethod
    def eval(cls, n, a, b, c, lx):
        if n.is_Number:
            if n.is_negative:
                raise ValueError("The index n must be nonnegative integer (got %r)" % n)

    def _eval_rewrite_as_hyper(self, n, a, b, c, lx):
        # It does not matter which root of l(x) = x^2 we take
        x = sqrt(lx)
        return (RisingFactorial(a + b, n) * RisingFactorial(a + c, n)
                * hyper([-n, a + I*x, a - I*x], [a + b, a + c], 1))

#----------------------------------------------------------------------------
# Wilson polynomials
#

class wilson(OrthogonalPolynomial):
    r"""
    The Wilson polynomial in x, :math:`W_n(a, b, c, d, x^2)`

    .. math::
        W_n(a, b, c, d, x^2) := (a + b)_n (a + c)_n (a + d)_n {}_4F_3\left(
        \begin{matrix}
        -n, n + a + b + c + d - 1, a + \imath x, a - \imath x \\
        a + b, a + c, a + d
        \end{matrix}
        \middle| 1 \right)

    Examples
    ========

    >>> from sympy import Symbol, wilson
    >>> x = Symbol('x')
    >>> n = Symbol("n")
    >>> a = Symbol("a")
    >>> b = Symbol("b")
    >>> c = Symbol("c")
    >>> d = Symbol("d")

    >>> W = wilson(n, a, b, c, d, x)
    >>> W
    wilson(n, a, b, c, d, x)

    >>> from sympy import hyper
    >>> W.rewrite(hyper)    # doctest:+SKIP
    RisingFactorial(a + b, n) * RisingFactorial(a + c, n) * RisingFactorial(a + d, n) *
    hyper((-n, a + b + c + d + n + 1, a + I*sqrt(x), a - I*sqrt(x)), (a + b, a + c, a + d), 1)

    See Also
    ========

    jacobi, gegenbauer,
    chebyshevt, chebyshevt_root, chebyshevu, chebyshevu_root,
    legendre, assoc_legendre,
    hermite,
    assoc_laguerre,
    sympy.polys.orthopolys.jacobi_poly
    sympy.polys.orthopolys.gegenbauer_poly
    sympy.polys.orthopolys.chebyshevt_poly
    sympy.polys.orthopolys.chebyshevu_poly
    sympy.polys.orthopolys.hermite_poly
    sympy.polys.orthopolys.legendre_poly
    sympy.polys.orthopolys.laguerre_poly

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Wilson_polynomials
    .. [2] http://dlmf.nist.gov/18.25
    .. [3] http://dlmf.nist.gov/18.26
    .. [4] http://www.encyclopediaofmath.org/index.php?title=Wilson_polynomials
    """

    @classmethod
    def eval(cls, n, a, b, c, d, lx):
        if n.is_Number:
            if n.is_negative:
                raise ValueError("The index n must be nonnegative integer (got %r)" % n)

    def _eval_rewrite_as_hyper(self, n, a, b, c, d, lx):
        # It does not matter which root of l(x) = x^2 we take
        x = sqrt(lx)
        return (RisingFactorial(a + b, n) * RisingFactorial(a + c, n) * RisingFactorial(a + d, n)
                * hyper([-n, n + a + b + c + d + 1, a + I*x, a - I*x], [a + b, a + c, a + d], 1))

#----------------------------------------------------------------------------
# Racah polynomials
#

class racah(OrthogonalPolynomial):
    r"""
    The Racah polynomial in x, :math:`R_n(\lambda(x), \alpha, \beta, \gamma, \delta, x)`

    .. math::
        R_n(\lambda(x), \alpha, \beta, \gamma, \delta, x) := {}_3F_2\left(
        \begin{matrix}
        -n, b + \alpha + \beta + 1, -x, x + \gamma + \delta + 1 \\
        \alpha + 1, \beta + \delta + 1,\gamma + 1
        \end{matrix}
        \middle| 1 \right)

    for :math:`n = 0, 1, 2, \ldots, N` and with
    :math:`\lambda(x) = x (x + \gamma + \delta + 1)`.

    Examples
    ========

    >>> from sympy import Symbol, racah
    >>> x = Symbol('x')
    >>> n = Symbol("n")
    >>> a = Symbol("a")
    >>> b = Symbol("b")
    >>> c = Symbol("c")
    >>> d = Symbol("d")

    >>> R = racah(n, a, b, c, d, x)
    >>> R
    racah(n, a, b, c, d, x)

    >>> from sympy import hyper
    >>> R.rewrite(hyper)    # doctest:+SKIP
    hyper((-n, a + b + n + 1, c/2 + d/2 +
    sqrt(c**2 + 2*c*d + 2*c + d**2 + 2*d + 4*x + 1)/2 + 1/2,
    c/2 + d/2 - sqrt(c**2 + 2*c*d + 2*c + d**2 + 2*d + 4*x + 1)/2
    + 1/2), (a + 1, b + d + 1, c + 1), 1)

    See Also
    ========

    jacobi, gegenbauer,
    chebyshevt, chebyshevt_root, chebyshevu, chebyshevu_root,
    legendre, assoc_legendre,
    hermite,
    assoc_laguerre,
    sympy.polys.orthopolys.jacobi_poly
    sympy.polys.orthopolys.gegenbauer_poly
    sympy.polys.orthopolys.chebyshevt_poly
    sympy.polys.orthopolys.chebyshevu_poly
    sympy.polys.orthopolys.hermite_poly
    sympy.polys.orthopolys.legendre_poly
    sympy.polys.orthopolys.laguerre_poly

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Racah_polynomials
    .. [2] http://dlmf.nist.gov/18.25
    .. [3] http://dlmf.nist.gov/18.26
    """

    @classmethod
    def eval(cls, n, a, b, c, d, lx):
        if n.is_Number:
            if n.is_negative:
                raise ValueError("The index n must be nonnegative integer (got %r)" % n)

    def _eval_rewrite_as_hyper(self, n, a, b, c, d, lx):
        # It does not matter which root of l(x) = x * (x + c + d + 1) we take
        x = -d/2 - c/2 - sqrt(d**2 + 2*d*c + 2*d + c**2 + 2*c + 4*lx + 1)/2 - S.Half
        return hyper([-n, n + a + b + 1, -x, x + c + d + 1], [a + 1, b + d + 1, c + 1], 1)

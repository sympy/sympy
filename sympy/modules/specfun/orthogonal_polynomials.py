from sympy.core.symbol import Symbol
from sympy.core.numbers import Rational, Real, pi
from sympy.core.function import DefinedFunction
#from sympy.modules.trigonometric import cos
from factorials import Function2
import decimal


# Simple implementation of Newton's method for root-finding
def _newton(h, x0, eps):
    x = x0
    prevdiff = 1
    while 1:
        new = x - h(x)
        diff = abs(x - new)
        if diff <= eps:
            break
        prevdiff = diff
        x = new
    return new

class _PolynomialSequence(DefinedFunction):
    nofargs = 2
    _x = Symbol('x')

    def _calc(self, n):
        raise NotImplementedError

    def poly(self, n, x):
        assert n.is_integer and n >= 0
        n = int(n)
        m = len(self._memo)
        if n < m:
            return self._memo[n]
        else:
            for i in range(m, n+1):
                L = self._calc(i)
                L = L.expand()
                self._memo[i] = L
            return self._memo[n]

    def _eval_apply(self, n, x):
        if isinstance(x, Legendre_zero) and x._args[0] == n:
            return 0
        if n.is_integer and n >= 0:
            for k in range(n):
                if x == self._zero_class(n, k):
                    return 0
            #return self.poly().subs(self._x, x)
            return self.poly(n, x).subs(self._x, x)
        #return self


class Legendre(_PolynomialSequence):
    """
    Usage
    =====
        legendre(n, x) - nth Legendre polynomial of x, P_n(x)

    Notes
    =====
        The Legendre polynomials are orthogonal on [-1, 1] with respect
        to the constant weight 1.

        For all n, P_n(1) = 1

        P_n is odd for odd n and even for even n

    Examples
    ========
        >>> x = Symbol('x')
        >>> legendre(3, x)
        -3/2*x+5/2*x**3

    See also
    ========
       External links
       --------------
         U{Wikipedia: Legendre polynomial<http://en.wikipedia.org/wiki/Legendre_polynomial>}
    """
    _memo = {}

    def _calc(self, n):
        if n == 0: return Rational(1)
        if n == 1: return self._x
        return ((2*n-1)*self._x*self._memo[n-1] - (n-1)*self._memo[n-2])/n

legendre = Legendre()


class Legendre_zero(DefinedFunction):
    """
    Usage
    =====
        legendre_zero(n, k) represents the kth zero (counting from zero)
        of the nth Legendre polynomial; that is, if 0 <= k < n,
        legendre(n, legendre_zero(n, k)) == 0.

        All zeros for a given Legendre polynomial are located symmetrically
        around 0 in the open interval (-1, 1). The zeros are indexed from
        left to right.

    Examples
    ========
        >>> legendre(5, legendre_zero(5, 3)) == 0
        True

    """

    def eval(self):
        n, k = self._args
        if n.is_odd and (n-1)/2 == k:
            return Rational(0)
        if n == 2 and k == 0: return -sqrt(Rational(1,3))
        if n == 2 and k == 1: return sqrt(Rational(1,3))
        if n == 3 and k == 0: return -sqrt(Rational(3,5))
        if n == 3 and k == 2: return sqrt(Rational(3,5))
        # We could use SymPy's polynomial root-finding code for higher-degree
        # polynomials, but it might not be helpful to do so by default
        # since the expressions grow extremely complicated
        return self

    def evalf(self, prec=10):
        # Increasing the precision is really just a matter of using
        # a lower epsilon; the problem is that numerical evaluation of
        # polynomials currently doesn't work as it should
        if prec > 10:
            raise NotImplementedError
        eps = 1e-10

        n, k = self._args
        assert 0 <= k < n

        L = lambda x: legendre(n, x)

        t = Symbol('t')
        Ldpol = legendre(n, t).diff(t)
        Ld = lambda x: Ldpol.subs(t, x)

        # Good initial estimate for use with Newton's method
        import math
        x = -math.cos(math.pi*(k+1-0.25)/(n+0.5))

        return _newton(lambda t: L(t)/Ld(t), x, eps)

legendre_zero = Legendre_zero()
legendre._zero_class = legendre_zero


class Chebyshev3(_PolynomialSequence):
    """
    Usage
    =====
        chebyshev(n, x) - nth Chebyshev polynomial (of the first
        kind) of x, T_n(x)

    Notes
    =====
        The Chebyshev polynomials are orthogonal on [-1, 1] with
        respect to the weight 1/sqrt(1-x**2).

    Examples
    ========
        >>> x = Symbol('x')
        >>> chebyshev(3, x)
        -3*x+4*x**3

    See also
    ========
       External links
       --------------
         U{Wikipedia: Chebyshev polynomial<http://en.wikipedia.org/wiki/Chebyshev_polynomial>}
    """
    _memo = {}

    def _calc(self, n):
        if n == 0: return Rational(1)
        if n == 1: return self._x
        return 2*self._x*self._memo[n-1] - self._memo[n-2]

#chebyshev = Chebyshev()


class Chebyshev_zero(DefinedFunction):
    """
    Usage
    =====
        chebyshev_zero(n, k) returns the kth zero (counting from zero)
        of the nth Chebyshev polynomial; that is, if 0 <= k < n,
        chebyshev(n, chebyshev_zero(n, k)) == 0.

        The n,k-th zero is given explicitly by cos(pi*(2*k+1)/(2*n)).
        Due to this simple form, chebyshev_zero always returns an
        explicit expression (unlike legendre_zero).

    Examples
    ========
        >>> chebyshev(5, chebyshev_zero(5, 3)) == 0
        True

    """
    nofargs = 2

    def eval(self):
        n, k = self._args
        return cos(pi*(2*k+1)/(2*n))


chebyshev_zero = Chebyshev_zero()
Chebyshev3._zero_class = chebyshev_zero

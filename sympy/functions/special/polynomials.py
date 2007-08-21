
from sympy.core.basic import Basic, S
from sympy.core.function import Function, DefinedFunction, Apply, Lambda

###############################################################################
############################ ORTHOGONAL POLYNOMIALS ###########################
###############################################################################

class OrthogonalPolynomial(Function):

    def __new__(cls, index, **assumptions):
        index = Basic.sympify(index)
        if index.is_negative:
            raise ValueError("%s index must be nonnegative integer (got %r)" % (cls.__name__, index))
        obj = Basic.__new__(cls, index, **assumptions)
        obj.name = obj.__class__.__name__
        return obj

    def _hashable_content(self):
        return (self.name, self.index)

    @property
    def index(self):
        return self._args[0]

    def torepr(self):
        return '%s(%r)' % (self.__class__.__name__, self.index)

    precedence = Basic.Apply_precedence

    def tostr(self, level=0):
        r = '%s(%s)' % (self.name, self.index)
        p = self.precedence
        if p<=level: r = '(%s)' % r
        return r

    def weight(self):
        raise NotImplementedError('%s.weight()' % (self))

###############################################################################
############# CHEBYSHEV POLYNOMIALS of FIRST and SECOND KIND ##################
###############################################################################

class Chebyshev(OrthogonalPolynomial):
    """
    Chebyshev polynomials of the first kind.
    """

    nofargs = 1
    _cache = {}
    kind = 1
    is_commutative = True

    def as_lambda(self):
        i = self.index
        x = Basic.Symbol('x',dummy=True)
        cls = self.__class__
        if not isinstance(i, Basic.Integer):
            return Lambda((2*x*cls(i-1)(x)-cls(i-2)(x)),x)
        try:
            return cls._cache[i]
        except KeyError:
            pass
        if i==0:
            r = Lambda(Basic.One(),x)
        elif i==1:
            r = Lambda(self.kind * x,x)
        elif i>1:
            r = Lambda((2*x*cls(i-1).as_lambda()(x) - cls(i-2).as_lambda()(x)).expand(), x)
        else:
            raise ValueError("%s index must be nonnegative integer (got %r)" % (self, i))
        cls._cache[i] = r
        return r

    def _eval_expand_basic(self, *args):
        return self.as_lambda()

    def _eval_apply(self, x):
        if isinstance(x, Basic.Number):
            return self.as_lambda()(x)

        if isinstance(x, ApplyChebyshev): # nesting property
            return self.__class__(self.index * x.func.index)(x.args[0])

    def weight(self):
        x = Basic.Symbol('x',dummy=True)
        return Lambda(Basic.Sqrt()(1-x**2),x)

class Chebyshev2(Chebyshev):
    """
    Chebyshev polynomials of the second kind.
    """

    _cache = {}
    kind = 2

    def _eval_apply(self, x):
        if isinstance(x, Basic.Number):
            return self.as_lambda()(x)

class ApplyChebyshev(Apply):

    def _eval_expand_basic(self, *args):
        return self.func.as_lambda()(self.args[0]).expand()

class ApplyChebyshev2(ApplyChebyshev):

    pass

Basic.singleton['Chebyshev'] = lambda: Chebyshev
Basic.singleton['Chebyshev2'] = lambda: Chebyshev2

#### TO BE REFACTORED

# Simple implementation of Newton's method for root-finding
def _newton(h, x0, eps=1e-10):
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
    _x = Basic.Symbol('x')

    def _calc(self, n):
        raise NotImplementedError

    def poly(self, n, x):
        assert n.is_integer and n >= 0
        n = int(n)
        m = len(self._memo)
        if n < m:
            return self._memo[n]
        else:
            for i in xrange(m, n+1):
                L = self._calc(i)
                L = L.expand()
                self._memo[i] = L
            return self._memo[n]

    def _eval_apply(self, n, x):
        #if isinstance(x, Apply) and x.func == self._zero_class and x.args[0] == n:
        #    return Basic.Zero()
        if n.is_integer and n >= 0:
            #for k in xrange(int(n)):
            #    if x == self._zero_class(n, k):
            #        return Basic.Zero()
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
        if n == 0: return Basic.One()
        if n == 1: return self._x
        return ((2*n-1)*self._x*self._memo[n-1] - (n-1)*self._memo[n-2])/n

Basic.singleton['legendre'] = lambda: Legendre

'''
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
    nofargs = 2

    def _eval_apply(self, n, k):
        if n.is_odd and (n-1)/2 == k:
            return Basic.Zero()
        if n == 2 and k == 0: return -sqrt(Rational(1,3))
        if n == 2 and k == 1: return sqrt(Rational(1,3))
        if n == 3 and k == 0: return -sqrt(Rational(3,5))
        if n == 3 and k == 2: return sqrt(Rational(3,5))
        # We could use SymPy's polynomial root-finding code for higher-degree
        # polynomials, but it might not be helpful to do so by default
        # since the expressions grow extremely complicated

class ApplyLegendre_zero(Apply):
    def _eval_evalf(self):
        # Increasing the precision is really just a matter of using
        # a lower epsilon; the problem is that numerical evaluation of
        # polynomials currently doesn't work as it should
        x = Symbol('x')
        t = Symbol('t')

        n, k = self.args
        assert 0 <= k < n

        #L = lambda x: legendre(n, x)
        L = Lambda(legendre(n, x), x)
        Ldpol = legendre(n, t).diff(t)
        Ld = Lambda(Ldpol.subs(t, x), x)

        # Good initial estimate for use with Newton's method
        import math
        x = -math.cos(math.pi*(k+1-0.25)/(n+0.5))

        return _newton(lambda t: L(t)/Ld(t), x)

legendre_zero = Legendre_zero()
Legendre._zero_class = legendre_zero
'''

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
        if n == 0: return Basic.One()
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

    def _eval_apply(self, n, k):
        return S.Cos(pi*(2*k+1)/(2*n))


#chebyshev_zero = Chebyshev_zero()
#Chebyshev3._zero_class = chebyshev_zero

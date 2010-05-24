"""Implementation of :class:`ExpressionDomain` class. """

from sympy.polys.domains.field import Field
from sympy.polys.domains.simpledomain import SimpleDomain
from sympy.polys.domains.characteristiczero import CharacteristicZero

from sympy.core import sympify
from sympy.polys.polyerrors import DomainError

class ExpressionDomain(Field, CharacteristicZero, SimpleDomain):
    """A class for arbitrary expressions. """

    is_EX = True

    class Expression(object):
        """An arbitrary expression. """

        __slots__ = ['ex']

        def __init__(self, ex):
            if not isinstance(ex, self.__class__):
                self.ex = sympify(ex)
            else:
                self.ex = ex.ex

        def __repr__(self):
            return 'EX(%s)' % repr(self.ex)

        def __str__(self):
            return 'EX(%s)' % str(self.ex)

        def __hash__(self):
            return hash((self.__class__.__name__, self.ex))

        def as_expr(self):
            return self.ex

        def numer(self):
            return EX(self.ex.as_numer_denom()[0])

        def denom(self):
            return EX(self.ex.as_numer_denom()[1])

        def simplify(self, ex):
            return self.__class__(ex.cancel())

        def __abs__(self):
            return self.__class__(abs(f.ex))

        def __neg__(self):
            return self.__class__(-self.ex)

        def __add__(self, g):
            return self.simplify(self.ex+self.__class__(g).ex)

        def __radd__(self, g):
            return self.simplify(self.__class__(g).ex+self.ex)

        def __sub__(self, g):
            return self.simplify(self.ex-self.__class__(g).ex)

        def __rsub__(self, g):
            return self.simplify(self.__class__(g).ex-self.ex)

        def __mul__(self, g):
            return self.simplify(self.ex*self.__class__(g).ex)

        def __rmul__(self, g):
            return self.simplify(self.__class__(g).ex*self.ex)

        def __pow__(self, n):
            return self.simplify(self.ex**n)

        def __div__(self, g):
            return self.simplify(self.ex/self.__class__(g).ex)

        def __rdiv__(self, g):
            return self.simplify(self.__class__(g).ex/self.ex)

        def __truediv__(self, g):
            return self.simplify(self.ex/self.__class__(g).ex)

        def __rtruediv__(self, g):
            return self.simplify(self.__class__(g).ex/self.ex)

        def __eq__(self, g):
            return self.ex == self.__class__(g).ex

        def __req__(self, g):
            return self.__class__(g).ex == self.ex

        def __ne__(self, g):
            return self.ex != self.__class__(g).ex

        def __rne__(self, g):
            return self.__class__(g).ex != self.ex

        def __nonzero__(self):
            return self.ex != 0

    dtype = Expression

    zero  = Expression(0)
    one   = Expression(1)

    rep   = 'EX'

    has_assoc_Ring         = False
    has_assoc_Field        = True

    def __init__(self):
        pass

    def to_sympy(self, a):
        """Convert `a` to a SymPy object. """
        return a.as_expr()

    def from_sympy(self, a):
        """Convert SymPy's expression to `dtype`. """
        return self.dtype(a)

    def from_ZZ_python(K1, a, K0):
        """Convert a Python `int` object to `dtype`. """
        return K1(K0.to_sympy(a))

    def from_QQ_python(K1, a, K0):
        """Convert a Python `Fraction` object to `dtype`. """
        return K1(K0.to_sympy(a))

    def from_ZZ_sympy(K1, a, K0):
        """Convert a SymPy `Integer` object to `dtype`. """
        return K1(K0.to_sympy(a))

    def from_QQ_sympy(K1, a, K0):
        """Convert a SymPy `Rational` object to `dtype`. """
        return K1(K0.to_sympy(a))

    def from_ZZ_gmpy(K1, a, K0):
        """Convert a GMPY `mpz` object to `dtype`. """
        return K1(K0.to_sympy(a))

    def from_QQ_gmpy(K1, a, K0):
        """Convert a GMPY `mpq` object to `dtype`. """
        return K1(K0.to_sympy(a))

    def from_RR_sympy(K1, a, K0):
        """Convert a SymPy `Real` object to `dtype`. """
        return K1(K0.to_sympy(a))

    def from_RR_mpmath(K1, a, K0):
        """Convert a mpmath `mpf` object to `dtype`. """
        return K1(K0.to_sympy(a))

    def from_PolynomialRing(K1, a, K0):
        """Convert a `DMP` object to `dtype`. """
        return K1(K0.to_sympy(a))

    def from_FractionField(K1, a, K0):
        """Convert a `DMF` object to `dtype`. """
        return K1(K0.to_sympy(a))

    def from_ExpressionDomain(K1, a, K0):
        """Convert a `EX` object to `dtype`. """
        return a

    def get_ring(self):
        """Returns a ring associated with `self`. """
        raise DomainError('there is no ring associated with %s' % self)

    def get_field(self):
        """Returns a field associated with `self`. """
        return self

    def is_positive(self, a):
        """Returns True if `a` is positive. """
        return a.ex.as_coeff_mul()[0].is_positive

    def is_negative(self, a):
        """Returns True if `a` is negative. """
        return a.ex.as_coeff_mul()[0].is_negative

    def is_nonpositive(self, a):
        """Returns True if `a` is non-positive. """
        return a.ex.as_coeff_mul()[0].is_nonpositive

    def is_nonnegative(self, a):
        """Returns True if `a` is non-negative. """
        return a.ex.as_coeff_mul()[0].is_nonnegative

    def numer(self, a):
        """Returns numerator of `a`. """
        return a.numer()

    def denom(self, a):
        """Returns denominator of `a`. """
        return a.denom()

"""Implementaton of :class:`GMPYRationalField` class. """

from sympy.polys.domains.rationalfield import RationalField

from sympy.polys.domains.groundtypes import (
    GMPYRationalType, SymPyRationalType,
    gmpy_numer, gmpy_denom, gmpy_factorial,
    gmpy_gcdex, gmpy_gcd, gmpy_lcm, gmpy_sqrt,
)

from sympy.polys.polyerrors import CoercionFailed

class GMPYRationalField(RationalField):
    """Rational field based on GMPY mpq class. """

    dtype = GMPYRationalType
    zero  = dtype(0)
    one   = dtype(1)
    alias = 'QQ_gmpy'

    def __init__(self):
        pass

    def to_sympy(self, a):
        """Convert `a` to a SymPy object. """
        return SymPyRationalType(int(gmpy_numer(a)),
                                 int(gmpy_denom(a)))

    def from_sympy(self, a):
        """Convert SymPy's Integer to `dtype`. """
        if a.is_Rational and a.q != 0:
            return GMPYRationalType(a.p, a.q)
        elif a.is_Real:
            from sympy.polys.domains import RR
            return GMPYRationalType(*RR.as_integer_ratio(a))
        else:
            raise CoercionFailed("expected `Rational` object, got %s" % a)

    def from_ZZ_python(K1, a, K0):
        """Convert a Python `int` object to `dtype`. """
        return GMPYRationalType(a)

    def from_QQ_python(K1, a, K0):
        """Convert a Python `Fraction` object to `dtype`. """
        return GMPYRationalType(a.numerator, a.denominator)

    def from_ZZ_sympy(K1, a, K0):
        """Convert a SymPy `Integer` object to `dtype`. """
        return GMPYRationalType(a.p)

    def from_QQ_sympy(K1, a, K0):
        """Convert a SymPy `Rational` object to `dtype`. """
        return GMPYRationalType(a.p, a.q)

    def from_ZZ_gmpy(K1, a, K0):
        """Convert a GMPY `mpz` object to `dtype`. """
        return GMPYRationalType(a)

    def from_QQ_gmpy(K1, a, K0):
        """Convert a GMPY `mpq` object to `dtype`. """
        return a

    def from_RR_sympy(K1, a, K0):
        """Convert a SymPy `Real` object to `dtype`. """
        return GMPYRationalType(*K0.as_integer_ratio(a))

    def from_RR_mpmath(K1, a, K0):
        """Convert a mpmath `mpf` object to `dtype`. """
        return GMPYRationalType(*K0.as_integer_ratio(a))

    def exquo(self, a, b):
        """Exact quotient of `a` and `b`, implies `__div__`.  """
        return GMPYRationalType(a.qdiv(b))

    def quo(self, a, b):
        """Quotient of `a` and `b`, implies `__div__`. """
        return GMPYRationalType(a.qdiv(b))

    def rem(self, a, b):
        """Remainder of `a` and `b`, implies nothing.  """
        return self.zero

    def div(self, a, b):
        """Division of `a` and `b`, implies `__div__`. """
        return GMPYRationalType(a.qdiv(b)), self.zero

    def numer(self, a):
        """Returns numerator of `a`. """
        return gmpy_numer(a)

    def denom(self, a):
        """Returns denominator of `a`. """
        return gmpy_denom(a)

    def factorial(self, a):
        """Returns factorial of `a`. """
        return GMPYRationalType(gmpy_factorial(int(a)))


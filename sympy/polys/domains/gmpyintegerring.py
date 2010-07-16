"""Implementaton of :class:`GMPYIntegerRing` class. """

from sympy.polys.domains.integerring import IntegerRing

from sympy.polys.domains.groundtypes import (
    GMPYIntegerType, SymPyIntegerType,
    gmpy_numer, gmpy_denom, gmpy_factorial,
    gmpy_gcdex, gmpy_gcd, gmpy_lcm, gmpy_sqrt,
)

from sympy.polys.polyerrors import CoercionFailed

class GMPYIntegerRing(IntegerRing):
    """Integer ring based on GMPY mpz class. """

    dtype = GMPYIntegerType
    zero  = dtype(0)
    one   = dtype(1)
    alias = 'ZZ_gmpy'

    def __init__(self):
        pass

    def to_sympy(self, a):
        """Convert `a` to a SymPy object. """
        return SymPyIntegerType(int(a))

    def from_sympy(self, a):
        """Convert SymPy's Integer to `dtype`. """
        if a.is_Integer:
            return GMPYIntegerType(a.p)
        elif a.is_Real and int(a) == a:
            return GMPYIntegerType(int(a))
        else:
            raise CoercionFailed("expected Integer object, got %s" % a)

    def from_ZZ_python(K1, a, K0):
        """Convert a Python `int` object to `dtype`. """
        return GMPYIntegerType(a)

    def from_QQ_python(K1, a, K0):
        """Convert a Python `Fraction` object to `dtype`. """
        if a.denominator == 1:
            return GMPYIntegerType(a.numerator)

    def from_ZZ_sympy(K1, a, K0):
        """Convert a SymPy `Integer` object to `dtype`. """
        return GMPYIntegerType(a.p)

    def from_QQ_sympy(K1, a, K0):
        """Convert a SymPy `Rational` object to `dtype`. """
        if a.q == 1:
            return GMPYIntegerType(a.p)

    def from_ZZ_gmpy(K1, a, K0):
        """Convert a GMPY `mpz` object to `dtype`. """
        return a

    def from_QQ_gmpy(K1, a, K0):
        """Convert a GMPY `mpq` object to `dtype`. """
        if a.denom() == 1:
            return a.numer()

    def from_RR_sympy(K1, a, K0):
        """Convert a SymPy `Real` object to `dtype`. """
        p, q = K0.as_integer_ratio(a)

        if q == 1:
            return GMPYIntegerType(p)

    def from_RR_mpmath(K1, a, K0):
        """Convert a mpmath `mpf` object to `dtype`. """
        p, q = K0.as_integer_ratio(a)

        if q == 1:
            return GMPYIntegerType(p)

    def gcdex(self, a, b):
        """Extended GCD of `a` and `b`. """
        h, s, t = gmpy_gcdex(a, b)
        return s, t, h

    def gcd(self, a, b):
        """Returns GCD of `a` and `b`. """
        return gmpy_gcd(a, b)

    def lcm(self, a, b):
        """Returns LCM of `a` and `b`. """
        return gmpy_lcm(a, b)

    def sqrt(self, a):
        """Returns square root of `a`. """
        return gmpy_sqrt(a)

    def factorial(self, a):
        """Returns factorial of `a`. """
        return gmpy_factorial(int(a))

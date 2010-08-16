"""Implementaton of :class:`GMPYIntegerRing` class. """

from sympy.polys.domains.integerring import IntegerRing

from sympy.polys.domains.groundtypes import (
    GMPYIntegerType, SymPyIntegerType,
    gmpy_numer, gmpy_denom, gmpy_factorial,
    gmpy_gcdex, gmpy_gcd, gmpy_lcm, gmpy_sqrt,
)

from sympy.polys.polyerrors import CoercionFailed

class GMPYIntegerRing(IntegerRing):
    """Integer ring based on GMPY's `mpz` type. """

    dtype = GMPYIntegerType
    zero  = dtype(0)
    one   = dtype(1)
    alias = 'ZZ_gmpy'

    def __init__(self):
        """Allow instantiation of this domain. """

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
            raise CoercionFailed("expected an integer, got %s" % a)

    def from_FF_python(K1, a, K0):
        """Convert `ModularInteger(int)` to GMPY's `mpz`. """
        return GMPYIntegerType(a.to_int())

    def from_ZZ_python(K1, a, K0):
        """Convert Python's `int` to GMPY's `mpz`. """
        return GMPYIntegerType(a)

    def from_QQ_python(K1, a, K0):
        """Convert Python's `Fraction` to GMPY's `mpz`. """
        if a.denominator == 1:
            return GMPYIntegerType(a.numerator)

    def from_FF_sympy(K1, a, K0):
        """Convert `ModularInteger(Integer)` to GMPY's `mpz`. """
        return GMPYIntegerType(a.to_int().p)

    def from_ZZ_sympy(K1, a, K0):
        """Convert SymPy's `Integer` to GMPY's `mpz`. """
        return GMPYIntegerType(a.p)

    def from_QQ_sympy(K1, a, K0):
        """Convert SymPy's `Rational` to GMPY's `mpz`. """
        if a.q == 1:
            return GMPYIntegerType(a.p)

    def from_FF_gmpy(K1, a, K0):
        """Convert `ModularInteger(mpz)` to GMPY's `mpz`. """
        return a.to_int()

    def from_ZZ_gmpy(K1, a, K0):
        """Convert GMPY's `mpz` to GMPY's `mpz`. """
        return a

    def from_QQ_gmpy(K1, a, K0):
        """Convert GMPY `mpq` to GMPY's `mpz`. """
        if a.denom() == 1:
            return a.numer()

    def from_RR_sympy(K1, a, K0):
        """Convert SymPy's `Real` to GMPY's `mpz`. """
        p, q = K0.as_integer_ratio(a)

        if q == 1:
            return GMPYIntegerType(p)

    def from_RR_mpmath(K1, a, K0):
        """Convert mpmath's `mpf` to GMPY's `mpz`. """
        p, q = K0.as_integer_ratio(a)

        if q == 1:
            return GMPYIntegerType(p)

    def gcdex(self, a, b):
        """Compute extended GCD of `a` and `b`. """
        h, s, t = gmpy_gcdex(a, b)
        return s, t, h

    def gcd(self, a, b):
        """Compute GCD of `a` and `b`. """
        return gmpy_gcd(a, b)

    def lcm(self, a, b):
        """Compute LCM of `a` and `b`. """
        return gmpy_lcm(a, b)

    def sqrt(self, a):
        """Compute square root of `a`. """
        return gmpy_sqrt(a)

    def factorial(self, a):
        """Compute factorial of `a`. """
        return gmpy_factorial(a)


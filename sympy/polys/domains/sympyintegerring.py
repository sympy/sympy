"""Implementaton of :class:`SymPyIntegerRing` class. """

from sympy.polys.domains.integerring import IntegerRing
from sympy.polys.domains.groundtypes import SymPyIntegerType

from sympy.polys.polyerrors import CoercionFailed

class SymPyIntegerRing(IntegerRing):
    """Integer ring based on SymPy's `Integer` type. """

    dtype = SymPyIntegerType
    zero  = dtype(0)
    one   = dtype(1)
    alias = 'ZZ_sympy'

    def __init__(self):
        """Allow instantiation of this domain. """

    def to_sympy(self, a):
        """Convert `a` to a SymPy object. """
        return a

    def from_sympy(self, a):
        """Convert SymPy's Integer to SymPy's `Integer`. """
        if a.is_Integer:
            return a
        elif a.is_Real and int(a) == a:
            return SymPyIntegerType(int(a))
        else:
            raise CoercionFailed("expected an integer, got %s" % a)

    def from_FF_python(K1, a, K0):
        """Convert `ModularInteger(int)` to SymPy's `Integer`. """
        return SymPyIntegerType(a.to_int())

    def from_ZZ_python(K1, a, K0):
        """Convert Python's `int` to SymPy's `Integer`. """
        return SymPyIntegerType(a)

    def from_QQ_python(K1, a, K0):
        """Convert Python's `Fraction` to SymPy's `Integer`. """
        if a.denominator == 1:
            return SymPyIntegerType(a.numerator)

    def from_FF_sympy(K1, a, K0):
        """Convert `ModularInteger(Integer)` to SymPy's `Integer`. """
        return a.to_int()

    def from_ZZ_sympy(K1, a, K0):
        """Convert SymPy's `Integer` to SymPy's `Integer`. """
        return a

    def from_QQ_sympy(K1, a, K0):
        """Convert SymPy's `Rational` to SymPy's `Integer`. """
        if a.q == 1:
            return SymPyIntegerType(a.p)

    def from_FF_gmpy(K1, a, K0):
        """Convert `ModularInteger(mpz)` to SymPy's `Integer`. """
        return SymPyIntegerType(int(a.to_int()))

    def from_ZZ_gmpy(K1, a, K0):
        """Convert GMPY's `mpz` to SymPy's `Integer`. """
        return SymPyIntegerType(int(a))

    def from_QQ_gmpy(K1, a, K0):
        """Convert GMPY's `mpq` to SymPy's `Integer`. """
        if a.denom() == 1:
            return SymPyIntegerType(int(a.numer()))

    def from_RR_sympy(K1, a, K0):
        """Convert SymPy's `Real` to SymPy's `Integer`. """
        p, q = K0.as_integer_ratio(a)

        if q == 1:
            return SymPyIntegerType(p)

    def from_RR_mpmath(K1, a, K0):
        """Convert mpmath's `mpf` to SymPy's `Integer`. """
        p, q = K0.as_integer_ratio(a)

        if q == 1:
            return SymPyIntegerType(p)

    def gcdex(self, a, b):
        """Compute extended GCD of `a` and `b`. """
        return a.gcdex(b)

    def gcd(self, a, b):
        """Compute GCD of `a` and `b`. """
        return a.gcd(b)

    def lcm(self, a, b):
        """Compute LCM of `a` and `b`. """
        return a.lcm(b)

    def sqrt(self, a):
        """Compute square root of `a`. """
        return a.sqrt()

    def factorial(self, a):
        """Compute factorial of `a`. """
        return a.factorial()

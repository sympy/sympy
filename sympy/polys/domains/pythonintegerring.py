"""Implementaton of :class:`PythonIntegerRing` class. """

from sympy.polys.domains.integerring import IntegerRing

from sympy.polys.domains.groundtypes import (
    PythonIntegerType, SymPyIntegerType, python_sqrt,
    python_factorial, python_gcdex, python_gcd, python_lcm,
)

from sympy.polys.polyerrors import CoercionFailed

class PythonIntegerRing(IntegerRing):
    """Integer ring based on Python int class. """

    dtype = PythonIntegerType
    zero  = dtype(0)
    one   = dtype(1)
    alias = 'ZZ_python'

    def __init__(self):
        pass

    def to_sympy(self, a):
        """Convert `a` to a SymPy object. """
        return SymPyIntegerType(a)

    def from_sympy(self, a):
        """Convert SymPy's Integer to `dtype`. """
        if a.is_Integer:
            return PythonIntegerType(a.p)
        elif a.is_Real and int(a) == a:
            return SymPyIntegerType(int(a))
        else:
            raise CoercionFailed("expected `Integer` object, got %s" % a)

    def from_ZZ_python(K1, a, K0):
        """Convert a Python `int` object to `dtype`. """
        return a

    def from_QQ_python(K1, a, K0):
        """Convert a Python `Fraction` object to `dtype`. """
        if a.denominator == 1:
            return a.numerator

    def from_ZZ_sympy(K1, a, K0):
        """Convert a SymPy `Integer` object to `dtype`. """
        return a.p

    def from_QQ_sympy(K1, a, K0):
        """Convert a SymPy `Rational` object to `dtype`. """
        if a.q == 1:
            return a.p

    def from_ZZ_gmpy(K1, a, K0):
        """Convert a GMPY `mpz` object to `dtype`. """
        return PythonIntegerType(a)

    def from_QQ_gmpy(K1, a, K0):
        """Convert a GMPY `mpq` object to `dtype`. """
        if a.denom() == 1:
            return PythonIntegerType(a.numer())

    def from_RR_sympy(K1, a, K0):
        """Convert a SymPy `Real` object to `dtype`. """
        p, q = K0.as_integer_ratio(a)

        if q == 1:
            return PythonIntegerType(p)

    def from_RR_mpmath(K1, a, K0):
        """Convert a mpmath `mpf` object to `dtype`. """
        p, q = K0.as_integer_ratio(a)

        if q == 1:
            return PythonIntegerType(p)

    def gcdex(self, a, b):
        """Extended GCD of `a` and `b`. """
        return python_gcdex(a, b)

    def gcd(self, a, b):
        """Returns GCD of `a` and `b`. """
        return python_gcd(a, b)

    def lcm(self, a, b):
        """Returns LCM of `a` and `b`. """
        return python_lcm(a, b)

    def sqrt(self, a):
        """Returns square root of `a`. """
        return PythonIntegerType(python_sqrt(a))

"""Implementaton of :class:`PythonRationalField` class. """

from sympy.polys.domains.rationalfield import RationalField

from sympy.polys.domains.groundtypes import PythonIntegerType
from sympy.polys.domains.groundtypes import PythonRationalType
from sympy.polys.domains.groundtypes import SymPyRationalType

from sympy.polys.polyerrors import CoercionFailed

class PythonRationalField(RationalField):
    """Rational field based on Python Fraction class. """

    dtype = PythonRationalType
    zero  = dtype(0)
    one   = dtype(1)
    alias = 'QQ_python'

    def __init__(self):
        pass

    def to_sympy(self, a):
        """Convert `a` to a SymPy object. """
        return SymPyRationalType(a.numerator, a.denominator)

    def from_sympy(self, a):
        """Convert SymPy's Rational to `dtype`. """
        if a.is_Rational and a.q != 0:
            return PythonRationalType(a.p, a.q)
        elif a.is_Real:
            from sympy.polys.domains import RR
            return PythonRationalType(*RR.as_integer_ratio(a))
        else:
            raise CoercionFailed("expected `Rational` object, got %s" % a)

    def from_ZZ_python(K1, a, K0):
        """Convert a Python `int` object to `dtype`. """
        return PythonRationalType(a)

    def from_QQ_python(K1, a, K0):
        """Convert a Python `Fraction` object to `dtype`. """
        return a

    def from_ZZ_sympy(K1, a, K0):
        """Convert a SymPy `Integer` object to `dtype`. """
        return PythonRationalType(a.p)

    def from_QQ_sympy(K1, a, K0):
        """Convert a SymPy `Rational` object to `dtype`. """
        return PythonRationalType(a.p, a.q)

    def from_ZZ_gmpy(K1, a, K0):
        """Convert a GMPY `mpz` object to `dtype`. """
        return PythonRationalType(PythonIntegerType(a))

    def from_QQ_gmpy(K1, a, K0):
        """Convert a GMPY `mpq` object to `dtype`. """
        return PythonRationalType(PythonIntegerType(a.numer()),
                                  PythonIntegerType(a.denom()))

    def from_RR_sympy(K1, a, K0):
        """Convert a SymPy `Real` object to `dtype`. """
        return PythonRationalType(*K0.as_integer_ratio(a))

    def from_RR_mpmath(K1, a, K0):
        """Convert a mpmath `mpf` object to `dtype`. """
        return PythonRationalType(*K0.as_integer_ratio(a))

    def numer(self, a):
        """Returns numerator of `a`. """
        return a.numerator

    def denom(self, a):
        """Returns denominator of `a`. """
        return a.denominator


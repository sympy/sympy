"""Implementation of :class:`PythonRationalField` class. """

from sympy.polys.domains.rationalfield import RationalField

from sympy.polys.domains.groundtypes import PythonInteger
from sympy.polys.domains.groundtypes import PythonRational
from sympy.polys.domains.groundtypes import SymPyRational

from sympy.polys.polyerrors import CoercionFailed

class PythonRationalField(RationalField):
    """Rational field based on Python rational number type. """

    dtype = PythonRational
    zero = dtype(0)
    one = dtype(1)
    alias = 'QQ_python'

    def __init__(self):
        pass

    def get_ring(self):
        """Returns ring associated with ``self``. """
        from sympy.polys.domains import PythonIntegerRing
        return PythonIntegerRing()

    def to_sympy(self, a):
        """Convert `a` to a SymPy object. """
        return SymPyRational(a.numerator, a.denominator)

    def from_sympy(self, a):
        """Convert SymPy's Rational to `dtype`. """
        if a.is_Rational:
            return PythonRational(a.p, a.q)
        elif a.is_Float:
            from sympy.polys.domains import RR
            return PythonRational(*RR.as_integer_ratio(a))
        else:
            raise CoercionFailed("expected `Rational` object, got %s" % a)

    def from_ZZ_python(K1, a, K0):
        """Convert a Python `int` object to `dtype`. """
        return PythonRational(a)

    def from_QQ_python(K1, a, K0):
        """Convert a Python `Fraction` object to `dtype`. """
        return a

    def from_ZZ_gmpy(K1, a, K0):
        """Convert a GMPY `mpz` object to `dtype`. """
        return PythonRational(PythonInteger(a))

    def from_QQ_gmpy(K1, a, K0):
        """Convert a GMPY `mpq` object to `dtype`. """
        return PythonRational(PythonInteger(a.numer()),
                              PythonInteger(a.denom()))

    def from_RR_mpmath(K1, a, K0):
        """Convert a mpmath `mpf` object to `dtype`. """
        return PythonRational(*K0.as_integer_ratio(a))

    def numer(self, a):
        """Returns numerator of `a`. """
        return a.numerator

    def denom(self, a):
        """Returns denominator of `a`. """
        return a.denominator

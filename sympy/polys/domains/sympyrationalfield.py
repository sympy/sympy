"""Implementaton of :class:`SymPyRationalField` class. """

from sympy.polys.domains.rationalfield import RationalField
from sympy.polys.domains.groundtypes import SymPyIntegerType
from sympy.polys.domains.groundtypes import SymPyRationalType

from sympy.polys.polyerrors import CoercionFailed

from sympy import (
    Rational as sympy_rat,
)

class SymPyRationalField(RationalField):
    """Rational field based on SymPy Rational class. """

    dtype = SymPyRationalType
    zero  = dtype(0)
    one   = dtype(1)
    alias = 'QQ_sympy'

    def __init__(self):
        pass

    def of_type(self, a):
        """
        Check if ``a`` is of type ``Rational``.

        Examples
        =======

        >>> from sympy import Rational, Real
        >>> from sympy.polys.domains import QQ_sympy
        >>> QQ_sympy().of_type(Rational(3, 2))
        True
        >>> QQ_sympy().of_type(2)
        False
        """
        return type(a) in [type(self.one), type(self.zero), type(sympy_rat(-1)),
                           type(sympy_rat(2)), type(sympy_rat(1, 2)),
                           type(sympy_rat(3, 2))]

    def to_sympy(self, a):
        """Convert `a` to a SymPy object. """
        return a

    def from_sympy(self, a):
        """Convert SymPy's Rational to `dtype`. """
        if a.is_Rational:
            return a
        elif a.is_Float:
            from sympy.polys.domains import RR
            return SymPyRationalType(*RR.as_integer_ratio(a))
        else:
            raise CoercionFailed("expected `Rational` object, got %s" % a)

    def from_ZZ_python(K1, a, K0):
        """Convert a Python `int` object to `dtype`. """
        return SymPyRationalType(a)

    def from_QQ_python(K1, a, K0):
        """Convert a Python `Fraction` object to `dtype`. """
        return SymPyRationalType(a.numerator, a.denominator)

    def from_ZZ_sympy(K1, a, K0):
        """Convert a SymPy `Integer` object to `dtype`. """
        return SymPyRationalType(a.p)

    def from_QQ_sympy(K1, a, K0):
        """Convert a SymPy `Rational` object to `dtype`. """
        return a

    def from_ZZ_gmpy(K1, a, K0):
        """Convert a GMPY `mpz` object to `dtype`. """
        return SymPyRationalType(int(a))

    def from_QQ_gmpy(K1, a, K0):
        """Convert a GMPY `mpq` object to `dtype`. """
        return SymPyRationalType(int(a.numer()),
                                 int(a.denom()))

    def from_RR_sympy(K1, a, K0):
        """Convert a SymPy `Float` object to `dtype`. """
        return SymPyRationalType(*K0.as_integer_ratio(a))

    def from_RR_mpmath(K1, a, K0):
        """Convert a mpmath `mpf` object to `dtype`. """
        return SymPyRationalType(*K0.as_integer_ratio(a))

    def numer(self, a):
        """Returns numerator of `a`. """
        return SymPyIntegerType(a.p)

    def denom(self, a):
        """Returns denominator of `a`. """
        return SymPyIntegerType(a.q)

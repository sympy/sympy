"""Implementation of :class:`SymPyRealDomain` class. """

from sympy.polys.domains.realdomain import RealDomain
from sympy.polys.domains.groundtypes import SymPyRealType

from sympy.core import S

class SymPyRealDomain(RealDomain):
    """Domain for real numbers based on SymPy Real type. """

    dtype = SymPyRealType
    zero  = dtype(0)
    one   = dtype(1)
    alias = 'RR_sympy'

    def __init__(self):
        pass

    def to_sympy(self, a):
        """Convert `a` to a SymPy object. """
        return a

    def from_sympy(self, a):
        """Convert SymPy's Integer to `dtype`. """
        b = a.evalf()

        if b.is_Real and b not in [S.Infinity, S.NegativeInfinity]:
            return b
        else:
            raise CoercionFailed("expected Real object, got %s" % a)

    def from_ZZ_python(K1, a, K0):
        """Convert a Python `int` object to `dtype`. """
        return SymPyRealType(a)

    def from_QQ_python(K1, a, K0):
        """Convert a Python `Fraction` object to `dtype`. """
        return SymPyRealType(a.numerator) / a.denominator

    def from_ZZ_sympy(K1, a, K0):
        """Convert a SymPy `Integer` object to `dtype`. """
        return SymPyRealType(a.p)

    def from_QQ_sympy(K1, a, K0):
        """Convert a SymPy `Rational` object to `dtype`. """
        return SymPyRealType(a.p) / a.q

    def from_ZZ_gmpy(K1, a, K0):
        """Convert a GMPY `mpz` object to `dtype`. """
        return SymPyRealType(int(a))

    def from_QQ_gmpy(K1, a, K0):
        """Convert a GMPY `mpq` object to `dtype`. """
        return SymPyRealType(int(a.numer())) / int(a.denom())

    def from_RR_sympy(K1, a, K0):
        """Convert a SymPy `Real` object to `dtype`. """
        return a

    def from_RR_mpmath(K1, a, K0):
        """Convert a mpmath `mpf` object to `dtype`. """
        return SymPyRealType(a)


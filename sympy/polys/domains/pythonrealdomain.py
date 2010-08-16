"""Implementation of :class:`PythonRealDomain` class. """

from sympy.polys.domains.realdomain import RealDomain
from sympy.polys.polyerrors import CoercionFailed

from sympy.core import S

class PythonRealDomain(RealDomain): # XXX: tmp solution
    """Float domain. """

    rep   = 'RR'

    is_FF = True

    dtype = float
    zero  = dtype(0)
    one   = dtype(1)
    alias = 'RR_python'

    def __init__(self):
        pass

    def normal(self, a):
        if abs(a) < 1e-15:
            return self.zero
        else:
            return self.dtype(a)

    def to_sympy(self, a):
        """Convert `a` to a SymPy object. """
        return sympy_mpf(a)

    def from_sympy(self, a):
        """Convert SymPy's Integer to `dtype`. """
        b = a.evalf()

        if b.is_Real and b not in [S.Infinity, S.NegativeInfinity]:
            return float(b)
        else:
            raise CoercionFailed("expected Real object, got %s" % a)

    def from_ZZ_python(K1, a, K0):
        """Convert a Python `int` object to `dtype`. """
        return K1.dtype(a)

    def from_QQ_python(K1, a, K0):
        """Convert a Python `Fraction` object to `dtype`. """
        return K1.dtype(a.numerator) / a.denominator

    def from_ZZ_sympy(K1, a, K0):
        """Convert a SymPy `Integer` object to `dtype`. """
        return K1.dtype(a.p)

    def from_QQ_sympy(K1, a, K0):
        """Convert a SymPy `Rational` object to `dtype`. """
        return K1.dtype(a.p) / a.q

    def from_ZZ_gmpy(K1, a, K0):
        """Convert a GMPY `mpz` object to `dtype`. """
        return K1.dtype(int(a))

    def from_QQ_gmpy(K1, a, K0):
        """Convert a GMPY `mpq` object to `dtype`. """
        return K1.dtype(int(a.numer())) / int(a.denom)

    def from_RR_sympy(K1, a, K0):
        """Convert a SymPy `Real` object to `dtype`. """
        return K1.dtype(a)

    def from_RR_mpmath(K1, a, K0):
        """Convert a mpmath `mpf` object to `dtype`. """
        return K1.dtype(a)

    def complex_domain(self):
        return CC


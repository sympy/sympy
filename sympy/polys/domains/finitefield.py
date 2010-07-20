"""Implementation of :class:`FiniteField` class. """

from sympy.polys.domains.field import Field
from sympy.polys.domains.simpledomain import SimpleDomain
from sympy.polys.domains.groundtypes import SymPyIntegerType
from sympy.polys.domains.modularinteger import ModularIntegerFactory

from sympy.polys.polyerrors import CoercionFailed

class FiniteField(Field, SimpleDomain):
    """General class for finite fields. """

    rep = 'FF'

    is_Numerical = True

    has_assoc_Ring         = False
    has_assoc_Field        = True

    dom = None
    mod = None

    def __init__(self, mod, symmetric=True):
        if mod <= 0:
            raise ValueError('modulus must be a positive integer, got %s' % mod)

        self.dtype = ModularIntegerFactory(mod, self.dom, symmetric)
        self.zero  = self.dtype(0)
        self.one   = self.dtype(1)
        self.mod   = mod

    def __str__(self):
        return 'GF(%s)' % self.mod

    def __hash__(self):
        return hash((self.__class__.__name__, self.dtype, self.mod, self.dom))

    def __eq__(self, other):
        """Returns `True` if two domains are equivalent. """
        return self.dtype.__class__.__name__ == other.dtype.__class__.__name__ and \
                self.mod == other.mod and self.dom == other.dom # XXX: this is a hack

    def __ne__(self, other):
        """Returns `False` if two domains are equivalent. """
        return self.dtype.__class__.__name__ != other.dtype.__class__.__name__ or \
                self.mod != other.mod or self.dom != other.dom # XXX: this is a hack

    def characteristic(self):
        """Return the characteristic of this domain. """
        return self.mod

    def get_ring(self):
        """Returns a ring associated with `self`. """
        return None

    def get_field(self):
        """Returns a field associated with `self`. """
        return self

    def to_sympy(self, a):
        """Convert `a` to a SymPy object. """
        return SymPyIntegerType(int(a))

    def from_sympy(self, a):
        """Convert SymPy's Integer to SymPy's `Integer`. """
        if a.is_Integer:
            return self.dtype(self.dom.dtype(int(a)))
        elif a.is_Real and int(a) == a:
            return self.dtype(self.dom.dtype(int(a)))
        else:
            raise CoercionFailed("expected an integer, got %s" % a)

    def from_FF_python(K1, a, K0=None):
        """Convert `ModularInteger(int)` to `dtype`. """
        return K1.dtype(K1.dom.from_ZZ_python(a.val, K0.dom))

    def from_ZZ_python(K1, a, K0=None):
        """Convert Python's `int` to `dtype`. """
        return K1.dtype(K1.dom.from_ZZ_python(a, K0))

    def from_QQ_python(K1, a, K0=None):
        """Convert Python's `Fraction` to `dtype`. """
        if a.denominator == 1:
            return K1.from_ZZ_python(a.numerator)

    def from_FF_sympy(K1, a, K0=None):
        """Convert `ModularInteger(Integer)` to `dtype`. """
        return K1.dtype(K1.dom.from_ZZ_sympy(a.val, K0.dom))

    def from_ZZ_sympy(K1, a, K0=None):
        """Convert SymPy's `Integer` to `dtype`. """
        return K1.dtype(K1.dom.from_ZZ_sympy(a, K0))

    def from_QQ_sympy(K1, a, K0=None):
        """Convert SymPy's `Rational` to `dtype`. """
        if a.q == 1:
            return K1.from_ZZ_python(a.p)

    def from_FF_gmpy(K1, a, K0=None):
        """Convert `ModularInteger(mpz)` to `dtype`. """
        return K1.dtype(K1.dom.from_ZZ_gmpy(a.val, K0.dom))

    def from_ZZ_gmpy(K1, a, K0=None):
        """Convert GMPY's `mpz` to `dtype`. """
        return K1.dtype(K1.dom.from_ZZ_gmpy(a, K0))

    def from_QQ_gmpy(K1, a, K0=None):
        """Convert GMPY's `mpq` to `dtype`. """
        if a.denom() == 1:
            return K1.from_ZZ_gmpy(a.numer())

    def from_RR_sympy(K1, a, K0=None):
        """Convert SymPy's `Real` to `dtype`. """
        p, q = K0.as_integer_ratio(a)

        if q == 1:
            return K1.dtype(self.dom.dtype(p))

    def from_RR_mpmath(K1, a, K0):
        """Convert mpmath's `mpf` to `dtype`. """
        p, q = K0.as_integer_ratio(a)

        if q == 1:
            return K1.dtype(self.dom.dtype(p))

    def gcdex(self, a, b):
        """Compute extended GCD of `a` and `b`. """
        return self.dtype(a.gcdex(b))

    def gcd(self, a, b):
        """Compute GCD of `a` and `b`. """
        return self.dtype(a.gcd(b))

    def lcm(self, a, b):
        """Compute LCM of `a` and `b`. """
        return self.dtype(a.lcm(b))

    def sqrt(self, a):
        """Compute square root of `a`. """
        return self.dtype(a.sqrt())

    def factorial(self, a):
        """Compute factorial of `a`. """
        return self.dtype(self.dom.factorial(a))


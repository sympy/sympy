"""Implementation of :class:`PolynomialRing` class. """

from sympy.polys.domains.ring import Ring
from sympy.polys.domains.compositedomain import CompositeDomain

class PolynomialRingNG(Ring, CompositeDomain):
    """A class for representing multivariate polynomial rings. """

    is_Poly      = True

    has_assoc_Ring         = True
    has_assoc_Field        = True

    def __init__(self, ring):
        self.dtype = ring.dtype
        self.ring  = ring

        self.zero = ring.zero
        self.one  = ring.one

        self.dom  = ring.domain
        self.gens = ring.gens

    def __str__(self):
        return str(self.dom) + '[' + ','.join(map(str, self.gens)) + ']'

    def __hash__(self):
        return hash((self.__class__.__name__, self.dtype, self.dom, self.gens))

    def __eq__(self, other):
        """Returns `True` if two domains are equivalent. """
        return isinstance(other, PolynomialRingNG) and self.dtype == other.dtype and self.ring == other.ring

    def __ne__(self, other):
        """Returns `False` if two domains are equivalent. """
        return not self.__eq__(other)

    def to_sympy(self, a):
        """Convert `a` to a SymPy object. """
        return a.as_expr()

    def from_sympy(self, a):
        """Convert SymPy's expression to `dtype`. """
        return self.ring.ring_new(a)

    def __call__(self, a):
        """Construct an element of `self` domain from `a`. """
        return self.ring.ring_new(a)

    def from_ZZ_python(K1, a, K0):
        """Convert a Python `int` object to `dtype`. """
        return K1(K1.dom.convert(a, K0))

    def from_QQ_python(K1, a, K0):
        """Convert a Python `Fraction` object to `dtype`. """
        return K1(K1.dom.convert(a, K0))

    def from_ZZ_gmpy(K1, a, K0):
        """Convert a GMPY `mpz` object to `dtype`. """
        return K1(K1.dom.convert(a, K0))

    def from_QQ_gmpy(K1, a, K0):
        """Convert a GMPY `mpq` object to `dtype`. """
        return K1(K1.dom.convert(a, K0))

    def from_RR_mpmath(K1, a, K0):
        """Convert a mpmath `mpf` object to `dtype`. """
        return K1(K1.dom.convert(a, K0))

    def from_PolynomialRingNG(K1, a, K0):
        """Convert a `DMP` object to `dtype`. """
        if K1 == K0:
            return a
        else:
            return # TODO

    def from_FractionFieldNG(K1, a, K0):
        """Convert a `DMF` object to `dtype`. """
        if K1.ring == K0.field.ring and K0.denom(a) == 1:
            return K1.ring.ring_new(K0.numer(a))
        else:
            return # TODO

    def get_field(self):
        """Returns a field associated with `self`. """
        return self.ring.to_field().to_domain()

    def poly_ring(self, *symbols): # TODO:, order=lex):
        """Returns a polynomial ring, i.e. `K[X]`. """
        from sympy.polys.rings import PolyRing
        return PolyRing(symbols, self.ring, order).to_domain()

    def frac_field(self, *symbols): # TODO:, order=lex):
        """Returns a fraction field, i.e. `K(X)`. """
        from sympy.polys.fields import FracField
        return FracField(symbols, self.ring, order).to_domain()

    def is_positive(self, a):
        """Returns True if `LC(a)` is positive. """
        return self.dom.is_positive(a.LC)

    def is_negative(self, a):
        """Returns True if `LC(a)` is negative. """
        return self.dom.is_negative(a.LC)

    def is_nonpositive(self, a):
        """Returns True if `LC(a)` is non-positive. """
        return self.dom.is_nonpositive(a.LC)

    def is_nonnegative(self, a):
        """Returns True if `LC(a)` is non-negative. """
        return self.dom.is_nonnegative(a.LC)

    def gcdex(self, a, b):
        """Extended GCD of `a` and `b`. """
        return a.gcdex(b)

    def gcd(self, a, b):
        """Returns GCD of `a` and `b`. """
        return a.gcd(b)

    def lcm(self, a, b):
        """Returns LCM of `a` and `b`. """
        return a.lcm(b)

    def factorial(self, a):
        """Returns factorial of `a`. """
        return self.dtype(self.dom.factorial(a))

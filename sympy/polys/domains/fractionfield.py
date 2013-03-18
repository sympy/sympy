"""Implementation of :class:`PolynomialRing` class. """

from sympy.polys.domains.field import Field
from sympy.polys.domains.compositedomain import CompositeDomain

class FractionField(Field, CompositeDomain):
    """A class for representing multivariate rational function fields. """

    is_Frac      = True

    has_assoc_Ring         = True
    has_assoc_Field        = True

    def __init__(self, field):
        self.dtype = field.dtype
        self.field = field

        self.dom  = field.domain
        self.gens = field.symbols

        self.zero = field.zero
        self.one  = field.one


    def new(self, element):
        return self.field.field_new(element)

    def __str__(self):
        return str(self.dom) + '(' + ','.join(map(str, self.gens)) + ')'

    def __hash__(self):
        return hash((self.__class__.__name__, self.dtype, self.dom, self.gens))

    def __eq__(self, other):
        """Returns `True` if two domains are equivalent. """
        return isinstance(other, FractionField) and self.dtype == other.dtype and self.field == other.field

    def __ne__(self, other):
        """Returns `False` if two domains are equivalent. """
        return not self.__eq__(other)

    def to_sympy(self, a):
        """Convert `a` to a SymPy object. """
        return a.as_expr()

    def from_sympy(self, a):
        """Convert SymPy's expression to `dtype`. """
        return self.field.field_new(a)

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

    def from_PolynomialRing(K1, a, K0):
        """Convert a `DMP` object to `dtype`. """
        if K1.field.ring == K0.ring:
            return K1.field.field_new(a)
        else:
            return # TODO

    def from_FractionField(K1, a, K0):
        """Convert a `DMF` object to `dtype`. """
        if K1 == K0:
            return a
        else:
            return # TODO

    def get_ring(self):
        """Returns a field associated with `self`. """
        return self.field.to_ring().to_domain()

    def poly_ring(self, *symbols): # TODO:, order=lex):
        """Returns a polynomial ring, i.e. `K[X]`. """
        from sympy.polys.rings import PolyRing
        return PolyRing(symbols, self.field, order).to_domain()

    def frac_field(self, *symbols): # TODO:, order=lex):
        """Returns a fraction field, i.e. `K(X)`. """
        from sympy.polys.fields import FracField
        return FracField(symbols, self.field, order).to_domain()

    def is_positive(self, a):
        """Returns True if `LC(a)` is positive. """
        return self.dom.is_positive(a.numer.LC)

    def is_negative(self, a):
        """Returns True if `LC(a)` is negative. """
        return self.dom.is_negative(a.numer.LC)

    def is_nonpositive(self, a):
        """Returns True if `LC(a)` is non-positive. """
        return self.dom.is_nonpositive(a.numer.LC)

    def is_nonnegative(self, a):
        """Returns True if `LC(a)` is non-negative. """
        return self.dom.is_nonnegative(a.numer.LC)

    def numer(self, a):
        """Returns numerator of ``a``. """
        return a.numer

    def denom(self, a):
        """Returns denominator of ``a``. """
        return a.denom

    def factorial(self, a):
        """Returns factorial of `a`. """
        return self.dtype(self.dom.factorial(a))

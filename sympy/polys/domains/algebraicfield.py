"""Implementation of :class:`AlgebraicField` class. """

from __future__ import print_function, division

from sympy.core.sympify import _sympify
from sympy.polys.domains.characteristiczero import CharacteristicZero
from sympy.polys.domains.field import Field
from sympy.polys.domains.simpledomain import SimpleDomain
from sympy.polys.polyclasses import ANP
from sympy.polys.polyerrors import CoercionFailed, DomainError, NotAlgebraic, IsomorphismFailed
from sympy.utilities import public

@public
class AlgebraicField(Field, CharacteristicZero, SimpleDomain):
    """A class for representing algebraic number fields. """

    dtype = ANP

    is_AlgebraicField = is_Algebraic = True
    is_Numerical = True

    has_assoc_Ring = False
    has_assoc_Field = True

    def __new__(cls, dom, *ext):
        if not dom.is_QQ:
            raise DomainError("ground domain must be a rational field")

        from sympy.polys.numberfields import to_number_field

        obj = super(AlgebraicField, cls).__new__(cls)
        obj.orig_ext = ext
        obj.ext = to_number_field(ext)
        obj.mod = obj.ext.minpoly.rep
        obj.domain = obj.dom = dom

        obj.ngens = 1
        obj.symbols = obj.gens = (obj.ext,)
        obj.unit = obj([dom(1), dom(0)])

        obj.zero = obj.dtype.zero(obj.mod.rep, dom)
        obj.one = obj.dtype.one(obj.mod.rep, dom)
        return obj

    @property
    def args(self):
        # XXX Workaround for pickling
        return (self.dom, *self.orig_ext)

    def new(self, element):
        return self.dtype(element, self.mod.rep, self.dom)

    def __str__(self):
        return str(self.dom) + '<' + str(self.ext) + '>'

    def __hash__(self):
        return hash((self.__class__.__name__, self.dtype, self.dom, self.ext))

    def __eq__(self, other):
        """Returns ``True`` if two domains are equivalent. """
        return isinstance(other, AlgebraicField) and \
            self.dtype == other.dtype and self.ext == other.ext

    def algebraic_field(self, *extension):
        r"""Returns an algebraic field, i.e. `\mathbb{Q}(\alpha, \ldots)`. """
        return AlgebraicField(self.dom, *((self.ext,) + extension))

    def to_sympy(self, a):
        """Convert ``a`` to a SymPy object. """
        from sympy.polys.numberfields import AlgebraicNumber
        return AlgebraicNumber(self.ext, a).as_expr()

    def from_sympy(self, a):
        """Convert SymPy's expression to ``dtype``. """
        try:
            return self([self.dom.from_sympy(a)])
        except CoercionFailed:
            pass

        from sympy.polys.numberfields import to_number_field

        try:
            return self(to_number_field(a, self.ext).native_coeffs())
        except (NotAlgebraic, IsomorphismFailed):
            raise CoercionFailed(
                "%s is not a valid algebraic number in %s" % (a, self))

    def from_ZZ_python(K1, a, K0):
        """Convert a Python ``int`` object to ``dtype``. """
        return K1(K1.dom.convert(a, K0))

    def from_QQ_python(K1, a, K0):
        """Convert a Python ``Fraction`` object to ``dtype``. """
        return K1(K1.dom.convert(a, K0))

    def from_ZZ_gmpy(K1, a, K0):
        """Convert a GMPY ``mpz`` object to ``dtype``. """
        return K1(K1.dom.convert(a, K0))

    def from_QQ_gmpy(K1, a, K0):
        """Convert a GMPY ``mpq`` object to ``dtype``. """
        return K1(K1.dom.convert(a, K0))

    def from_RealField(K1, a, K0):
        """Convert a mpmath ``mpf`` object to ``dtype``. """
        return K1(K1.dom.convert(a, K0))

    def get_ring(self):
        """Returns a ring associated with ``self``. """
        raise DomainError('there is no ring associated with %s' % self)

    def is_positive(self, a):
        """Returns True if ``a`` is positive. """
        return self.dom.is_positive(a.LC())

    def is_negative(self, a):
        """Returns True if ``a`` is negative. """
        return self.dom.is_negative(a.LC())

    def is_nonpositive(self, a):
        """Returns True if ``a`` is non-positive. """
        return self.dom.is_nonpositive(a.LC())

    def is_nonnegative(self, a):
        """Returns True if ``a`` is non-negative. """
        return self.dom.is_nonnegative(a.LC())

    def numer(self, a):
        """Returns numerator of ``a``. """
        return a

    def denom(self, a):
        """Returns denominator of ``a``. """
        return self.one

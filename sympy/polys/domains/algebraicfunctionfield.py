"""Implementation of :class:`AlgebraicFunctionField` class. """

from sympy.polys.domains.field import Field
from sympy.polys.domains.simpledomain import SimpleDomain
from sympy.polys.domains.characteristiczero import CharacteristicZero

from sympy.polys.polyclasses import ANP
from sympy.polys.polyerrors import (
    CoercionFailed, DomainError, NotAlgebraic, IsomorphismFailed,
)


class AlgebraicFunctionField(Field, CharacteristicZero, SimpleDomain):
    """A class for representing algebraic function fields. """

    dtype = ANP

#    Should this be True?
#    is_Numerical = True
    is_Algebraic = True

    has_assoc_Ring = False
    has_assoc_Field = True

    def __init__(self, dom, *ext):
        from sympy import Dummy
        from sympy.polys.polytools import Poly
        from sympy.polys.functionfields import minpoly

        if not (dom.is_Frac and dom.dom.is_QQ):
            raise DomainError("ground domain must be a rational function field")

        if not ext:
            raise ValueError("empty extension")
        elif len(ext) > 1:
            raise NotImplementedError("only simple extensions are supported")
        self.ext = ext[0]

        Y = Dummy('Y')
        mp = Poly(minpoly(self.ext, Y, dom), Y, domain=dom)
        self.mod = mp.rep

        self.dom = dom

        self.gens = (self.ext,)
        self.unit = self([dom(1), dom(0)])

        self.zero = self.dtype.zero(self.mod.rep, dom)
        self.one = self.dtype.one(self.mod.rep, dom)

    def new(self, element):
        return self.dtype(element, self.mod.rep, self.dom)

    def __str__(self):
        return str(self.dom) + '<' + str(self.ext) + '>'

    def __hash__(self):
        return hash((self.__class__.__name__, self.dtype, self.dom, self.ext))

    def __eq__(self, other):
        """Returns ``True`` if two domains are equivalent. """
        return isinstance(other, AlgebraicFunctionField) and \
            self.dom == self.dom and self.ext == other.ext

    def algebraic_function_field(self, minpoly, *extension):
        r"""Returns an algebraic function field, i.e. `\mathbb{Q}(x)(\alpha, \dots)`. """
        return AlgebraicFunctionField(self.dom, *((self.ext,) + extension))

    def to_sympy(self, a):
        """Convert ``a`` to a SymPy object. """
        raise NotImplementedError

    def from_sympy(self, a):
        """Convert SymPy's expression to ``dtype``. """
        raise NotImplementedError

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

    def from_RR_mpmath(K1, a, K0):
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

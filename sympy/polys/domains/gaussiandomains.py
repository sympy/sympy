"""Domains of Gaussian type."""

from sympy.core.basic import Basic
from sympy.core.numbers import Rational, I
from sympy.polys.polyerrors import CoercionFailed
from sympy.polys.domains import ZZ, QQ
from sympy.polys.domains.algebraicfield import AlgebraicField
from sympy.polys.domains.domainelement import DomainElement
from sympy.polys.domains.field import Field
from sympy.polys.domains.ring import Ring


class GaussianElement(DomainElement):
    """Base class for elements of Gaussian type domains."""
    base = None  # base ring
    _parent = None

    __slots__ = ('x', 'y')

    def __init__(self, x, y=0):
        conv = self.base.convert
        self.x = conv(x)
        self.y = conv(y)

    def as_expr(self):
        return Rational(self.x) + I * Rational(self.y)

    def parent(self):
        return self._parent

    def __hash__(self):
        return hash((self.x, self.y))

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.x == other.x and self.y == other.y
        else:
            return NotImplemented

    def __neg__(self):
        return self.__class__(-self.x, -self.y)

    def __repr__(self):
        return "%s(%s, %s)" % (type(self).__name__, self.x, self.y)

    def __str__(self):
        return str(self._parent.to_sympy(self))

    @classmethod
    def _get_xy(cls, other):
        if not isinstance(other, cls):
            try:
                other = cls._parent.convert(other)
            except CoercionFailed:
                return None, None
        return other.x, other.y

    def __add__(self, other):
        x, y = self._get_xy(other)
        if x is not None:
            return self.__class__(self.x + x, self.y + y)
        else:
            return NotImplemented

    __radd__ = __add__

    def __sub__(self, other):
        x, y = self._get_xy(other)
        if x is not None:
            return self.__class__(self.x - x, self.y - y)
        else:
            return NotImplemented

    def __rsub__(self, other):
        x, y = self._get_xy(other)
        if x is not None:
            return self.__class__(x - self.x, y - self.y)
        else:
            return NotImplemented

    def __mul__(self, other):
        x, y = self._get_xy(other)
        if x is not None:
            return self.__class__(self.x*x - self.y*y,
                                  self.x*y + self.y*x)
        else:
            return NotImplemented

    __rmul__ = __mul__

    def __pow__(self, exp):
        if exp == 0:
            return self.__class__(1, 0)
        if exp < 0:
            self, exp = 1/self, -exp
        prod = self
        for n in range(exp-1):
            prod *= self
        return prod

    def __bool__(self):
        return bool(self.x) or bool(self.y)

    __nonzero__ = __bool__  # for Python 2

    def quadrant(self):
        """Return quadrant index 0-3.

        0 is included in quadrant 0.
        """
        if self.y > 0:
            return 0 if self.x > 0 else 1
        elif self.y < 0:
            return 2 if self.x < 0 else 3
        else:
            return 0 if self.x >= 0 else 2


class GaussianInteger(GaussianElement):
    base = ZZ

    def __truediv__(self, other):
        """Return a Gaussian rational."""
        return QQ_I.convert(self)/other

    __div__ = __truediv__

    def __rtruediv__(self, other):
        return other/QQ_I.convert(self)

    __rdiv__ = __rtruediv__

    def __divmod__(self, other):
        if not other:
            raise ZeroDivisionError('divmod({}, 0)'.format(self))
        x, y = self._get_xy(other)
        if x is None:
            return NotImplemented

        # multiply self and other by x - I*y
        # self/other == (a + I*b)/c
        a, b = self.x*x + self.y*y, -self.x*y + self.y*x
        c = x*x + y*y

        # find integers qx and qy such that
        # |a - qx*c| <= c/2 and |b - qy*c| <= c/2
        qx = (2*a + c) // (2*c)  # -c <= 2*a - qx*2*c < c
        qy = (2*b + c) // (2*c)

        q = GaussianInteger(qx, qy)
        # |self/other - q| < 1 since
        # |a/c - qx|**2 + |b/c - qy|**2 <= 1/4 + 1/4 < 1

        return q, self - q*other  # |r| < |other|


    def __rdivmod__(self, other):
        try:
            other = self._parent.convert(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return other.__divmod__(self)

    def __floordiv__(self, other):
        qr = self.__divmod__(other)
        return qr if qr is NotImplemented else qr[0]

    def __rfloordiv__(self, other):
        qr = self.__rdivmod__(other)
        return qr if qr is NotImplemented else qr[0]

    def __mod__(self, other):
        qr = self.__divmod__(other)
        return qr if qr is NotImplemented else qr[1]

    def __rmod__(self, other):
        qr = self.__rdivmod__(other)
        return qr if qr is NotImplemented else qr[1]


class GaussianRational(GaussianElement):
    base = QQ

    def __truediv__(self, other):
        if not other:
            raise ZeroDivisionError('{} / 0'.format(self))
        x, y = self._get_xy(other)
        if x is None:
            return NotImplemented
        c = x*x + y*y

        return GaussianRational((self.x*x + self.y*y)/c,
                                (-self.x*y + self.y*x)/c)

    __floordiv__ = __div__ = __truediv__

    def __rtruediv__(self, other):
        try:
            other = self._parent.convert(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return other.__truediv__(self)

    __rfloordiv__ = __rdiv__ = __rtruediv__

    def __mod__(self, other):
        try:
            other = self._parent.convert(other)
        except CoercionFailed:
            return NotImplemented
        if not other:
            raise ZeroDivisionError('{} % 0'.format(self))
        else:
            return self._parent.zero  # XXX always 0?

    def __rmod__(self, other):
        try:
            other = self._parent.convert(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return other.__mod__(self)

    def __divmod__(self, other):
        return self.__truediv__(other), self.__mod__(other)

    def __rdivmod__(self, other):
        return self.__rtruediv__(other), self.__rmod__(other)


class GaussianDomain():
    """Base class for Gaussian domains."""
    base = None  # base domain, ZZ or QQ

    has_assoc_Ring = True
    has_assoc_Field = True

    def to_sympy(self, a):
        """Convert ``a`` to a SymPy object. """
        conv = self.base.to_sympy
        return conv(a.x) + I*conv(a.y)

    def from_sympy(self, a):
        """Convert a SymPy object to ``self.dtype``."""
        r, b = a.as_coeff_Add()
        x = self.base.from_sympy(r)  # may raise CoercionFailed
        if not b:
            return self.new(x, 0)
        r, b = b.as_coeff_Mul()
        y = self.base.from_sympy(r)
        if b is I:
            return self.new(x, y)
        else:
            raise CoercionFailed("{} is not Gaussian".format(a))

    def convert(self, element, base=None):
        """Convert ``element`` to ``self.dtype``.

        Raises CoercionFailed on failure.
        """
        if isinstance(element, self.dtype):
            return element
        elif isinstance(element, GaussianElement):
            return self.new(element.x, element.y)
        elif isinstance(element, Basic):
            return self.from_sympy(element)
        else:  # convertible to base type or failure
            return self.new(element, 0)

    def inject(self, *gens):
        """Inject generators into this domain. """
        return self.poly_ring(*gens)

    def as_AlgebraicField(self):
        """Get equivalent domain as an ``AlgebraicField``. """
        return AlgebraicField(self.base, I)

    # Override the negative etc handlers because this isn't an ordered domain.

    def is_negative(self, element):
        """Returns ``False`` for any ``GaussianElement``. """
        return False

    def is_positive(self, element):
        """Returns ``False`` for any ``GaussianElement``. """
        return False

    def is_nonnegative(self, element):
        """Returns ``False`` for any ``GaussianElement``. """
        return False

    def is_nonpositive(self, element):
        """Returns ``False`` for any ``GaussianElement``. """
        return False


class GaussianIntegerRing(GaussianDomain, Ring):
    """Ring of Gaussian integers."""
    base = ZZ
    dtype = GaussianInteger
    zero = dtype(0, 0)
    one = dtype(1, 0)
    imag_unit = dtype(0, 1)
    units = (one, imag_unit, -one, -imag_unit)  # powers of i

    rep = 'ZZ_I'

    is_GaussianRing = True
    is_ZZ_I = True

    def __init__(self):  # override Domain.__init__
        """For constructing ZZ_I."""

    def get_ring(self):
        """Returns a ring associated with ``self``. """
        return self

    def get_field(self):
        """Returns a field associated with ``self``. """
        return QQ_I

    def normalize(self, d, *args):
        """Return first quadrant element associated with ``d``.

        Also multiply the other arguments by the same power of i.
        """
        unit = self.units[-d.quadrant()]  # - for inverse power
        d *= unit
        args = tuple(a*unit for a in args)
        return (d,) + args if args else d

    def gcd(self, a, b):
        while b:
            a, b = b, a % b
        return self.normalize(a)

    def lcm(self, a, b):
        return (a * b) // self.gcd(a, b)

ZZ_I = GaussianInteger._parent = GaussianIntegerRing()


class GaussianRationalField(GaussianDomain, Field):
    """Field of Gaussian rational numbers."""
    base = QQ
    dtype = GaussianRational
    zero = dtype(0, 0)
    one = dtype(1, 0)

    rep = 'QQ_I'

    is_GaussianField = True
    is_QQ_I = True

    def __init__(self):  # override Domain.__init__
        """For constructing QQ_I."""

    def get_ring(self):
        """Returns a ring associated with ``self``. """
        return ZZ_I

    def get_field(self):
        """Returns a field associated with ``self``. """
        return self

    def denom(self, a):
        """Get the denominator of ``a``."""
        ZZ = self.base.get_ring()
        QQ = self.base
        ZZ_I = self.get_ring()
        denom_ZZ = ZZ.lcm(QQ.denom(a.x), QQ.denom(a.y))
        return ZZ_I(denom_ZZ, ZZ.zero)

QQ_I = GaussianRational._parent = GaussianRationalField()

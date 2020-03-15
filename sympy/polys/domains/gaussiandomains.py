"""Domains of Gaussian type."""

from __future__ import print_function, division

from sympy.core.basic import Basic
from sympy.core.numbers import I
from sympy.polys.polyerrors import CoercionFailed
from sympy.polys.domains import ZZ, QQ
from sympy.polys.domains.domainelement import DomainElement
from sympy.polys.domains.field import Field
from sympy.polys.domains.ring import Ring


class GaussianElement(DomainElement):
    """Base class for elements of Gaussian type domains."""
    base = None  # base ring
    _parent = None

    __slots__ = ('x', 'y')

    def __init__(self, x, y):
        conv = self.base.convert
        self.x = conv(x)
        self.y = conv(y)

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

    __rmul__ = __mul__

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
            raise ZeroDivisionError('divmod(%s, 0)' % self)
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
        qr = self._rdivmod__(self, other)
        return qr if qr is NotImplemented else qr[0]

    def __mod__(self, other):
        qr = self.__divmod__(other)
        return qr if qr is NotImplemented else qr[1]

    def __rmod__(self, other):
        qr = self._rdivmod__(self, other)
        return qr if qr is NotImplemented else qr[1]


class GaussianRational(GaussianElement):
    base = QQ

    def __truediv__(self, other):
        if not other:
            raise ZeroDivisionError('%s / 0' % self)
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
            raise ZeroDivisionError('%s % 0' % self)
        else:
            return self._parent.zero

    def _rmod__(self, other):
        try:
            other = self._parent.convert(other)
        except CoercionFailed:
            return NotImplemented
        else:
            other.__mod__(self)

    def __divmod__(self, other):
        return self.__truediv__(other), self.__mod__(other)

    def _rdivmod__(self, other):
        return self._rtruediv__(other), self.__rmod__(other)


class GaussianDomain(object):
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
            raise CoercionFailed("%s is not Gaussian" % a)

    def convert(self, element):
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


class GaussianIntegerRing(GaussianDomain, Ring):
    """Ring of Gaussian integers."""
    base = ZZ
    dtype = GaussianInteger
    zero = dtype(0, 0)
    one = dtype(1, 0)
    imag_unit = dtype(0, 1)
    units = (one, imag_unit, -one, -imag_unit)  # powers of i

    rep = 'ZZ_I'

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

ZZ_I = GaussianInteger._parent = GaussianIntegerRing()


class GaussianRationalField(GaussianDomain, Field):
    """Field of Gaussian rational numbers."""
    base = QQ
    dtype = GaussianRational
    zero = dtype(0, 0)
    one = dtype(1, 0)

    rep = 'QQ_I'

    def __init__(self):  # override Domain.__init__
        """For constructing QQ_I."""

    def get_ring(self):
        """Returns a ring associated with ``self``. """
        return ZZ_I

    def get_field(self):
        """Returns a field associated with ``self``. """
        return self

QQ_I = GaussianRational._parent = GaussianRationalField()

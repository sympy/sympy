"""Domains of Gaussian type."""

from sympy.core.numbers import I
from sympy.polys.polyerrors import CoercionFailed
from sympy.polys.domains import ZZ, QQ
from sympy.polys.domains.algebraicfield import AlgebraicField
from sympy.polys.domains.domain import Domain
from sympy.polys.domains.domainelement import DomainElement
from sympy.polys.domains.field import Field
from sympy.polys.domains.ring import Ring


class GaussianElement(DomainElement):
    """Base class for elements of Gaussian type domains."""
    base = None  # type: Domain
    _parent = None  # type: Domain

    __slots__ = ('x', 'y')

    def __init__(self, x, y=0):
        conv = self.base.convert
        self.x = conv(x)
        self.y = conv(y)

    @classmethod
    def new(cls, x, y):
        """Create a new GaussianElement of the same domain."""
        return cls(x, y)

    def parent(self):
        """The domain that this is an element of (ZZ_I or QQ_I)"""
        return self._parent

    def __hash__(self):
        return hash((self.x, self.y))

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.x == other.x and self.y == other.y
        else:
            return NotImplemented

    def __lt__(self, other):
        if not isinstance(other, GaussianElement):
            return NotImplemented
        return [self.y, self.x] < [other.y, other.x]

    def __neg__(self):
        return self.new(-self.x, -self.y)

    def __repr__(self):
        return "%s(%s, %s)" % (self._parent.rep, self.x, self.y)

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
            return self.new(self.x + x, self.y + y)
        else:
            return NotImplemented

    __radd__ = __add__

    def __sub__(self, other):
        x, y = self._get_xy(other)
        if x is not None:
            return self.new(self.x - x, self.y - y)
        else:
            return NotImplemented

    def __rsub__(self, other):
        x, y = self._get_xy(other)
        if x is not None:
            return self.new(x - self.x, y - self.y)
        else:
            return NotImplemented

    def __mul__(self, other):
        x, y = self._get_xy(other)
        if x is not None:
            return self.new(self.x*x - self.y*y, self.x*y + self.y*x)
        else:
            return NotImplemented

    __rmul__ = __mul__

    def __pow__(self, exp):
        if exp == 0:
            return self.new(1, 0)
        if exp < 0:
            self, exp = 1/self, -exp
        if exp == 1:
            return self
        pow2 = self
        prod = self if exp % 2 else self._parent.one
        exp //= 2
        while exp:
            pow2 *= pow2
            if exp % 2:
                prod *= pow2
            exp //= 2
        return prod

    def __bool__(self):
        return bool(self.x) or bool(self.y)

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

    def __rdivmod__(self, other):
        try:
            other = self._parent.convert(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return other.__divmod__(self)

    def __rtruediv__(self, other):
        try:
            other = QQ_I.convert(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return other.__truediv__(self)

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


class GaussianInteger(GaussianElement):
    base = ZZ

    def __truediv__(self, other):
        """Return a Gaussian rational."""
        return QQ_I.convert(self)/other

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


class GaussianRational(GaussianElement):
    base = QQ

    def __truediv__(self, other):
        """Return a Gaussian rational."""
        if not other:
            raise ZeroDivisionError('{} / 0'.format(self))
        x, y = self._get_xy(other)
        if x is None:
            return NotImplemented
        c = x*x + y*y

        return GaussianRational((self.x*x + self.y*y)/c,
                                (-self.x*y + self.y*x)/c)

    def __divmod__(self, other):
        try:
            other = self._parent.convert(other)
        except CoercionFailed:
            return NotImplemented
        if not other:
            raise ZeroDivisionError('{} % 0'.format(self))
        else:
            return self/other, QQ_I.zero


class GaussianDomain():
    """Base class for Gaussian domains."""
    dom = None  # type: Domain

    is_Numerical = True
    is_Exact = True

    has_assoc_Ring = True
    has_assoc_Field = True

    def to_sympy(self, a):
        """Convert ``a`` to a SymPy object. """
        conv = self.dom.to_sympy
        return conv(a.x) + I*conv(a.y)

    def from_sympy(self, a):
        """Convert a SymPy object to ``self.dtype``."""
        r, b = a.as_coeff_Add()
        x = self.dom.from_sympy(r)  # may raise CoercionFailed
        if not b:
            return self.new(x, 0)
        r, b = b.as_coeff_Mul()
        y = self.dom.from_sympy(r)
        if b is I:
            return self.new(x, y)
        else:
            raise CoercionFailed("{} is not Gaussian".format(a))

    def inject(self, *gens):
        """Inject generators into this domain. """
        return self.poly_ring(*gens)

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

    def from_ZZ_gmpy(K1, a, K0):
        """Convert a GMPY mpz to ``self.dtype``."""
        return K1(a)

    def from_ZZ_python(K1, a, K0):
        """Convert a ZZ_python element to ``self.dtype``."""
        return K1(a)

    def from_QQ_gmpy(K1, a, K0):
        """Convert a GMPY mpq to ``self.dtype``."""
        return K1(a)

    def from_QQ_python(K1, a, K0):
        """Convert a QQ_python element to ``self.dtype``."""
        return K1(a)

    def from_AlgebraicField(K1, a, K0):
        """Convert an element from ZZ<I> or QQ<I> to ``self.dtype``."""
        if K0.ext.args[0] == I:
            return K1.from_sympy(K0.to_sympy(a))


class GaussianIntegerRing(GaussianDomain, Ring):
    """Ring of Gaussian integers."""
    dom = ZZ
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
        """Greatest common divisor of a and b over ZZ_I."""
        while b:
            a, b = b, a % b
        return self.normalize(a)

    def lcm(self, a, b):
        """Least common multiple of a and b over ZZ_I."""
        return (a * b) // self.gcd(a, b)

    def from_GaussianIntegerRing(K1, a, K0):
        """Convert a ZZ_I element to ZZ_I."""
        return a

    def from_GaussianRationalField(K1, a, K0):
        """Convert a QQ_I element to ZZ_I."""
        return K1.new(ZZ.convert(a.x), ZZ.convert(a.y))

ZZ_I = GaussianInteger._parent = GaussianIntegerRing()


class GaussianRationalField(GaussianDomain, Field):
    """Field of Gaussian rational numbers."""
    dom = QQ
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

    def as_AlgebraicField(self):
        """Get equivalent domain as an ``AlgebraicField``. """
        return AlgebraicField(self.dom, I)

    def numer(self, a):
        """Get the numerator of ``a``."""
        ZZ_I = self.get_ring()
        return ZZ_I.convert(a * self.denom(a))

    def denom(self, a):
        """Get the denominator of ``a``."""
        ZZ = self.dom.get_ring()
        QQ = self.dom
        ZZ_I = self.get_ring()
        denom_ZZ = ZZ.lcm(QQ.denom(a.x), QQ.denom(a.y))
        return ZZ_I(denom_ZZ, ZZ.zero)

    def from_GaussianIntegerRing(K1, a, K0):
        """Convert a ZZ_I element to QQ_I."""
        return K1.new(a.x, a.y)

    def from_GaussianRationalField(K1, a, K0):
        """Convert a QQ_I element to QQ_I."""
        return a

QQ_I = GaussianRational._parent = GaussianRationalField()

"""Implementation of :class:`RealField` class. """


from sympy.core.numbers import Float
from mpmath.libmp.libmpf import prec_to_dps
from sympy.polys.domains.field import Field
from sympy.polys.domains.simpledomain import SimpleDomain
from sympy.polys.domains.characteristiczero import CharacteristicZero
from sympy.polys.domains.mpelements import MPContext, RealElement
from sympy.polys.polyerrors import CoercionFailed
from sympy.utilities import public

@public
class RealField(Field, CharacteristicZero, SimpleDomain):
    """Real numbers up to the given precision. """

    rep = 'RR'

    is_RealField = is_RR = True

    is_Exact = False
    is_Numerical = True
    is_PID = False

    has_assoc_Ring = False
    has_assoc_Field = True

    _default_precision = 53

    @property
    def has_default_precision(self):
        return self.precision == self._default_precision

    @property
    def precision(self):
        return self._context.prec

    @property
    def dps(self):
        return self._context.dps

    @property
    def tolerance(self):
        return self._context.tolerance

    def __init__(self, prec=_default_precision, dps=None, tol=None):
        context = MPContext(prec, dps, tol, True)
        context._parent = self
        self._context = context

        self.dtype = context.mpf
        self.zero = self.dtype(0)
        self.one = self.dtype(1)

    def __eq__(self, other):
        return (isinstance(other, RealField)
           and self.precision == other.precision
           and self.tolerance == other.tolerance)

    def __hash__(self):
        return hash((self.__class__.__name__, self.dtype, self.precision, self.tolerance))

    def to_sympy(self, element):
        """Convert ``element`` to SymPy number. """
        _mpf = element._mpf_
        dps = prec_to_dps(_mpf[-1])
        return Float(_mpf, dps)

    def from_sympy(self, expr):
        """Convert SymPy's number to ``dtype``. """
        if isinstance(expr, Float):
            # I don't know how else to get the converted
            # expr to look like the original. I tried
            # passing (expr._mpf_, prec=expr._prec) but
            # that did not pass the added tests
            return RealElement(str(expr))
        number = expr.evalf(n=self.dps)
        if number.is_Number:
            return self.dtype(number)
        else:
            raise CoercionFailed("expected real number, got %s" % expr)

    def from_ZZ(self, element, base):
        return self.dtype(element)

    def from_ZZ_python(self, element, base):
        return self.dtype(element)

    def from_QQ(self, element, base):
        return self.dtype(element.numerator) / element.denominator

    def from_QQ_python(self, element, base):
        return self.dtype(element.numerator) / element.denominator

    def from_ZZ_gmpy(self, element, base):
        return self.dtype(int(element))

    def from_QQ_gmpy(self, element, base):
        return self.dtype(int(element.numerator)) / int(element.denominator)

    def from_AlgebraicField(self, element, base):
        return self.from_sympy(base.to_sympy(element).evalf(self.dps))

    def from_RealField(self, element, base):
        if self == base:
            return element
        else:
            return self.dtype(element)

    def from_ComplexField(self, element, base):
        if not element.imag:
            return self.dtype(element.real)

    def to_rational(self, element, limit=True):
        """Convert a real number to rational number. """
        return self._context.to_rational(element, limit)

    def get_ring(self):
        """Returns a ring associated with ``self``. """
        return self

    def get_exact(self):
        """Returns an exact domain associated with ``self``. """
        from sympy.polys.domains import QQ
        return QQ

    def gcd(self, a, b):
        """Returns GCD of ``a`` and ``b``. """
        return self.one

    def lcm(self, a, b):
        """Returns LCM of ``a`` and ``b``. """
        return a*b

    def almosteq(self, a, b, tolerance=None):
        """Check if ``a`` and ``b`` are almost equal. """
        return self._context.almosteq(a, b, tolerance)


RR = RealField()

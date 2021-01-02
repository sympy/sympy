"""Rational number type based on Python integers. """


import operator

from sympy.core.numbers import Rational, Integer
from sympy.core.sympify import converter
from sympy.polys.polyutils import PicklableWithSlots
from sympy.polys.domains.domainelement import DomainElement
from sympy.printing.defaults import DefaultPrinting
from sympy.utilities import public

@public
class PythonRational(DefaultPrinting, PicklableWithSlots, DomainElement):
    """
    Rational number type based on Python integers.

    This is the domain element ``dtype`` for the :py:class:`~.Domain`
    :ref:`QQ` representing rational numbers if ``gmpy`` is not installed.

    This was originally needed for compatibility with older Python versions
    which don't support Fraction. It should probably be removed now and
    replaced by the stdlib Fraction implementation. If ``gmpy`` is installed
    then ``gmpy.mpq`` is used which is much faster than
    :py:class:`PythonRational`. Otherwise this implementation is slower than
    ``Fraction`` so it should be just be removed.

    Examples
    ========

    >>> from sympy.polys.domains import PythonRational

    >>> PythonRational(1)
    1
    >>> PythonRational(2, 3)
    2/3
    >>> PythonRational(14, 10)
    7/5

    """

    __slots__ = ('p', 'q')

    def parent(self):
        from sympy.polys.domains import PythonRationalField
        return PythonRationalField()

    def __init__(self, p, q=1, _gcd=True):
        from sympy.polys.domains.groundtypes import python_gcd as gcd
        if isinstance(p, Integer):
            p = p.p
        elif isinstance(p, Rational):
            p, q = p.p, p.q

        if not q:
            raise ZeroDivisionError('rational number')
        elif q < 0:
            p, q = -p, -q

        if not p:
            self.p = 0
            self.q = 1
        elif p == 1 or q == 1:
            self.p = p
            self.q = q
        else:
            if _gcd:
                x = gcd(p, q)

                if x != 1:
                    p //= x
                    q //= x

            self.p = p
            self.q = q

    @classmethod
    def new(cls, p, q):
        obj = object.__new__(cls)
        obj.p = p
        obj.q = q
        return obj

    def __hash__(self):
        if self.q == 1:
            return hash(self.p)
        else:
            return hash((self.p, self.q))

    def __int__(self):
        p, q = self.p, self.q
        if p < 0:
            return -(-p//q)
        return p//q

    def __float__(self):
        return float(self.p)/self.q

    def __abs__(self):
        return self.new(abs(self.p), self.q)

    def __pos__(self):
        return self.new(+self.p, self.q)

    def __neg__(self):
        return self.new(-self.p, self.q)

    def __add__(self, other):
        from sympy.polys.domains.groundtypes import python_gcd as gcd
        if isinstance(other, PythonRational):
            ap, aq, bp, bq = self.p, self.q, other.p, other.q
            g = gcd(aq, bq)
            if g == 1:
                p = ap*bq + aq*bp
                q = bq*aq
            else:
                q1, q2 = aq//g, bq//g
                p, q = ap*q2 + bp*q1, q1*q2
                g2 = gcd(p, g)
                p, q = (p // g2), q * (g // g2)
        elif isinstance(other, int):
            p = self.p + self.q*other
            q = self.q
        else:
            return NotImplemented

        return self.__class__(p, q, _gcd=False)

    def __radd__(self, other):
        if not isinstance(other, int):
            return NotImplemented

        p = self.p + self.q*other
        q = self.q

        return self.__class__(p, q, _gcd=False)

    def __sub__(self, other):
        from sympy.polys.domains.groundtypes import python_gcd as gcd
        if isinstance(other, PythonRational):
            ap, aq, bp, bq = self.p, self.q, other.p, other.q
            g = gcd(aq, bq)
            if g == 1:
                p = ap*bq - aq*bp
                q = bq*aq
            else:
                q1, q2 = aq//g, bq//g
                p, q = ap*q2 - bp*q1, q1*q2
                g2 = gcd(p, g)
                p, q = (p // g2), q * (g // g2)
        elif isinstance(other, int):
            p = self.p - self.q*other
            q = self.q
        else:
            return NotImplemented

        return self.__class__(p, q, _gcd=False)

    def __rsub__(self, other):
        if not isinstance(other, int):
            return NotImplemented

        p = self.q*other - self.p
        q = self.q

        return self.__class__(p, q, _gcd=False)

    def __mul__(self, other):
        from sympy.polys.domains.groundtypes import python_gcd as gcd
        if isinstance(other, PythonRational):
            ap, aq, bp, bq = self.p, self.q, other.p, other.q
            x1 = gcd(ap, bq)
            x2 = gcd(bp, aq)
            p, q = ((ap//x1)*(bp//x2), (aq//x2)*(bq//x1))
        elif isinstance(other, int):
            x = gcd(other, self.q)
            p = self.p*(other//x)
            q = self.q//x
        else:
            return NotImplemented

        return self.__class__(p, q, _gcd=False)

    def __rmul__(self, other):
        from sympy.polys.domains.groundtypes import python_gcd as gcd
        if not isinstance(other, int):
            return NotImplemented

        x = gcd(self.q, other)
        p = self.p*(other//x)
        q = self.q//x

        return self.__class__(p, q, _gcd=False)

    def __truediv__(self, other):
        from sympy.polys.domains.groundtypes import python_gcd as gcd
        if isinstance(other, PythonRational):
            ap, aq, bp, bq = self.p, self.q, other.p, other.q
            x1 = gcd(ap, bp)
            x2 = gcd(bq, aq)
            p, q = ((ap//x1)*(bq//x2), (aq//x2)*(bp//x1))
        elif isinstance(other, int):
            x = gcd(other, self.p)
            p = self.p//x
            q = self.q*(other//x)
        else:
            return NotImplemented

        return self.__class__(p, q, _gcd=False)

    def __rtruediv__(self, other):
        from sympy.polys.domains.groundtypes import python_gcd as gcd
        if not isinstance(other, int):
            return NotImplemented

        x = gcd(self.p, other)
        p = self.q*(other//x)
        q = self.p//x

        return self.__class__(p, q)

    def __mod__(self, other):
        return self.__class__(0)

    def __divmod__(self, other):
        return (self//other, self % other)

    def __pow__(self, exp):
        p, q = self.p, self.q

        if exp < 0:
            p, q, exp = q, p, -exp

        return self.__class__(p**exp, q**exp, _gcd=False)

    def __bool__(self):
        return self.p != 0

    def __eq__(self, other):
        if isinstance(other, PythonRational):
            return self.q == other.q and self.p == other.p
        elif isinstance(other, int):
            return self.q == 1 and self.p == other
        else:
            return NotImplemented

    def __ne__(self, other):
        return not self == other

    def _cmp(self, other, op):
        try:
            diff = self - other
        except TypeError:
            return NotImplemented
        else:
            return op(diff.p, 0)

    def __lt__(self, other):
        return self._cmp(other, operator.lt)

    def __le__(self, other):
        return self._cmp(other, operator.le)

    def __gt__(self, other):
        return self._cmp(other, operator.gt)

    def __ge__(self, other):
        return self._cmp(other, operator.ge)

    @property
    def numer(self):
        return self.p

    @property
    def denom(self):
        return self.q

    numerator = numer
    denominator = denom


def sympify_pythonrational(arg):
    return Rational(arg.p, arg.q)
converter[PythonRational] = sympify_pythonrational

"""Rational number type based on Python integers. """

from sympy.core.numbers import igcd

import operator

class PythonRationalType(object):
    """
    Rational number type based on Python integers.

    **Examples**

    >>> from sympy.polys.domains import PythonRationalType

    >>> PythonRationalType(1)
    1/1
    >>> PythonRationalType(2, 3)
    2/3
    >>> PythonRationalType(14, 10)
    7/5

    """

    __slots__ = ['p', 'q']

    def __init__(self, p, q=None):
        if q is None:
            self.p = p
            self.q = 1
        else:
            if not q:
                raise ZeroDivisionError('rational number')
            elif q < 0:
                p, q = -p, -q

            g = igcd(p, q)

            self.p = p//g
            self.q = q//g

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

    def __repr__(self):
        return "%s(%d, %d)" % (self.__class__.__name__, self.p, self.q)

    def __str__(self):
        return "%d/%d" % (self.p, self.q)

    def __int__(self):
        return int(float(self.p)/self.q)

    def __float__(self):
        return float(self.p)/self.q

    def __abs__(self):
        return self.new(abs(self.p), self.q)

    def __pos__(self):
        return self.new(+self.p, self.q)

    def __neg__(self):
        return self.new(-self.p, self.q)

    def __add__(self, other):
        if isinstance(other, PythonRationalType):
            p = self.p*other.q + self.q*other.p
            q = self.q*other.q
        elif isinstance(other, (int, long)):
            p = self.p + self.q*other
            q = self.q
        else:
            return NotImplemented

        return self.__class__(p, q)

    def __radd__(self, other):
        if not isinstance(other, (int, long)):
            return NotImplemented

        p = self.p + self.q*other
        q = self.q

        return self.__class__(p, q)

    def __sub__(self, other):
        if isinstance(other, PythonRationalType):
            p = self.p*other.q - self.q*other.p
            q = self.q*other.q
        elif isinstance(other, (int, long)):
            p = self.p - self.q*other
            q = self.q
        else:
            return NotImplemented

        return self.__class__(p, q)

    def __rsub__(self, other):
        if not isinstance(other, (int, long)):
            return NotImplemented

        p = self.q*other - self.p
        q = self.q

        return self.__class__(p, q)

    def __mul__(self, other):
        if isinstance(other, PythonRationalType):
            p = self.p*other.p
            q = self.q*other.q
        elif isinstance(other, (int, long)):
            p = self.p*other
            q = self.q
        else:
            return NotImplemented

        return self.__class__(p, q)

    def __rmul__(self, other):
        if not isinstance(other, (int, long)):
            return NotImplemented

        p = self.p*other
        q = self.q

        return self.__class__(p, q)

    def __div__(self, other):
        if isinstance(other, PythonRationalType):
            p = self.p*other.q
            q = self.q*other.p
        elif isinstance(other, (int, long)):
            p = self.p
            q = self.q*other
        else:
            return NotImplemented

        return self.__class__(p, q)

    __truediv__ = __div__

    def __rdiv__(self, other):
        if not isinstance(other, (int, long)):
            return NotImplemented

        p = self.q*other
        q = self.p

        return self.__class__(p, q)

    __rtruediv__ = __rdiv__

    def __mod__(self, other):
        return self.__class__(0)

    def __divmod__(self, other):
        return (self//other, self%other)

    def __pow__(self, exp):
        p, q = self.p, self.q

        if exp < 0:
            p, q, exp = q, p, -exp

        return self.new(p**exp, q**exp)

    def __nonzero__(self):
        return self.p != 0

    def __eq__(self, other):
        if isinstance(other, PythonRationalType):
            return self.q == other.q and self.p == other.p
        elif isinstance(other, (int, long)):
            return self.q == 1 and self.p == other
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

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

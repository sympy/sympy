"""
Interval arithmetic with correct rounding.

"""

from numerics_float import *
from basic import sympify


class _Rounding:
    def begin(self): self.stored = Float.getmode()
    def end(self): Float.setmode(self.stored)
    def down(self): Float.setmode(ROUND_FLOOR)
    def up(self): Float.setmode(ROUND_CEILING)

rounding = _Rounding()


class Interval(Real):
    """
    An Interval represents the set of all real numbers between two
    endpoints a and b. If the interval is closed (which is what is
    implemented here), the endpoints are themselves considered part
    of the interval. That is, the interval between a and b, denoted
    by [a, b], is the set of points x satisfying a <= x <= b.
    """

    def __new__(cls, a, b=None):
        """
        Interval(a) creates an exact interval (width 0)
        Interval(a, b) creates the interval [a, b]
        """
        a = sympify(a)
        if b is None:
            if isinstance(a, Interval):
                return a
            else:
                return Interval(a, a)
        else:
            b = sympify(b)
            assert a <= b, "endpoints must be properly ordered"
            self = object.__new__(cls)
            self.a, self.b = a, b
            return self

    def __repr__(self):
        return "Interval(%r, %r)" % (self.a, self.b)

    def __str__(self):
        return '[%s, %s]' % (self.a, self.b)

    def __hash__(self):
        return hash((self.a, self.b))

    def __eq__(self, other):
        """Two intervals are considered equal if all endpoints are equal"""
        other = Interval(other)
        return (self.a, self.b) == (other.a, other.b)

    def __contains__(self, x):
        """Return True if x is contained in the interval, otherwise False."""
        return (self.a <= x) and (x <= self.b)

    def __neg__(self):
        return Interval(-self.b, -self.a)

    def __add__(l, r):
        r = Interval(r)
        rounding.begin()
        rounding.down()
        a = l.a + r.a
        rounding.up()
        b = l.b + r.b
        rounding.end()
        return Interval(a, b)

    __radd__ = __add__

    def __sub__(l, r):
        return l + (-r)

    def __rsub__(r, l):
        return -(r - l)

    def __mul__(l, r):
        r = Interval(r)
        rounding.begin()
        rounding.down()
        xd, yd, zd, wd = l.a*r.a, l.a*r.b, l.b*r.a, l.b*r.b
        rounding.up()
        xu, yu, zu, wu = l.a*r.a, l.a*r.b, l.b*r.a, l.b*r.b
        rounding.end()
        return Interval(min(xd,yd,zd,wd), max(xu,yu,zu,wu))

    __rmul__ = __mul__

    def __div__(l, r):
        r = Interval(r)
        if 0 in r:
            raise ZeroDivisionError, "cannot divide by interval containing 0"
        rounding.begin()
        rounding.down()
        xd, yd, zd, wd = l.a/r.a, l.a/r.b, l.b/r.a, l.b/r.b
        rounding.up()
        xu, yu, zu, wu = l.a/r.a, l.a/r.b, l.b/r.a, l.b/r.b
        rounding.end()
        return Interval(min(xd,yd,zd,wd), max(xu,yu,zu,wu))

    def __rdiv__(r, l):
        return Interval(l) / r

    __truediv__ = __div__
    __rtruediv__ = __rdiv__
    __floordiv__ = __div__
    __rfloordiv__ = __rdiv__

    def mid(self):
        return (self.a+self.b)/2

    @staticmethod
    def from_absolute_error(x, error=0):
        return Interval(x-error, x+error)

    @staticmethod
    def from_relative_error(x, error=0):
        # XXX: round
        return Interval(x*(1-error), x*(1+error))

    def absolute_error(self):
        rounding.begin()
        rounding.up()
        p = abs((self.b-self.a)/2)
        rounding.end()
        return p

    def relative_error(self):
        return self.absolute_error() / self.mid()


from __future__ import print_function, division

from sympy.core import Add, Mul
from sympy.core.sympify import sympify
from sympy.core.sets import Interval, FiniteSet
from sympy.core.expr import Expr
from sympy.core.numbers import Number, oo, _sympifyit
from sympy.core.singleton import S


class IV(Interval, Expr):
    """
    Class Object to used as Interval for interval arithmetics
    It differs from the ivf of mpmaths in sense that it can handle symbolic
    objects like `pi`, `E` and `log(2)`. Do not use it for floating point
    interval arithmetic calculations use mpmath.iv instead.
    """

    def _intersect(self, other):
        res = Interval(self.start, self.end).intersect(Interval(
            other.start, other.end))
        return IV(res.start, res.end)

    @property
    def delta(self):
        return self.end - self.start

    @property
    def mid(self):
        return (self.start + self.end)/2

    @_sympifyit('other', NotImplemented)
    def _eval_power(self, other):
        return self.__pow__(other)

    @_sympifyit('other', NotImplemented)
    def __add__(self, other):
        if isinstance(other, IV):
            return IV(Add(self.start, other.start), Add(self.end, other.end))
        elif Number.is_real:
            return IV(Add(self.start, other), Add(self.end, other))
        return NotImplemented
    __radd__ = __add__

    def __neg__(self):
        return IV(-self.end, -self.start)

    @_sympifyit('other', NotImplemented)
    def __sub__(self, other):
        if isinstance(other, IV):
            return IV(Add(self.start, -other.end), Add(self.end, -other.start))
        elif Number.is_real:
            return IV(Add(self.start, -other), Add(self.end, - other))
        return NotImplemented

    @_sympifyit('other', NotImplemented)
    def __rsub__(self, other):
        return - self.__sub__(other)

    @_sympifyit('other', NotImplemented)
    def __mul__(self, other):
        other = sympify(other)
        from sympy.functions.elementary.miscellaneous import Min, Max
        if isinstance(other, IV):
            return IV(Min(Mul(self.start, other.start),
                          Mul(self.start, other.end),
                          Mul(self.end, other.start),
                          Mul(self.end, other.end)),
                      Max(Mul(self.start, other.start),
                          Mul(self.start, other.end),
                          Mul(self.end, other.start),
                          Mul(self.end, other.end)))

        elif Number.is_real:
            if other > S.Zero:
                return IV(Mul(self.start, other), Mul(self.end, other))
            elif other < S.Zero:
                return IV(Mul(self.end, other), Mul(self.start, other))
            else:
                if self.start > -oo and self.end < oo:
                    return FiniteSet(0)
                else:
                    return S.EmptySet()
        return NotImplemented

    __rmul__ = __mul__

    @_sympifyit('other', NotImplemented)
    def __div__(self, other):
        from sympy.functions.elementary.miscellaneous import Min, Max
        if isinstance(other, IV):
            if other.subset(FiniteSet(S.Zero)):
                return IV(-oo, oo)
            else:
                return IV(Min(Mul(self.start, 1/other.start),
                              Mul(self.start, 1/other.end),
                              Mul(self.end, 1/other.start),
                              Mul(self.end, 1/other.end)),
                          Max(Mul(self.start, 1/other.start),
                              Mul(self.start, 1/other.end),
                              Mul(self.end, 1/other.start),
                              Mul(self.end, 1/other.end)))

        elif other.is_real:
            if other > S.Zero:
                return IV(self.start/other, self.end/other)
            elif other < S.Zero:
                return IV(self.end/other, self.start/other)
            else:
                return NotImplemented
        else:
            return NotImplemented

    __truediv__ = __div__

    @_sympifyit('other', NotImplemented)
    def __rdiv__(self, other):
        from sympy.functions.elementary.miscellaneous import Min, Max
        if other.is_real:
            if self.subset(FiniteSet(S.Zero)):
                return IV(-oo, oo)
            else:
                return IV(Min(other/self.start, other/self.end),
                          Max(other/self.start, other/self.end))
        else:
            return NotImplemented

    __rtruediv__ = __rdiv__

    @_sympifyit('other', NotImplemented)
    def __pow__(self, other):
        from sympy.functions.elementary.miscellaneous import Max
        if Number.is_real:
            if other.is_Integer and other > 0:
                if other % 2 == 0:
                    if self.start >= 0:
                        return IV(self.start**other, self.end**other)
                    elif self.end < 0:
                        return IV(self.end**other, self.start**other)
                    else:
                        return IV(0, Max(self.start**other, self.end**other))
                else:
                    return IV(self.start**other, self.end**other)
        return NotImplemented

    def __abs__(self):
        if self.end < S.Zero:
            return self.__neg__()
        elif self.start < S.Zero:
            return IV(S.Zero, self.end)
        else:
            return self

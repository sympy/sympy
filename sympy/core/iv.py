from __future__ import print_function, division

from sympy.core import Add, Mul
from sympy.core.sympify import sympify
from sympy.core.sets import Interval, FiniteSet
from sympy.core.expr import AtomicExpr
from sympy.core.numbers import oo, _sympifyit
from sympy.core.singleton import S


class IV(Interval, AtomicExpr):
    """
    IV represents an interval `[a, b]` that is an arbitary real number greater
    than equal to `a` and less than equal to `b`
    `[a,b] = \{x \in \mathbb{R} \,|\, a \le x \le b\}`
    where `a` and `b` are real numbers.

    Operation on intervals are defined as

    `[a, b] + [c, d] = \{ x+y| a\le x \le b, c \le y \le d}`

    `[a, b] - [c, d] = \{ x-y| a\le x \le b, c \le y \le d}`

    `[a, b] * [c, d] = \{ x*y | a \le x \le b, c \le y \le d\}`

    `[a, b] / [c, d] = \{ x/y | a \le x \le b, c \le y \le d\}`

    `[a, b]^n = \{ x^n | a \le x \le b}`

    Examples
    ========

    >>> from sympy import IV, sin, exp, log, pi, E
    >>> from sympy.abc import x

    >>> IV(0, 1) + IV(1, 2)
    [0, 3]

    >>> IV(0, 1) - IV(0, 2)
    [-1, 1]

    >>> IV(-2, 3)*IV(-1, 1)
    [-3, 2]

    >>> IV(1, 2)*IV(3, 5)
    [1/5, 2/3]

    Note: `[a, b]^2` is not same as `[a, b]*[a, b]`

    >>> IV(-1, 1)**2
    [0, 1]

    Some elementary functions can also take intervals as input.
    A function `f` evaluated for some interval `[a, b]` is defined as
    `f([a, b]) = \{ f(x) | a \le x \le b \}`

    >>> sin(IV(pi/6, pi/3))
    [1/2, sqrt(3)/2]

    >>> exp(IV(0, 1))
    [1, E]

    >>> log(IV(1, E))
    [0, 1]

    Some symbol in an experssion can be substituted for an Interval.
    But it doesn't necessarily evaluate the Interval for that expression.

    >>> (x**2 + 2*x + 1).subs(x, IV(-1, 1))
    [-1, 4]

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Interval_arithmetic

    Notes
    =====

    Do not use ``IV`` for floating point interval arithmetic calculations
    use ``mpmath.iv`` instead.
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
        elif other.is_real:
            return IV(Add(self.start, other), Add(self.end, other))
        return NotImplemented
    __radd__ = __add__

    def __neg__(self):
        return IV(-self.end, -self.start)

    @_sympifyit('other', NotImplemented)
    def __sub__(self, other):
        if isinstance(other, IV):
            return IV(Add(self.start, -other.end), Add(self.end, -other.start))
        elif other.is_real:
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

        elif other.is_real:
            if other.is_positive:
                return IV(Mul(self.start, other), Mul(self.end, other))
            elif other.is_negative:
                return IV(Mul(self.end, other), Mul(self.start, other))
            elif other.is_zero:
                if self.start is not -oo and self.end is not oo:
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
            if other.is_positive:
                return IV(self.start/other, self.end/other)
            elif other.is_negative:
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
        if other.is_real:
            if other.is_Integer and other.is_positive:
                if other % 2 == 0:
                    if self.start.is_positive:
                        return IV(self.start**other, self.end**other)
                    elif self.end.is_negative:
                        return IV(self.end**other, self.start**other)
                    else:
                        return IV(0, Max(self.start**other, self.end**other))
                else:
                    return IV(self.start**other, self.end**other)
        return NotImplemented

    def __abs__(self):
        if self.end.is_negative:
            return self.__neg__()
        elif self.start.is_negative:
            return IV(S.Zero, self.end)
        else:
            return self

from sympy.core import Basic, Expr
from sympy.core.numbers import oo
from sympy.core.symbol import symbols
from sympy.multipledispatch import Dispatcher
from sympy.sets.setexpr import set_mul
from sympy.sets.sets import FiniteSet, Union, Interval, Set
from sympy.sets.fancysets import ImageSet
from typing import Any


_x, _y = symbols("x y")


_set_mul: Dispatcher = Dispatcher('_set_mul')
_set_div: Dispatcher = Dispatcher('_set_div')


@_set_mul.register(Basic, Basic)
def _(x, y) -> FiniteSet | Interval | Any | ImageSet | Union:
    return None

@_set_mul.register(Set, Set)
def _(x, y) -> FiniteSet | Interval | Any | ImageSet | Union:
    return None

@_set_mul.register(Expr, Expr)
def _(x, y) -> FiniteSet | Interval | Any | ImageSet | Union:
    return x*y

@_set_mul.register(Interval, Interval)
def _(x, y) -> FiniteSet | Interval | Any | ImageSet | Union:
    """
    Multiplications in interval arithmetic
    https://en.wikipedia.org/wiki/Interval_arithmetic
    """
    # TODO: some intervals containing 0 and oo will fail as 0*oo returns nan.
    comvals = (
        (x.start * y.start, bool(x.left_open or y.left_open)),
        (x.start * y.end, bool(x.left_open or y.right_open)),
        (x.end * y.start, bool(x.right_open or y.left_open)),
        (x.end * y.end, bool(x.right_open or y.right_open)),
    )
    # TODO: handle symbolic intervals
    minval, minopen = min(comvals)
    maxval, maxopen = max(comvals)
    return Interval(
        minval,
        maxval,
        minopen,
        maxopen
    )

@_set_div.register(Basic, Basic)
def _(x, y) -> FiniteSet | Interval | Any | ImageSet | Union:
    return None

@_set_div.register(Expr, Expr)
def _(x, y) -> FiniteSet | Interval | Any | ImageSet | Union:
    return x/y

@_set_div.register(Set, Set)
def _(x, y) -> FiniteSet | Interval | Any | ImageSet | Union:
    return None

@_set_div.register(Interval, Interval)
def _(x, y) -> FiniteSet | Interval | Any | ImageSet | Union:
    """
    Divisions in interval arithmetic
    https://en.wikipedia.org/wiki/Interval_arithmetic
    """
    if (y.start*y.end).is_negative:
        return Interval(-oo, oo)
    if y.start == 0:
        s2 = oo
    else:
        s2 = 1/y.start
    if y.end == 0:
        s1 = -oo
    else:
        s1 = 1/y.end
    return set_mul(x, Interval(s1, s2, y.right_open, y.left_open))

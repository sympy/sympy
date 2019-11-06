from sympy import symbols, S
from sympy.core import Basic, Expr
from sympy.core.numbers import Infinity, NegativeInfinity
from sympy.multipledispatch import dispatch
from sympy.sets import Interval, FiniteSet


_x, _y = symbols("x y")


@dispatch(Basic, Basic)
def _set_add(x, y):
    return None

@dispatch(Expr, Expr)
def _set_add(x, y):
    return x+y

@dispatch(Interval, Interval)
def _set_add(x, y):
    """
    Additions in interval arithmetic
    https://en.wikipedia.org/wiki/Interval_arithmetic
    """
    return Interval(x.start + y.start, x.end + y.end,
        x.left_open or y.left_open, x.right_open or y.right_open)

@dispatch(Interval, Infinity)
def _set_add(x, y):
    if x.start is S.NegativeInfinity:
        return Interval(-oo, oo)
    return FiniteSet({S.Infinity})

@dispatch(Interval, NegativeInfinity)
def _set_add(x, y):
    if x.end is S.Infinity:
        return Interval(-oo, oo)
    return FiniteSet({S.NegativeInfinity})


@dispatch(Basic, Basic)
def _set_sub(x, y):
    return None

@dispatch(Expr, Expr)
def _set_sub(x, y):
    return x-y

@dispatch(Interval, Interval)
def _set_sub(x, y):
    """
    Subtractions in interval arithmetic
    https://en.wikipedia.org/wiki/Interval_arithmetic
    """
    return Interval(x.start - y.end, x.end - y.start,
        x.left_open or y.right_open, x.right_open or y.left_open)

@dispatch(Interval, Infinity)
def _set_sub(x, y):
    if self.start is S.NegativeInfinity:
        return Interval(-oo, oo)
    return FiniteSet(-oo)

@dispatch(Interval, NegativeInfinity)
def _set_sub(x, y):
    if self.start is S.NegativeInfinity:
        return Interval(-oo, oo)
    return FiniteSet(-oo)

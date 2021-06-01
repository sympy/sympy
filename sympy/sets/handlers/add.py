from sympy import symbols, S, oo

from sympy.core import Basic, Expr
from sympy.core.numbers import Infinity, NegativeInfinity
from sympy.multipledispatch import dispatch
from sympy.sets import Interval, FiniteSet



# XXX: The functions in this module are clearly not tested and are broken in a
# number of ways.

_x, _y = symbols("x y")


@dispatch(Basic, Basic)  # type: ignore # noqa:F811
def _set_add(x, y): # noqa:F811
    return None


@dispatch(Expr, Expr)  # type: ignore # noqa:F811
def _set_add(x, y): # noqa:F811
    return x+y


@dispatch(Interval, Interval)  # type: ignore # noqa:F811
def _set_add(x, y): # noqa:F811
    """
    Additions in interval arithmetic
    https://en.wikipedia.org/wiki/Interval_arithmetic
    """
    return Interval(x.start + y.start, x.end + y.end,
                    x.left_open or y.left_open, x.right_open or y.right_open)


@dispatch(Interval, Infinity)  # type: ignore # noqa:F811
def _set_add(x, y): # noqa:F811
    if x.start is S.NegativeInfinity:
        return Interval(-oo, oo)
    return FiniteSet({S.Infinity})

@dispatch(Interval, NegativeInfinity)  # type: ignore # noqa:F811
def _set_add(x, y): # noqa:F811
    if x.end is S.Infinity:
        return Interval(-oo, oo)
    return FiniteSet({S.NegativeInfinity})


@dispatch(Basic, Basic)  # type: ignore
def _set_sub(x, y): # noqa:F811
    return None


@dispatch(Expr, Expr)  # type: ignore # noqa:F811
def _set_sub(x, y): # noqa:F811
    return x-y


@dispatch(Interval, Interval)  # type: ignore # noqa:F811
def _set_sub(x, y): # noqa:F811
    """
    Subtractions in interval arithmetic
    https://en.wikipedia.org/wiki/Interval_arithmetic
    """
    return Interval(x.start - y.end, x.end - y.start,
                    x.left_open or y.right_open, x.right_open or y.left_open)


@dispatch(Interval, Infinity)  # type: ignore # noqa:F811
def _set_sub(x, y): # noqa:F811
    if x.start is S.NegativeInfinity:
        return Interval(-oo, oo)
    return FiniteSet(-oo)

@dispatch(Interval, NegativeInfinity)  # type: ignore # noqa:F811
def _set_sub(x, y): # noqa:F811
    if x.start is S.NegativeInfinity:
        return Interval(-oo, oo)
    return FiniteSet(-oo)

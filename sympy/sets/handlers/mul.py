from sympy import Set, symbols
from sympy.core import Basic, Expr
from sympy.multipledispatch import dispatch
from sympy.sets import Interval

_x, _y = symbols("x y")


@dispatch(Basic, Basic)  # type: ignore # noqa:F811
def _set_mul(x, y): # noqa:F811
    return None

@dispatch(Set, Set)  # type: ignore # noqa:F811
def _set_mul(x, y): # noqa:F811
    return None

@dispatch(Expr, Expr)  # type: ignore # noqa:F811
def _set_mul(x, y): # noqa:F811
    return x*y

@dispatch(Interval, Interval)  # type: ignore # noqa:F811
def _set_mul(x, y): # noqa:F811
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

@dispatch(Basic, Basic)  # type: ignore # noqa:F811
def _set_div(x, y): # noqa:F811
    return None

@dispatch(Expr, Expr)  # type: ignore # noqa:F811
def _set_div(x, y): # noqa:F811
    return x/y

@dispatch(Set, Set)  # type: ignore # noqa:F811 # noqa:F811
def _set_div(x, y): # noqa:F811
    return None

@dispatch(Interval, Interval)  # type: ignore # noqa:F811
def _set_div(x, y): # noqa:F811
    """
    Divisions in interval arithmetic
    https://en.wikipedia.org/wiki/Interval_arithmetic
    """
    from sympy.sets.setexpr import set_mul
    from sympy import oo
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

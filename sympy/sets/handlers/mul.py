from sympy.multipledispatch import dispatch, Dispatcher
from sympy.core import Basic, Expr, Function, Add, Mul, Pow, Dummy, Integer
from sympy import Min, Max, Set, sympify, Lambda, symbols, exp, log, S, oo
from sympy.sets import (imageset, Interval, FiniteSet, Union, ImageSet,
    ProductSet, EmptySet, Intersection)
from sympy.core.function import FunctionClass
from sympy.logic.boolalg import And, Or, Not, true, false


_x, _y = symbols("x y")


@dispatch(Basic, Basic)
def _mul_sets(x, y):
    return None

@dispatch(Set, Set)
def _mul_sets(x, y):
    return None

@dispatch(Expr, Expr)
def _mul_sets(x, y):
    return x*y

@dispatch(Interval, Interval)
def _mul_sets(x, y):
    """
    Multiplications in interval arithmetic
    https://en.wikipedia.org/wiki/Interval_arithmetic
    """
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
    return SetExpr(Interval(start, end))


@dispatch(Basic, Basic)
def _div_sets(x, y):
    return None

@dispatch(Expr, Expr)
def _div_sets(x, y):
    return x/y

@dispatch(Set, Set)
def _div_sets(x, y):
    return None

@dispatch(Interval, Interval)
def _div_sets(x, y):
    """
    Divisions in interval arithmetic
    https://en.wikipedia.org/wiki/Interval_arithmetic
    """
    from sympy.sets.setexpr import mul_sets
    from sympy import oo
    if (y.start*y.end).is_negative:
        return Interval(-oo, oo)
    if y.start == 0:
        s2 = oo
    else:
        s2 = 1/y.start
    if y.end == 0:
        s1 = oo
    else:
        s1 = 1/y.end
    return mul_sets(x, Interval(s1, s2, y.right_open, y.left_open))

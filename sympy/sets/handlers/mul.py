from sympy.multipledispatch import dispatch, Dispatcher
from sympy.core import Basic, Expr, Function, Add, Mul, Pow, Dummy, Integer
from sympy import Min, Max, Set, sympify, Lambda, symbols, exp, log, S
from sympy.sets import (imageset, Interval, FiniteSet, Union, ImageSet,
    ProductSet, EmptySet, Intersection)
from sympy.core.function import FunctionClass
from sympy.logic.boolalg import And, Or, Not, true, false


_x, _y = symbols("x y")


@dispatch(Set, Set)
def mul_sets(x, y):
    return ImageSet(Lambda((_x, _y), (_x * _y)), x, y)


@dispatch(Expr, Expr)
def mul_sets(x, y):
    return x*y


@dispatch(Interval, Interval)
def mul_sets(x, y):
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


@dispatch(Expr, Expr)
def div_sets(x, y):
    return x/y


@dispatch(Set, Set)
def div_sets(x, y):
    return ImageSet(Lambda((_x, _y), (_x / _y)), x, y)


@dispatch(Interval, Interval)
def div_sets(x, y):
    """
    Divisions in interval arithmetic
    https://en.wikipedia.org/wiki/Interval_arithmetic
    """
    if (y.start*y.end).is_negative:
        from sympy import oo
        return Interval(-oo, oo)
    return mul_sets(x, Interval(1/y.end, 1/y.start, y.right_open, y.left_open))

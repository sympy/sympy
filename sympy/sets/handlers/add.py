from sympy.multipledispatch import dispatch, Dispatcher
from sympy.core import Basic, Expr, Function, Add, Mul, Pow, Dummy, Integer
from sympy import Min, Max, Set, sympify, Lambda, symbols, exp, log, S
from sympy.sets import (imageset, Interval, FiniteSet, Union, ImageSet,
    ProductSet, EmptySet, Intersection)
from sympy.core.function import FunctionClass
from sympy.logic.boolalg import And, Or, Not, true, false


_x, _y = symbols("x y")


@dispatch(Set, Set)
def add_sets(x, y):
    return ImageSet(Lambda((_x, _y), (_x+_y)), x, y)


@dispatch(Expr, Expr)
def add_sets(x, y):
    return x+y


@dispatch(Interval, Interval)
def add_sets(x, y):
    """
    Additions in interval arithmetic
    https://en.wikipedia.org/wiki/Interval_arithmetic
    """
    return Interval(x.start + y.start, x.end + y.end,
        x.left_open or y.left_open, x.right_open or y.right_open)


@dispatch(Expr, Expr)
def sub_sets(x, y):
    return x-y


@dispatch(Set, Set)
def sub_sets(x, y):
    return ImageSet(Lambda((_x, _y), (_x - _y)), x, y)


@dispatch(Interval, Interval)
def sub_sets(x, y):
    """
    Subtractions in interval arithmetic
    https://en.wikipedia.org/wiki/Interval_arithmetic
    """
    return Interval(x.start - y.end, x.end - y.start,
        x.left_open or y.right_open, x.right_open or y.left_open)

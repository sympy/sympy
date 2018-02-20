from sympy.multipledispatch import dispatch, Dispatcher
from sympy.core import Basic, Expr, Function, Add, Mul, Pow, Dummy, Integer
from sympy import Min, Max, Set, sympify, Lambda, symbols, exp, log, S
from sympy.sets import (imageset, Interval, FiniteSet, Union, ImageSet,
    ProductSet, EmptySet, Intersection)
from sympy.core.function import FunctionClass
from sympy.logic.boolalg import And, Or, Not, true, false


_x, _y = symbols("x y")


@dispatch(Set, Set)
def pow_sets(x, y):
    return ImageSet(Lambda((_x, _y), (_x ** _y)), x, y)


@dispatch(Expr, Expr)
def pow_sets(x, y):
    return x**y


@dispatch(Interval, Integer)
def pow_sets(x, y):
    """
    Powers in interval arithmetic
    https://en.wikipedia.org/wiki/Interval_arithmetic
    """
    exponent = sympify(exponent)
    if exponent.is_odd:
        return Interval(x.start**exponent, x.end**exponent, x.left_open, x.right_open)
    if exponent.is_even:
        if (x.start*x.end).is_negative:
            if -x.start > x.end:
                left_limit = x.start
                left_open = x.right_open
            else:
                left_limit = x.end
                left_open = x.left_open
            return Interval(S.Zero, left_limit ** exponent, S.Zero not in x, left_open)
        elif x.start.is_negative and x.end.is_negative:
            return Interval(x.end**exponent, x.start**exponent, x.right_open, x.left_open)
        else:
            return Interval(x.start**exponent, x.end**exponent, x.left_open, x.right_open)

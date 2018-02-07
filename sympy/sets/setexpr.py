from sympy.core import Basic, Expr, Function, Add, Mul, Pow, Dummy, Integer
from sympy import sift, latex, Min, Max, Set, sympify, Lambda, symbols
from sympy.sets import imageset, Interval, FiniteSet, Union, ImageSet, ProductSet
from sympy.core.decorators import call_highest_priority, _sympifyit
from sympy.multipledispatch import dispatch, Dispatcher

from itertools import count


d = Dummy('d')


class SetExpr(Expr):
    """An expression that can take on values of a set

    >>> from sympy import Interval, FiniteSet
    >>> from sympy.sets.setexpr import SetExpr

    >>> a = SetExpr(Interval(0, 5))
    >>> b = SetExpr(FiniteSet(1, 10))
    >>> (a + b).set
    Union(Interval(1, 6), Interval(10, 15))
    >>> (2*a + b).set
    Interval(1, 20)
    """
    _op_priority = 11.0

    def __new__(cls, setarg):
        return Expr.__new__(cls, setarg)

    set = property(lambda self: self.args[0])

    def _latex(self, printer):
        return printer._print(self.set)

    def _pretty(self, printer):
        return printer._print(self.set)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__radd__')
    def __add__(self, other):
        return apply_operation(add_sets, self, other)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__add__')
    def __radd__(self, other):
        return apply_operation(add_sets, other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return apply_operation(mul_sets, self, other)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return apply_operation(mul_sets, other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return apply_operation(sub_sets, self, other)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return apply_operation(sub_sets, other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        return apply_operation(pow_sets, self, other)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        return apply_operation(pow_sets, other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rdiv__')
    def __div__(self, other):
        return apply_operation(div_sets, self, other)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__div__')
    def __rdiv__(self, other):
        return apply_operation(div_sets, other, self)

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def _eval_func(self, func):
        return SetExpr(imageset(func, self.set))


_x, _y = symbols("x y")


@dispatch(Set, Set)
def add_sets(x, y):
    return ImageSet(Lambda((_x, _y), (_x+_y)), ProductSet(x, y))


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
    return ImageSet(Lambda((_x, _y), (_x - _y)), ProductSet(x, y))


@dispatch(Interval, Interval)
def sub_sets(x, y):
    """
    Subtractions in interval arithmetic
    https://en.wikipedia.org/wiki/Interval_arithmetic
    """
    return Interval(x.start - y.end, x.end - y.start,
        x.left_open or y.right_open, x.right_open or y.left_open)


@dispatch(Set, Set)
def mul_sets(x, y):
    return ImageSet(Lambda((_x, _y), (_x * _y)), ProductSet(x, y))


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
    return ImageSet(Lambda((_x, _y), (_x / _y)), ProductSet(x, y))


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


@dispatch(Set, Set)
def pow_sets(x, y):
    return ImageSet(Lambda((_x, _y), (_x ** _y)), ProductSet(x, y))


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


@dispatch(FiniteSet, FiniteSet, Dispatcher)
def dispatch_on_operation(x, y, op):
    return FiniteSet(*[op(i, j) for i in x for j in y])

@dispatch(Set, FiniteSet, Dispatcher)
def dispatch_on_operation(x, y, op):
    return Union(*[dispatch_on_operation(x, i, op) for i in y])

@dispatch(FiniteSet, Set, Dispatcher)
def dispatch_on_operation(x, y, op):
    return Union(*[dispatch_on_operation(i, y, op) for i in x])

@dispatch(Basic, Basic, Dispatcher)
def dispatch_on_operation(x, y, op):
    return op(x, y)

@dispatch(Expr, Set, Dispatcher)
def dispatch_on_operation(x, y, op):
    return ImageSet(Lambda(d, op(x, d)), y).doit()

@dispatch(Set, Expr, Dispatcher)
def dispatch_on_operation(x, y, op):
    return ImageSet(Lambda(d, op(d, y)), x).doit()


def apply_operation(op, x, y):
    if isinstance(x, SetExpr):
        x = x.set
    if isinstance(y, SetExpr):
        y = y.set
    out = dispatch_on_operation(x, y, op)
    assert isinstance(out, Set)
    return SetExpr(out)

from sympy.core import Basic, Expr, Function, Add, Mul, Pow, Dummy, Integer
from sympy import latex, Min, Max, Set, sympify, Lambda, symbols
from sympy.sets import (imageset, Interval, FiniteSet, Union, ImageSet,
        ProductSet)
from sympy.core.decorators import call_highest_priority, _sympifyit
from sympy.multipledispatch import dispatch, Dispatcher
from sympy.sets.dispatchers import (add_sets, sub_sets, mul_sets, div_sets,
        pow_sets, function_sets)


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
        return r"SetExpr\left({0}\right)".format(printer._print(self.set))

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
        # TODO: this could be implemented straight into `imageset`:
        return SetExpr(function_sets(func, self.set))


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

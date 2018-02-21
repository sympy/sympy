from sympy.core import Basic, Expr, Function, Add, Mul, Pow, Dummy, Integer
from sympy import latex, Min, Max, Set, sympify, Lambda, symbols
from sympy.sets import (imageset, Interval, FiniteSet, Union, ImageSet,
        ProductSet)
from sympy.core.decorators import call_highest_priority, _sympifyit
from sympy.utilities.iterables import sift
from sympy.multipledispatch import dispatch, Dispatcher


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
        return _setexpr_apply_operation(add_sets, self, other)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__add__')
    def __radd__(self, other):
        return _setexpr_apply_operation(add_sets, other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        return _setexpr_apply_operation(mul_sets, self, other)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        return _setexpr_apply_operation(mul_sets, other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rsub__')
    def __sub__(self, other):
        return _setexpr_apply_operation(sub_sets, self, other)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__sub__')
    def __rsub__(self, other):
        return _setexpr_apply_operation(sub_sets, other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        return _setexpr_apply_operation(pow_sets, self, other)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        return _setexpr_apply_operation(pow_sets, other, self)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__rdiv__')
    def __div__(self, other):
        return _setexpr_apply_operation(div_sets, self, other)

    @_sympifyit('other', NotImplemented)
    @call_highest_priority('__div__')
    def __rdiv__(self, other):
        return _setexpr_apply_operation(div_sets, other, self)

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def _eval_func(self, func):
        # TODO: this could be implemented straight into `imageset`:
        res = function_sets(func, self.set)
        if res is None:
            return SetExpr(ImageSet(func, self.set))
        return SetExpr(res)


def _setexpr_apply_operation(op, x, y):
    if isinstance(x, SetExpr):
        x = x.set
    if isinstance(y, SetExpr):
        y = y.set
    out = op(x, y)
    return SetExpr(out)


def _handle_finite_sets(op, x, y, commutative):
    # Handle finite sets:
    fs_args, other = sift([x, y], lambda x: isinstance(x, FiniteSet), binary=True)
    if len(fs_args) == 2:
        return FiniteSet(*[op(i, j) for i in fs_args[0] for j in fs_args[1]])
    elif len(fs_args) == 1:
        sets = [_apply_operation(op, other[0], i, commutative) for i in fs_args[0]]
        return Union(*sets)
    else:
        return None

def _apply_operation(op, x, y, commutative):
    out = _handle_finite_sets(op, x, y, commutative)
    if out is None:
        out = op(x, y)

    if out is None and commutative:
        out = op(y, x)
    if out is None:
        _x, _y = symbols("x y")
        if isinstance(x, Set) and not isinstance(y, Set):
            out = ImageSet(Lambda(d, op(d, y)), x).doit()
        elif not isinstance(x, Set) and isinstance(y, Set):
            out = ImageSet(Lambda(d, op(x, d)), y).doit()
        else:
            out = ImageSet(Lambda((_x, _y), op(_x, _y)), x, y)
    return out

def add_sets(x, y):
    from sympy.sets.handlers.add import _add_sets
    return _apply_operation(_add_sets, x, y, commutative=True)

def sub_sets(x, y):
    from sympy.sets.handlers.add import _sub_sets
    return _apply_operation(_sub_sets, x, y, commutative=False)

def mul_sets(x, y):
    from sympy.sets.handlers.mul import _mul_sets
    return _apply_operation(_mul_sets, x, y, commutative=True)

def div_sets(x, y):
    from sympy.sets.handlers.mul import _div_sets
    return _apply_operation(_div_sets, x, y, commutative=False)

def pow_sets(x, y):
    from sympy.sets.handlers.power import _pow_sets
    return _apply_operation(_pow_sets, x, y, commutative=False)

def function_sets(f, x):
    from sympy.sets.handlers.functions import _function_sets
    return _function_sets(f, x)

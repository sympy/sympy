"""Composite operator."""

from .operator import OperatorExpr, Operator
from sympy.core.add import Add
from sympy.core.compatibility import iterable
from sympy.core.containers import Tuple
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.core.singleton import S

class CompositeOperator(OperatorExpr):
    """
    Composite operator.

    Examples
    ========

    >>> from sympy import sin, cos
    >>> from sympy.abc import x

    >>> sin@cos
    CompositeOperator(sin, cos)
    >>> (sin@cos)(x)
    sin(cos(x))
    >>> (2@sin)(x)
    2
    >>> (cos@3)(x)
    cos(3)

    Parameters
    ==========
    f : instance of Basic or BasicMeta
            Operator which will be evaluated.

    gs : tuple of Basic or BasicMeta
            Operators which will be passed as arguments of f.
    """
    def __new__(cls, f, *gs, **kwargs):
        f = Operator(f)
        gs = [Operator(g) for g in gs]
        gs = Tuple(*gs)
        obj = super(cls, CompositeOperator).__new__(cls, f, *gs, **kwargs)
        obj.f = f
        obj.gs = gs
        return obj

    def _eval_operation(self, *args, **kwargs):
        kwargs.pop('evaluate', None)
        f_args = []
        for g in self.gs:
            f_args.append(g(*args, evaluate=True, **kwargs))
        return self.f(*f_args, evaluate=True, **kwargs)

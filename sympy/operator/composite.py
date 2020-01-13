"""Composite operator."""

from .operator import OpExpr, Op
from sympy.core.containers import Tuple

class CompositeOp(OpExpr):
    """
    Composite operator.

    Examples
    ========

    >>> from sympy import sin, cos
    >>> from sympy.abc import x

    >>> sin@cos
    CompositeOp(sin, cos)
    >>> (sin@cos)(x)
    sin(cos(x))
    >>> (2@sin)(x)
    2
    >>> (cos@3)(x)
    cos(3)

    Parameters
    ==========
    f : instance of Expr or FunctionClass
            Op which will be evaluated.

    gs : tuple of Expr or FunctionClass
            Op which will be passed as arguments of f.
    """
    def __new__(cls, f, *gs, **kwargs):
        f = Op(f)
        gs = [Op(g) for g in gs]
        gs = Tuple(*gs)
        obj = super(cls, CompositeOp).__new__(cls, f, *gs, **kwargs)
        obj.f = f
        obj.gs = gs
        return obj

    def _eval_operation(self, *args, **kwargs):
        kwargs.pop('evaluate', None)
        f_args = []
        for g in self.gs:
            f_args.append(g(*args, evaluate=True, **kwargs))
        return self.f(*f_args, evaluate=True, **kwargs)

"""Operated operators."""
from .operator import OperatorExpr, Operator
from sympy.core.add import Add
from sympy.core.compatibility import iterable
from sympy.core.containers import Tuple
from sympy.core.mul import Mul
from sympy.core.power import Pow
from sympy.core.singleton import S

class CompositeOperator(OperatorExpr):
    def __new__(cls, f, *gs, **kwargs):
        f, = cls._process_basic_basicmeta(f)
        gs = cls._process_basic_basicmeta(*gs)
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

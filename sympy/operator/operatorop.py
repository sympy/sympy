"""Operated operators."""
from .operator import OperatorExpr, Operator
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.power import Pow

class OperatorAdd(OperatorExpr, Add):
    identity = Operator(0)

    def __new__(cls, *args, **kwargs):
        new_args = cls._sep_basic_basicmeta(*args)
        obj = super(cls, OperatorAdd).__new__(cls, *new_args, **kwargs)
        if not isinstance(obj, OperatorExpr):
            return Operator(obj)
        return obj

    def _eval_operation(self, *args, **kwargs):
        args = [a(*args, evaluate=True, **kwargs) if isinstance(a, OperatorExpr) else a for a in self.args]
        return Add(*args)

class OperatorMul(OperatorExpr, Mul):
    identity = Operator(1)

    def __new__(cls, *args, **kwargs):
        new_args = cls._sep_basic_basicmeta(*args)
        obj = super(cls, OperatorMul).__new__(cls, *new_args, **kwargs)
        if not isinstance(obj, OperatorExpr):
            return Operator(obj)
        return obj


    def _eval_operation(self, *args, **kwargs):
        args = [a(*args, evaluate=True, **kwargs) if isinstance(a, OperatorExpr) else a for a in self.args]
        return Mul(*args)

class OperatorPow(OperatorExpr, Pow):
    def __new__(cls, b, e, **kwargs):
        new_b, new_e = cls._sep_basic_basicmeta(b, e)
        obj = super(cls, OperatorPow).__new__(cls, new_b, new_e, **kwargs)
        if not isinstance(obj, OperatorExpr):
            return Operator(obj)
        return obj

    def _eval_operation(self, *args, **kwargs):
        args = [a(*args, evaluate=True, **kwargs) if isinstance(a, OperatorExpr) else a for a in self.args]
        return Pow(*args)

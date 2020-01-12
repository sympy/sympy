"""Contains ``OperatorExpr``, which is a superclass for every operator classes,
and ``Operator`` class, which is a atomic operator class."""
from sympy.core.basic import Basic
from sympy.core.core import BasicMeta
from sympy.core.expr import Expr
from sympy.core.singleton import S
from sympy.core.sympify import sympify

class OperatorExpr(Expr):

    _op_priority = 20 # This is some random big value

    def __call__(self, *args, evaluate=True, **kwargs):
        args = sympify(args)
        return AppliedOperator(self, args, evaluate=evaluate, **kwargs)

    def _eval_operation(self, *args, **kwargs):
        raise NotImplementedError()

    def __pos__(self):
        return self

    def __neg__(self):
        from .operatorop import OperatorMul
        return OperatorMul(-1, self)

    def __add__(self, other):
        from .operatorop import OperatorAdd
        return OperatorAdd(self, other)

    def __radd__(self, other):
        from .operatorop import OperatorAdd
        return OperatorAdd(other, self)

    def __sub__(self, other):
        from .operatorop import OperatorAdd
        return OperatorAdd(self, -other)

    def __rsub__(self, other):
        from .operatorop import OperatorAdd
        return OperatorAdd(-self, other)

    def __mul__(self, other):
        from .operatorop import OperatorMul
        return OperatorMul(self, other)

    def __rmul__(self, other):
        from .operatorop import OperatorMul
        return OperatorMul(other, self)

    def __pow__(self, other):
        from .operatorop import OperatorPow
        return OperatorPow(self, other)

    def __rpow__(self, other):
        from .operatorop import OperatorPow
        return OperatorPow(other, self)

    def __div__(self, other):
        from .operatorop import OperatorMul, OperatorPow
        return OperatorMul(self, OperatorPow(other, S.NegativeOne))

    def __rdiv__(self, other):
        from .operatorop import OperatorMul, OperatorPow
        return OperatorMul(other, OperatorPow(self, S.NegativeOne))

    __truediv__ = __div__
    __rtruediv__ = __rdiv__


    @classmethod
    def _sep_basic_basicmeta(cls, *args):
        args = [sympify(a) for a in args]
        newargs = []
        for a in args:
            if isinstance(a, BasicMeta):
                newargs.append(Operator(a))
            elif isinstance(a, Basic):
                if isinstance(a, Operator):
                    if isinstance(a.operator, BasicMeta):
                        newargs.append(a)
                    elif isinstance(a.operator, Basic):
                        newargs.append(a.operator)
                else:
                    newargs.append(a)
            else:
                newargs.append(a)
        return newargs


class Operator(OperatorExpr):
    def __new__(cls, operator, **kwargs):
        operator = sympify(operator)
        if isinstance(operator, OperatorExpr):
            return operator

        obj = super(cls, Operator).__new__(cls, operator)
        obj.operator = operator
        return obj

    def _eval_operation(self, *args, **kwargs):
        if isinstance(self.operator, Basic):
            return self.operator
        elif isinstance(self.operator, BasicMeta):
            return self.operator(*args, **kwargs)

class AppliedOperator(OperatorExpr):
    def __new__(cls, operator, args, evaluate=False, **kwargs):
        if not isinstance(operator, OperatorExpr):
            raise TypeError("%s must be instance of OperatorExpr." % operator)
        if isinstance(operator, Operator):
            if isinstance(operator.operator, BasicMeta):
                return operator.operator(*args, evaluate=evaluate, **kwargs)
            elif isinstance(operator.operator, Basic):
                return operator.operator

        if evaluate:
            obj = operator._eval_operation(*args, **kwargs)
        else:
            obj = super(cls, AppliedOperator).__new__(cls, operator, args, **kwargs)
            obj.operator = obj.args[0]
            obj.arguments = obj.args[1]

        return obj

    def doit(self, **kwargs):
        return self.func(*self.args, evaluate=True)

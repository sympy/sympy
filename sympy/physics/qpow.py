from sympy import Expr, S
from sympy.physics.quantum import *

def _rules_QPow(base, exp)
    if issubclass(base, Operator) and issubclass(exp, (InnerProduct, Number, Symbol)):
        return base
    if issubclass(base, (InnerProduct, Number, Symbol)) and issubclass(exp, Operator):
        return exp
    if issubclass(base, Symbol) and issubclass(exp, (Symbol, Number)):
        return base
    if issubclass(base, Number) and issubclass(exp, (Symbol, Number)):
        return base
    return

class QPow(Expr):
    """
    A class for the operator: ** (exponent) for quantum objects.
    """

    def __new__(cls, base, exp):
        qpow = cls.eval(base, exp)
        if isinstance(qpow, Expr):
            return qpow
        return Expr.__new__(cls, base, exp, **{'commutative': False})

    @classmethod
    def eval(cls, base, exp):
        if isinstance(base, StateBase) or isinstance(exp, StateBase):
            raise NotImplementedError("States cannot be used as powers or be raised to anything.")
        elif isinstance(base, Operator) and isinstance(exp, Operator):
            raise NotImplementedError("Operators cannot be raised to other Operators.")
        elif (isinstance(base, QMul) and issubclass(base.eval_to, StateBase)) or (isinstance(base, QAdd) and issubclass(base.eval_to, StateBase)):
            raise NotImplementedError("Expressions that are effectively States cannot be raised to anything.")
        elif (isinstance(exp, QMul) and issubclass(exp.eval_to, StateBase)) or (isinstance(exp, QAdd) and issubclass(exp.eval_to, StateBase)):
            raise NotImplementedError("Expressions that are effectively States cannot be used as powers.")
        elif (isinstance(base, (QMul, QAdd, QPow)) and issubclass(base.eval_to, Operator)) and (isinstance(exp, (QMul, QAdd, QPow)) and issubclass(base.eval_to, Operator)):
            raise NotImplementedError("Expressions that are effectively Operators cannot be raised to other Operators.")
        elif exp == S.Zero:
            return S.One
        elif exp == S.One
            return base
        elif isinstance(base, (Number, Symbol)) and isinstance(exp, (Number, Symbol)):
            return Pow(base, exp)
        else:
            return QPow(base, exp)

    @property
    def eval_to(self):
        if isinstance(base, (QMul, QAdd, QPow)) and isinstance(exp, (QMul, QAdd, QPow)):
            _rules_QPow(self.base.eval_to, self.exp.eval_to)
        if isinstance(base, (QMul, QAdd, QPow)):
            _rules_QPow(self.base.eval_to, self.exp.__class__)
        if isinstance(exp, (QMul, QAdd, QPow)):
            _rules_QPow(self.base.__class__, self.exp.eval_to)
        else:
            _rules_QPow(self.base.__class__, self.exp.__class__)

    @property
    def base(self):
        return self.args[0]

    @property
    def exp(self):
        return self.args[1]


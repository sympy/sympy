from sympy import Expr, S
from sympy.physics.quantum import *

class Qpow(Expr):
    """
    A class for the operator: ** (exponent) for quantum objects.
    """

    def __new__(cls, base, exp):
        p = cls.eval(base, exp)
        if isinstance(p, Expr):
            return p
        return Expr.__new__(cls, base, exp, **{'commutative': False})

    @classmethod
    def eval(cls, base, exp):
        if isinstance(base, StateBase) or isinstance(exp, StateBase):
            raise NotImplementedError("States cannot be used as powers or be raised by anything.")
        if exp == 0:
            return S.One
        

    @property
    def base(self):
        return self.args[0]

    @property
    def exp(self):
        return self.args[1]


from sympy import Expr, Pow, S, Number, Symbol, sympify
from sympy.physics.quantum import *
from sympy.physics.qmul import QMul
from sympy.physics.qadd import QAdd
from sympy.physics.qoperations import QAssocOp
from sympy.core.decorators import call_highest_priority

class QPow(Expr):
    """
    A class for the operator: ** (exponent) for quantum objects.
    """

    __slots__ = ['evaluates', 'hilbert_space']

    def __new__(cls, base, exp):
        base = sympify(base)
        exp = sympify(exp)
        return cls._rules_QPow(base, exp)
        
    @classmethod
    def _rules_QPow(cls, base, exp):
        from sympy.physics.quantum import StateBase, Operator, OuterProduct, KetBase, BraBase, OuterProduct, InnerProduct
        if not isinstance(base, (Operator, OuterProduct, StateBase, QAssocOp, QPow)):
            if not isinstance(exp,(Operator, OuterProduct, StateBase, QAssocOp, QPow)):
                return Pow(base, exp)
            elif issubclass(exp.evaluates, (Operator, OuterProduct, InnerProduct)):
                ret = Expr.__new__(cls, base, exp)
                ret.hilbert_space = exp.hilbert_space
                ret.evaluates = exp.evaluates 
                return ret                        
        elif not isinstance(exp,(Operator, OuterProduct, StateBase, QAssocOp, QPow)):
            if issubclass(base.evaluates, (Operator, OuterProduct, InnerProduct)):
                if exp == S.Zero:
                    return S.One
                elif exp == S.One:
                    return base
                    
                ret = Expr.__new__(cls, base, exp)
                ret.hilbert_space = base.hilbert_space
                ret.evaluates = base.evaluates 
                return ret
        elif issubclass(exp.evaluates, InnerProduct) and issubclass(base.evaluates, (OuterProduct, InnerProduct, Operator)):
                ret = Expr.__new__(cls, base, exp)
                ret.hilbert_space = base.hilbert_space
                ret.evaluates = base.evaluates 
                return ret
        elif issubclass(base.evaluates, InnerProduct) and issubclass(exp.evaluates, (OuterProduct, InnerProduct, Operator)):           
                ret = Expr.__new__(cls, base, exp)
                ret.hilbert_space = exp.hilbert_space
                ret.evaluates = exp.evaluates 
                return ret

        #make a pretty error message if you have left everything a mess   
        if hasattr(exp, 'evaluates'):
            expname = exp.evaluates.__name__
        else:
            expname = exp.__class__.__name__
            
        if hasattr(base, 'evaluates'):
            basename = base.evaluates.__name__
        else:
            basename = base.__class__.__name__
            
        raise Exception("Can't do (%s)**(%s)" % (basename, expname))

    @property
    def base(self):
        return self.args[0]

    @property
    def exp(self):
        return self.args[1]

    @call_highest_priority('__rmul__')
    def __mul__(self, other):
        from sympy.physics.qmul import QMul
        return QMul(self, other)

    @call_highest_priority('__mul__')
    def __rmul__(self, other):
        from sympy.physics.qmul import QMul
        return QMul(other, self)

    @call_highest_priority('__radd__')
    def __add__(self, other):
        from sympy.physics.qadd import QAdd    
        return QAdd(self, other)

    @call_highest_priority('__add__')
    def __radd__(self, other):
        from sympy.physics.qadd import QAdd       
        return QAdd(other, self)

    @call_highest_priority('__rpow__')
    def __pow__(self, other):
        from sympy.physics.qpow import QPow
        return QPow(self, other)

    @call_highest_priority('__pow__')
    def __rpow__(self, other):
        from sympy.physics.qpow import QPow
        return QPow(other, self)


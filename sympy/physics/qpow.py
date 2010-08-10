from sympy import Expr, Pow, S, Number, Symbol, sympify
from sympy.physics.quantum import *
from sympy.physics.qmul import QMul
from sympy.physics.qadd import QAdd
from sympy.physics.qoperations import QAssocOp
from sympy.core.decorators import call_highest_priority
from sympy.physics.quantumbasic import QuantumError, QuantumBasic
from sympy.printing.str import sstr

class QPow(QuantumBasic):
    """
    A class for the operator: ** (exponent) for quantum objects.
    """

    def __new__(cls, base, exp):
        base = sympify(base)
        exp = sympify(exp)
        return cls._rules_QPow(base, exp)
        
    @classmethod
    def _rules_QPow(cls, base, exp):
        if not isinstance(base, QuantumBasic):
            if not isinstance(exp, QuantumBasic):
                return Pow(base, exp)
            elif issubclass(exp.evaluates, (Operator, OuterProduct, InnerProduct)):
                ret = Expr.__new__(cls, base, exp)
                ret.hilbert_space = exp.hilbert_space
                ret.evaluates = exp.evaluates 
                return ret                        
        elif not isinstance(exp, QuantumBasic):
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

    def _eval_dagger(self):
        return QPow(Dagger(self.base), Dagger(self.exp))

    def _sympystr(self, printer, *args):
        return '(' + sstr(self.base) + ')' '**' + sstr(self.exp)

    def _pretty(self, printer, *args):
        return printer._print(self.args[0], *args)**printer._print(self.args[1], *args)
         
    


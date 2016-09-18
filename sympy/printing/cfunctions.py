import math
from sympy.core import S
from sympy.core.symbol import Dummy
from sympy.core.function import Lambda
from sympy.functions.elementary.exponential import exp, log
from sympy.core.function import Function, ArgumentIndexError

_dummy = Dummy()

class expm1(Function):
    nargs = 1
    
    _imp_ = Lambda(_dummy, exp(_dummy) - S.One)

    def _eval_expand_func(self, **hints):
        return exp(self.args[0]) - S.One
    
    def _eval_rewrite_as_exp(self, arg):
        return exp(arg) - S.One
    
    _eval_rewrite_as_tractable = _eval_rewrite_as_exp
    
    @classmethod
    def eval(cls, arg):
        exp_arg = exp.eval(arg)
        if exp_arg is not None:
            return exp_arg - S.One
    
    def _eval_is_real(self):
        return self.args[0].is_real
    
    def _eval_is_finite(self):
        return self.args[0].is_finite
    

class log1p(Function):
    nargs = 1
    
    _imp_ = Lambda(_dummy, log(_dummy + S.One))

    def _eval_expand_func(self, **hints):
        return log(self.args[0] + S.One)
    
    def _eval_rewrite_as_log(self, arg):
        return log(arg + S.One)
    
    _eval_rewrite_as_tractable = _eval_rewrite_as_log
    
    @classmethod
    def eval(cls, arg):
        print(arg, arg.is_Float)
        if not arg.is_Float:  # not safe to add 1 to Float
            return log.eval(arg + S.One)
        # if arg.is_Float and arg.is_real and arg > -1:
        #     return S(math.log1p(arg))
        # else:
    
    def _eval_is_real(self):
        return (self.args[0] + S.One).is_nonnegative
    
    def _eval_is_finite(self):
        if (self.args[0] + S.One).is_zero:
            return False
        return self.args[0].is_finite
    
    def _eval_is_positive(self):
        return self.args[0].is_positive
    
    def _eval_is_zero(self):
        return self.args[0].is_zero
    
    def _eval_is_nonnegative(self):
        return self.args[0].is_nonnegative

_Two = S.One*2
class exp2(Function):
    nargs = 1
    _imp_ = Lambda(_dummy, _Two**_dummy)

    def _eval_expand_func(self, **hints):
        return _Two**self.args[0]

    @classmethod
    def eval(cls, arg):
        if arg.is_number:
            return _Two**arg

    
class log2(Function):
    nargs = 1
    _imp_ = Lambda(_dummy, log(_dummy)/log(_Two))
    
    @classmethod
    def eval(cls, arg):
        if arg.is_number:
            result = log.eval(arg, base=_Two)
            if result.is_Atom:
                return result

    def _eval_expand_func(self, **hints):
        return log(self.args[0])/log(_Two)

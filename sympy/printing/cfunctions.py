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
        return log(arg + 1)
    
    _eval_rewrite_as_tractable = _eval_rewrite_as_log
    
    @classmethod
    def eval(cls, arg):
        return log.eval(arg + S.One)
    
    def _eval_is_real(self):
        return (self.args[0] - S.One).is_positive
    
    def _eval_is_finite(self):
        if (self.args[0] - S.One).is_zero:
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

    
class log2(Function):
    nargs = 1
    _imp_ = Lambda(_dummy, log(_dummy)/log(_Two))
    
    def _eval_expand_func(self, **hints):
        return log(self.args[0])/log(_Two)

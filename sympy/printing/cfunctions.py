"""
Functions with corresponding implementations in C.

The functions defined in this module allows the user to express functions such as ``expm1``
as a SymPy function for symbolic manipulation.

"""

import math
from sympy.core import S, Rational
from sympy.core.function import Lambda
from sympy.core.power import Pow
from sympy.core.symbol import Dummy
from sympy.functions import gamma, exp, log
from sympy.core.function import Function, ArgumentIndexError

_dummy = Dummy()

class expm1(Function):
    nargs = 1

    _imp_ = Lambda(_dummy, exp(_dummy) - S.One)

    def fdiff(self, argindex=1):
        """
        Returns the first derivative of this function.
        """
        if argindex == 1:
            return exp(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

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


    def fdiff(self, argindex=1):
        """
        Returns the first derivative of this function.
        """
        if argindex == 1:
            return S.One/(self.args[0] + S.One)
        else:
            raise ArgumentIndexError(self, argindex)


    def _eval_expand_func(self, **hints):
        return log(self.args[0] + S.One)

    def _eval_rewrite_as_log(self, arg):
        return log(arg + S.One)

    _eval_rewrite_as_tractable = _eval_rewrite_as_log

    @classmethod
    def eval(cls, arg):
        if not arg.is_Float:  # not safe to add 1 to Float
            return log.eval(arg + S.One)
        elif arg.is_number:
            return log.eval(Rational(arg) + S.One)

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
    _imp_ = Lambda(_dummy, Pow(_Two, _dummy))


    def fdiff(self, argindex=1):
        """
        Returns the first derivative of this function.
        """
        if argindex == 1:
            return self*log(_Two)
        else:
            raise ArgumentIndexError(self, argindex)


    def _eval_expand_func(self, **hints):
        return Pow(_Two, self.args[0])

    @classmethod
    def eval(cls, arg):
        if arg.is_number:
            return Pow(_Two, arg)


class log2(Function):
    nargs = 1
    _imp_ = Lambda(_dummy, log(_dummy)/log(_Two))


    def fdiff(self, argindex=1):
        """
        Returns the first derivative of this function.
        """
        if argindex == 1:
            return S.One/(log(_Two)*self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)


    @classmethod
    def eval(cls, arg):
        if arg.is_number:
            result = log.eval(arg, base=_Two)
            if result.is_Atom:
                return result
        elif arg.is_Pow and arg.base == _Two:
            return arg.exp

    def _eval_expand_func(self, **hints):
        return log(self.args[0])/log(_Two)


_dummy0, _dummy1, _dummy2 = Dummy(), Dummy(), Dummy()


class fma(Function):
    nargs = 3
    _imp_ = Lambda((_dummy0, _dummy1, _dummy2), _dummy0*_dummy1 + _dummy2)


class log10(Function):
    nargs = 1
    _imp_ = Lambda(_dummy, log(_dummy)/log(10))


class Cbrt(Function):  # 'cbrt' already defined in sympy.functions.elementary.miscellaneous
    nargs = 1
    _imp_ = Lambda(_dummy, Pow(_dummy, Rational(1, 3)))


class hypot(Function):
    nargs = 2
    _imp_ = Lambda((_dummy0, _dummy1), Pow(Pow(_dummy0, _Two) + Pow(_dummy1, _Two), Rational(1, _Two)))


class lgamma(Function):
    nargs = 1
    _imp_ = Lambda(_dummy, log(gamma(_dummy)))

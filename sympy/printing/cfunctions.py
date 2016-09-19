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

def _expm1(x):
    return exp(x) - S.One


class expm1(Function):
    nargs = 1

    def fdiff(self, argindex=1):
        """
        Returns the first derivative of this function.
        """
        if argindex == 1:
            return exp(*self.args)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_expand_func(self, **hints):
        return _expm1(*self.args)

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


def _log1p(x):
    return log(x + S.One)


class log1p(Function):
    nargs = 1


    def fdiff(self, argindex=1):
        """
        Returns the first derivative of this function.
        """
        if argindex == 1:
            return S.One/(self.args[0] + S.One)
        else:
            raise ArgumentIndexError(self, argindex)


    def _eval_expand_func(self, **hints):
        return _log1p(*self.args)

    def _eval_rewrite_as_log(self, arg):
        return _log1p(arg)

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

_Two = S(2)

def _exp2(x):
    return Pow(_Two, x)

class exp2(Function):
    nargs = 1


    def fdiff(self, argindex=1):
        """
        Returns the first derivative of this function.
        """
        if argindex == 1:
            return self*log(_Two)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_Pow(self, arg):
        return _exp2(arg)

    _eval_rewrite_as_tractable = _eval_rewrite_as_Pow

    def _eval_expand_func(self, **hints):
        return _exp2(*self.args)

    @classmethod
    def eval(cls, arg):
        if arg.is_number:
            return _exp2(arg)


def _log2(x):
    return log(x)/log(_Two)


class log2(Function):
    nargs = 1

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
        return _log2(*self.args)

    def _eval_rewrite_as_log(self, arg):
        return _log2(arg)

    _eval_rewrite_as_tractable = _eval_rewrite_as_log


_dummy0, _dummy1, _dummy2 = Dummy(), Dummy(), Dummy()


def _fma(x, y, z):
    return x*y + z


class fma(Function):
    nargs = 3

    def _eval_expand_func(self, **hints):
        return _fma(*self.args)

    def _eval_rewrite_as_tractable(self, arg):
        return _fma(arg)


_Ten = S(10)


def _log10(x):
    return log(x)/log(_Ten)


class log10(Function):
    nargs = 1

    def fdiff(self, argindex=1):
        """
        Returns the first derivative of this function.
        """
        if argindex == 1:
            return S.One/(log(_Ten)*self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)


    @classmethod
    def eval(cls, arg):
        if arg.is_number:
            result = log.eval(arg, base=_Ten)
            if result.is_Atom:
                return result
        elif arg.is_Pow and arg.base == _Ten:
            return arg.exp

    def _eval_expand_func(self, **hints):
        return _log10(*self.args)

    def _eval_rewrite_as_log(self, arg):
        return _log10(arg)

    _eval_rewrite_as_tractable = _eval_rewrite_as_log


def _Cbrt(x):
    return Pow(x, Rational(1, 3))


class Cbrt(Function):  # 'cbrt' already defined in sympy.functions.elementary.miscellaneous
    nargs = 1

    def _eval_expand_func(self, **hints):
        return _Cbrt(*self.args)

    def _eval_rewrite_as_Pow(self, arg):
        return _Cbrt(arg)

    _eval_rewrite_as_tractable = _eval_rewrite_as_Pow



def _hypot(x, y):
    return Pow(Pow(x, 2) + Pow(y, 2), Rational(1, 2))

class hypot(Function):
    nargs = 2

    def _eval_expand_func(self, **hints):
        return _hypot(*self.args)

    def _eval_rewrite_as_Pow(self, arg):
        return _hypot(arg)

    _eval_rewrite_as_tractable = _eval_rewrite_as_Pow



def _lgamma(x):
    return log(gamma(x))


class lgamma(Function):
    nargs = 1

    def _eval_expand_func(self, **hints):
        return _lgamma(*self.args)

    def _eval_rewrite_as_gamma(self, arg):
        return _lgamma(arg)

    _eval_rewrite_as_tractable = _eval_rewrite_as_gamma

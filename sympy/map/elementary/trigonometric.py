from sympy import (
    S, pi, Interval, Rational
)
from sympy.core.logic import fuzzy_not
from sympy.map import Map, AppliedMap, isappliedmap

__all__ = [
    "Sin", "Cos", "Tan",
]

# Details are not implemented.
# This module should import nothing from trigonometric.py.
# Instead, all should be implemented under sympy.map.

###############################################################################
############################### Basic classes #################################
###############################################################################

class TrigonometricMap(Map):
    """Base class for trigonometric functions. """

    _singularities = (S.ComplexInfinity,)

    def __new__(cls, domain=S.Complexes, **kwargs):
        return super().__new__(cls, domain)

    @property
    def domain(self):
        return self.args[0]

    @property
    def codomain(self):
        domain = self.domain
        if not domain.is_subset(S.Reals):
            return S.ComplexesField
        return S.RealsField

    def __call__(self, *args, evaluate=False, **kwargs):
        return AppliedTrigonometricMap(self, args, evaluate=evaluate, **kwargs)

class AppliedTrigonometricMap(AppliedMap):
    """Base class for applied trigonometric functions. """

    @property
    def _singularities(self):
        return self.map._singularities

    def _eval_is_rational(self):
        s = self.map(*self.arguments, evaluate=True)
        if (s.func == self.func) and (s.map == self.map):
            x = s.arguments[0]
            if x.is_rational and fuzzy_not(x.is_zero):
                return False
        else:
            return s.is_rational

    def _eval_is_algebraic(self):
        s = self.map(*self.arguments, evaluate=True)
        if (s.func == self.func) and (s.map == self.map):
            x = s.arguments[0]
            if fuzzy_not(x.is_zero) and x.is_algebraic:
                return False
            pi_coeff = _pi_coeff(x)
            if pi_coeff is not None and pi_coeff.is_rational:
                return True
        else:
            return s.is_algebraic

    def _eval_expand_complex(self, deep=True, **hints):
        re_part, im_part = self.as_real_imag(deep=deep, **hints)
        return re_part + im_part*S.ImaginaryUnit

    def _as_real_imag(self, deep=True, **hints):
        x = self.arguments[0]
        if x.is_extended_real:
            if deep:
                hints['complex'] = False
                return (x.expand(deep, **hints), S.Zero)
            else:
                return (x, S.Zero)
        if deep:
            re, im = x.expand(deep, **hints).as_real_imag()
        else:
            re, im = x.as_real_imag()
        return (re, im)

    def period(self, symbol=None):
        return self._period(self.map.period, symbol)

    def _period(self, general_period, symbol=None):

        x = self.arguments[0]
        f = expand_mul(x)
        if symbol is None:
            symbol = tuple(f.free_symbols)[0]

        if not f.has(symbol):
            return S.Zero

        if f == symbol:
            return general_period

        if symbol in f.free_symbols:
            if f.is_Mul:
                g, h = f.as_independent(symbol)
                if h == symbol:
                    return general_period/abs(g)

            if f.is_Add:
                a, h = f.as_independent(symbol)
                g, h = h.as_independent(symbol, as_Add=False)
                if h == symbol:
                    return general_period/abs(g)

        raise NotImplementedError("Use the periodicity function instead.")

###############################################################################
########################## TRIGONOMETRIC FUNCTIONS ############################
###############################################################################

class Sin(TrigonometricMap):
    """
    The sine function.

    """
    name = 'sin'
    latex_name = '\\sin'
    period = 2*pi

    def _eval_range(self):
        domain = self.domain
        if not domain.is_subset(S.Reals):
            return S.Complexes
        return Interval(-1, 1)

class Cos(TrigonometricMap):
    """
    The cosine function.

    """
    name = 'cos'
    latex_name = '\\cos'
    period = 2*pi

    def _eval_range(self):
        domain = self.domain
        if not domain.is_subset(S.Reals):
            return S.Complexes
        return Interval(-1, 1)

class Tan(TrigonometricMap):
    """
    The tangent function

    """
    name = 'tan'
    latex_name = '\\tan'
    period = pi

    def _eval_range(self):
        domain = self.domain
        if not domain.is_subset(S.Reals):
            return S.Complexes
        return S.Reals

from sympy.core.function import expand_mul
from sympy.functions.elementary.trigonometric import _pi_coeff, _peeloff_pi

from sympy import (
    S, pi, Interval,
)
from sympy.map import Map, AppliedMap,

__all__ = [
    "Sine", "Cosine", "Tangent", "Cotangent",
    "Secant", "Cosecant",
]

# Details are not implemented.
# This module should import nothing from trigonometric.py.
# Instead, all should be implemented under sympy.map.

# NOTE: every functions will not be evaluated by default.
# NOTE: eschew old assumptions

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

    def _eval_period(self, applied, symbol=None):

        general_period = self.period

        x = applied.arguments[0]
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

class AppliedTrigonometricMap(AppliedMap):
    """Base class for applied trigonometric functions. """

    @property
    def _singularities(self):
        return self.map._singularities

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
        return self.map._eval_period(self, symbol)

###############################################################################
########################## TRIGONOMETRIC FUNCTIONS ############################
###############################################################################

class Sine(TrigonometricMap):
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

class Cosine(TrigonometricMap):
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

class Tangent(TrigonometricMap):
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

class Cotangent(TrigonometricMap):
    """
    The cotangent function

    """
    name = 'cot'
    latex_name = '\\cot'
    period = pi

    def _eval_range(self):
        domain = self.domain
        if not domain.is_subset(S.Reals):
            return S.Complexes
        return S.Reals

###############################################################################
##################### RECIPROCAL TRIGONOMETRIC FUNCTIONS ######################
###############################################################################

class ReciprocalTrigonometricMap(TrigonometricMap):
    """Base class for reciprocal functions of trigonometric functions. """

    @property
    def _reciprocal_of(self):
        # mandatory, to be defined in subclass
        return

    def _eval_period(self, applied, symbol):
        f = expand_mul(applied.arguments[0])
        return self._reciprocal_of(f, evaluate=True).period(symbol)

class Secant(ReciprocalTrigonometricMap):
    """
    The secant function

    """

    @property
    def _reciprocal_of(self):
        return Cosine(self.domain)

    def _eval_range(self):
        domain = self.domain
        if not domain.is_subset(S.Reals):
            return S.Complexes
        return S.Reals - Interval(-1, 1)

class Cosecant(ReciprocalTrigonometricMap):
    """
    The secant function

    """

    @property
    def _reciprocal_of(self):
        return Sine(self.domain)

    def _eval_range(self):
        domain = self.domain
        if not domain.is_subset(S.Reals):
            return S.Complexes
        return S.Reals - Interval(-1, 1)

from sympy.core.function import expand_mul

from sympy import S, pi, Interval, Add, Rational
from sympy.core.sympify import _sympify
from sympy.core.logic import fuzzy_not
from sympy.map import Map, AppliedMap, isappliedmap
# Belows will be substituted to Map object
from sympy.functions.elementary.exponential import log, exp
from sympy.functions.elementary.hyperbolic import cosh, sinh

__all__ = [
    'Sine', 'Cosine', 'Tangent', 'Cotangent',
    'Secant', 'Cosecant',
    'ArcSine', 'ArcCosine', 'ArcTangent', 'ArcCotangent',
    'ArcSecant', 'ArcCosecant',
]

def _peeloff_pi(arg):
    """
    Split ARG into two parts, a "rest" and a multiple of pi/2.
    This assumes ARG to be an Add.
    The multiple of pi returned in the second position is always a Rational.

    Examples
    ========

    >>> from sympy.map.elementary.trigonometric import _peeloff_pi as peel
    >>> from sympy import pi
    >>> from sympy.abc import x, y
    >>> peel(x + pi/2)
    (x, pi/2)
    >>> peel(x + 2*pi/3 + pi*y)
    (x + pi*y + pi/6, pi/2)

    """
    pi_coeff = S.Zero
    rest_terms = []
    for a in Add.make_args(arg):
        K = a.coeff(S.Pi)
        if K and K.is_rational:
            pi_coeff += K
        else:
            rest_terms.append(a)

    if pi_coeff is S.Zero:
        return arg, S.Zero

    m1 = (pi_coeff % S.Half)*S.Pi
    m2 = pi_coeff*S.Pi - m1
    final_coeff = m2 / S.Pi
    if final_coeff.is_integer or ((2*final_coeff).is_integer
        and final_coeff.is_even is False):
            return Add(*(rest_terms + [m1])), m2
    return arg, S.Zero

def _pi_coeff(arg, cycles=1):
    """
    When arg is a Number times pi (e.g. 3*pi/2) then return the Number
    normalized to be in the range [0, 2], else None.

    When an even multiple of pi is encountered, if it is multiplying
    something with known parity then the multiple is returned as 0 otherwise
    as 2.

    Examples
    ========

    >>> from sympy.map.elementary.trigonometric import _pi_coeff as coeff
    >>> from sympy import pi, Dummy
    >>> from sympy.abc import x
    >>> coeff(3*x*pi)
    3*x
    >>> coeff(11*pi/7)
    11/7
    >>> coeff(-11*pi/7)
    3/7
    >>> coeff(4*pi)
    0
    >>> coeff(5*pi)
    1
    >>> coeff(5.0*pi)
    1
    >>> coeff(5.5*pi)
    3/2
    >>> coeff(2 + pi)

    >>> coeff(2*Dummy(integer=True)*pi)
    2
    >>> coeff(2*Dummy(even=True)*pi)
    0

    """
    arg = _sympify(arg)
    if arg is S.Pi:
        return S.One
    elif not arg:
        return S.Zero
    elif arg.is_Mul:
        cx = arg.coeff(S.Pi)
        if cx:
            c, x = cx.as_coeff_Mul()  # pi is not included as coeff
            if c.is_Float:
                # recast exact binary fractions to Rationals
                f = abs(c) % 1
                if f != 0:
                    p = -int(round(log(f, 2).evalf()))
                    m = 2**p
                    cm = c*m
                    i = int(cm)
                    if i == cm:
                        c = Rational(i, m)
                        cx = c*x
                else:
                    c = Rational(int(c))
                    cx = c*x
            if x.is_integer:
                c2 = c % 2
                if c2 == 1:
                    return x
                elif not c2:
                    if x.is_even is not None:  # known parity
                        return S.Zero
                    return S(2)
                else:
                    return c2*x
            return cx
    elif arg.is_zero:
        return S.Zero

###############################################################################
############################### Basic classes #################################
###############################################################################

class TrigonometricMap(Map):
    """Base class for trigonometric functions. """

    singularities = (S.ComplexInfinity,)

    reciprocal_of = None       # mandatory, to be defined in subclass

    # optional, to be defined in subclasses:
    is_even = None  # type: FuzzyBool
    is_odd = None  # type: FuzzyBool

    def __new__(cls, domain=S.Complexes, **kwargs):
        return super().__new__(cls, domain)

    @property
    def domain(self):
        return self.args[0]

    @property
    def codomain(self):
        domain = self.domain
        if not domain.is_subset(S.Reals):
            return S.Complexes
        return S.Reals

    def __call__(self, *args, evaluate=False, **kwargs):
        return AppliedTrigonometricMap(self, args, evaluate=evaluate, **kwargs)

    def _applied_is_rational(self, expr):
        s = expr.doit(deep=False)
        if isappliedmap(s, self):
            args = s.arguments
            if args[0].is_rational and fuzzy_not(args[0].is_zero):
                return False
        else:
            return s.is_rational

    def _appled_is_algebraic(self, expr):
        s = expr.doit(deep=False)
        if isappliedmap(s, self):
            args = s.arguments
            if fuzzy_not(args[0].is_zero) and args[0].is_algebraic:
                return False
            pi_coeff = _pi_coeff(args[0])
            if pi_coeff is not None and pi_coeff.is_rational:
                return True
        else:
            return s.is_algebraic

    def _applied_as_real_imag(self, expr, **hints):
        deep = kwargs.get('deep', True)
        args = expr.arguments
        if args[0].is_extended_real:
            if deep:
                hints['complex'] = False
                return (args[0].expand(deep, **hints), S.Zero)
            else:
                return (args[0], S.Zero)
        if deep:
            re, im = args[0].expand(deep, **hints).as_real_imag()
        else:
            re, im = args[0].as_real_imag()
        return re, im

    def _applied_period(self, expr, symbol):
        general_period = self.period
        args = expr.arguments
        f = expand_mul(args[0])
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

    unbranched = True

    @property
    def _singularities(self):
        return self.map.singularities

    def _eval_is_rational(self):
        return self.map._applied_is_rational(self)

    def _eval_is_algebraic(self):
        return self.map._applied_is_algebraic(self)

    def _eval_expand_complex(self, deep=True, **hints):
        hints.update(deep=deep)
        re_part, im_part = self.map._applied_expand_complex(self, **hints)
        return re_part + im_part*S.ImaginaryUnit

    def as_real_imag(self, deep=True, **hints):
        hints.update(deep=deep)
        return self.map._applied_as_real_imag(self, **hints)

    def period(self, symbol=None):
        return self.map._applied_period(self, symbol)

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

    is_even = False
    is_odd = True

    @property
    def reciprocal_of(self):
        return Secant(self.domain)

    def _eval_range(self):
        domain = self.domain
        if not domain.is_subset(S.Reals):
            return S.Complexes
        return Interval(-1, 1)

    def _applied_as_real_imag(self, expr, **hints):
        cos = Cosine(self.domain)
        re, im = super()._applied_as_real_imag(self, expr, **hints)
        return (
            self(re, evaluate=True)*cosh(im, evaluate=True),
            cos(re, evaluate=True)*sinh(im, evaluate=True)
        )

class Cosine(TrigonometricMap):
    """
    The cosine function.

    """
    name = 'cos'
    latex_name = '\\cos'
    period = 2*pi

    is_even = True
    is_odd = False

    @property
    def reciprocal_of(self):
        return Secant(self.domain)

    def _eval_range(self):
        domain = self.domain
        if not domain.is_subset(S.Reals):
            return S.Complexes
        return Interval(-1, 1)

    def _applied_as_real_imag(self, expr, **hints):
        sin = Sine(self.domain)
        re, im = super()._applied_as_real_imag(self, expr, **hints)
        return (
            self(re, evaluate=True)*cosh(im, evaluate=True),
            -sin(re, evaluate=True)*sinh(im, evaluate=True)
        )

class Tangent(TrigonometricMap):
    """
    The tangent function

    """
    name = 'tan'
    latex_name = '\\tan'
    period = pi

    is_even = False
    is_odd = True

    @property
    def reciprocal_of(self):
        return Cotangent(self.domain)

    def _eval_range(self):
        domain = self.domain
        if not domain.is_subset(S.Reals):
            return S.Complexes
        return S.Reals

    def _applied_as_real_imag(self, expr, **hints):
        cos = Cosine(self.domain)
        sin = Sine(self.domain)
        re, im = super()._applied_as_real_imag(self, expr, **hints)
        if im:
            denom = cos(2*re, evaluate=True) + cosh(2*im, evaluate=True)
            return (
                sin(2*re, evaluate=True)/denom,
                sinh(2*im, evaluate=True)/denom
            )
        else:
            return (
                self(re, evaluate=True),
                S.Zero
            )

class Cotangent(TrigonometricMap):
    """
    The cotangent function

    """
    name = 'cot'
    latex_name = '\\cot'
    period = pi

    is_even = False
    is_odd = True

    @property
    def reciprocal_of(self):
        return Tangent(self.domain)

    def _eval_range(self):
        domain = self.domain
        if not domain.is_subset(S.Reals):
            return S.Complexes
        return S.Reals

    def _applied_as_real_imag(self, expr, **hints):
        cos = Cosine(self.domain)
        sin = Sine(self.domain)
        re, im = super()._applied_as_real_imag(self, expr, **hints)
        if im:
            denom = cos(2*re, evaluate=True) - cosh(2*im, evaluate=True)
            return (
                -sin(2*re, evaluate=True)/denom,
                sinh(2*im, evaluate=True)/denom
            )
        else:
            return (
                self(re, evaluate=True),
                S.Zero
            )

###############################################################################
##################### RECIPROCAL TRIGONOMETRIC FUNCTIONS ######################
###############################################################################

class ReciprocalTrigonometricMap(TrigonometricMap):
    """Base class for reciprocal functions of trigonometric functions. """

    def _applied_period(self, expr, symbol):
        f = expand_mul(expr.arguments[0])
        return self.reciprocal_of(f, evaluate=True).period(symbol)

    def _applied_as_real_imag(self, expr, **hints):
        arg = expr.arguments[0]
        denom = self.reciprocal_of(arg, evaluate=True)
        frac = 1/denom
        return frac.as_real_imag(**hints)

class Secant(ReciprocalTrigonometricMap):
    """
    The secant function

    """

    @property
    def reciprocal_of(self):
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
    def reciprocal_of(self):
        return Sine(self.domain)

    def _eval_range(self):
        domain = self.domain
        if not domain.is_subset(S.Reals):
            return S.Complexes
        return S.Reals - Interval(-1, 1)

###############################################################################
########################### TRIGONOMETRIC INVERSES ############################
###############################################################################

class InverseTrigonometricMap(Map):

    @staticmethod
    def _asin_table():
        # Only keys with could_extract_minus_sign() == False
        # are actually needed.
        return {
            sqrt(3)/2: S.Pi/3,
            sqrt(2)/2: S.Pi/4,
            1/sqrt(2): S.Pi/4,
            sqrt((5 - sqrt(5))/8): S.Pi/5,
            sqrt(2)*sqrt(5 - sqrt(5))/4: S.Pi/5,
            sqrt((5 + sqrt(5))/8): S.Pi*Rational(2, 5),
            sqrt(2)*sqrt(5 + sqrt(5))/4: S.Pi*Rational(2, 5),
            S.Half: S.Pi/6,
            sqrt(2 - sqrt(2))/2: S.Pi/8,
            sqrt(S.Half - sqrt(2)/4): S.Pi/8,
            sqrt(2 + sqrt(2))/2: S.Pi*Rational(3, 8),
            sqrt(S.Half + sqrt(2)/4): S.Pi*Rational(3, 8),
            (sqrt(5) - 1)/4: S.Pi/10,
            (1 - sqrt(5))/4: -S.Pi/10,
            (sqrt(5) + 1)/4: S.Pi*Rational(3, 10),
            sqrt(6)/4 - sqrt(2)/4: S.Pi/12,
            -sqrt(6)/4 + sqrt(2)/4: -S.Pi/12,
            (sqrt(3) - 1)/sqrt(8): S.Pi/12,
            (1 - sqrt(3))/sqrt(8): -S.Pi/12,
            sqrt(6)/4 + sqrt(2)/4: S.Pi*Rational(5, 12),
            (1 + sqrt(3))/sqrt(8): S.Pi*Rational(5, 12)
        }

    @staticmethod
    def _atan_table():
        # Only keys with could_extract_minus_sign() == False
        # are actually needed.
        return {
            sqrt(3)/3: S.Pi/6,
            1/sqrt(3): S.Pi/6,
            sqrt(3): S.Pi/3,
            sqrt(2) - 1: S.Pi/8,
            1 - sqrt(2): -S.Pi/8,
            1 + sqrt(2): S.Pi*Rational(3, 8),
            sqrt(5 - 2*sqrt(5)): S.Pi/5,
            sqrt(5 + 2*sqrt(5)): S.Pi*Rational(2, 5),
            sqrt(1 - 2*sqrt(5)/5): S.Pi/10,
            sqrt(1 + 2*sqrt(5)/5): S.Pi*Rational(3, 10),
            2 - sqrt(3): S.Pi/12,
            -2 + sqrt(3): -S.Pi/12,
            2 + sqrt(3): S.Pi*Rational(5, 12)
        }

    @staticmethod
    def _acsc_table():
        # Keys for which could_extract_minus_sign()
        # will obviously return True are omitted.
        return {
            2*sqrt(3)/3: S.Pi/3,
            sqrt(2): S.Pi/4,
            sqrt(2 + 2*sqrt(5)/5): S.Pi/5,
            1/sqrt(Rational(5, 8) - sqrt(5)/8): S.Pi/5,
            sqrt(2 - 2*sqrt(5)/5): S.Pi*Rational(2, 5),
            1/sqrt(Rational(5, 8) + sqrt(5)/8): S.Pi*Rational(2, 5),
            2: S.Pi/6,
            sqrt(4 + 2*sqrt(2)): S.Pi/8,
            2/sqrt(2 - sqrt(2)): S.Pi/8,
            sqrt(4 - 2*sqrt(2)): S.Pi*Rational(3, 8),
            2/sqrt(2 + sqrt(2)): S.Pi*Rational(3, 8),
            1 + sqrt(5): S.Pi/10,
            sqrt(5) - 1: S.Pi*Rational(3, 10),
            -(sqrt(5) - 1): S.Pi*Rational(-3, 10),
            sqrt(6) + sqrt(2): S.Pi/12,
            sqrt(6) - sqrt(2): S.Pi*Rational(5, 12),
            -(sqrt(6) - sqrt(2)): S.Pi*Rational(-5, 12)
        }

class ArcSine(InverseTrigonometricMap):
    pass

class ArcCosine(InverseTrigonometricMap):
    pass

class ArcTangent(InverseTrigonometricMap):
    pass

class ArcCotangent(InverseTrigonometricMap):
    pass

class ArcSecant(InverseTrigonometricMap):
    pass

class ArcCosecant(InverseTrigonometricMap):
    pass

from sympy.core.function import expand_mul

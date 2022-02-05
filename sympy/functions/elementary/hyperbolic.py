from sympy.core.logic import FuzzyBool

from sympy.core import S, sympify, cacheit, pi, I, Rational
from sympy.core.add import Add
from sympy.core.function import Function, ArgumentIndexError
from sympy.functions.combinatorial.factorials import factorial, RisingFactorial
from sympy.functions.elementary.exponential import exp, log, match_real_imag
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.integers import floor

from sympy.core.logic import fuzzy_or, fuzzy_and

def _rewrite_hyperbolics_as_exp(expr):
    expr = sympify(expr)
    return expr.xreplace({h: h.rewrite(exp)
        for h in expr.atoms(HyperbolicFunction)})


###############################################################################
########################### HYPERBOLIC FUNCTIONS ##############################
###############################################################################


class HyperbolicFunction(Function):
    """
    Base class for hyperbolic functions.

    See Also
    ========

    sinh, cosh, tanh, coth
    """

    unbranched = True


def _peeloff_ipi(arg):
    r"""
    Split ARG into two parts, a "rest" and a multiple of $I\pi$.
    This assumes ARG to be an ``Add``.
    The multiple of $I\pi$ returned in the second position is always a ``Rational``.

    Examples
    ========

    >>> from sympy.functions.elementary.hyperbolic import _peeloff_ipi as peel
    >>> from sympy import pi, I
    >>> from sympy.abc import x, y
    >>> peel(x + I*pi/2)
    (x, 1/2)
    >>> peel(x + I*2*pi/3 + I*pi*y)
    (x + I*pi*y + I*pi/6, 1/2)
    """
    ipi = S.Pi*S.ImaginaryUnit
    for a in Add.make_args(arg):
        if a == ipi:
            K = S.One
            break
        elif a.is_Mul:
            K, p = a.as_two_terms()
            if p == ipi and K.is_Rational:
                break
    else:
        return arg, S.Zero

    m1 = (K % S.Half)
    m2 = K - m1
    return arg - m2*ipi, m2


class sinh(HyperbolicFunction):
    r"""
    ``sinh(x)`` is the hyperbolic sine of ``x``.

    The hyperbolic sine function is $\frac{e^x - e^{-x}}{2}$.

    Examples
    ========

    >>> from sympy import sinh
    >>> from sympy.abc import x
    >>> sinh(x)
    sinh(x)

    See Also
    ========

    cosh, tanh, asinh
    """

    def fdiff(self, argindex=1):
        """
        Returns the first derivative of this function.
        """
        if argindex == 1:
            return cosh(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        """
        Returns the inverse of this function.
        """
        return asinh

    @classmethod
    def eval(cls, arg):
        from sympy.functions.elementary.trigonometric import sin

        arg = sympify(arg)

        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.Infinity
            elif arg is S.NegativeInfinity:
                return S.NegativeInfinity
            elif arg.is_zero:
                return S.Zero
            elif arg.is_negative:
                return -cls(-arg)
        else:
            if arg is S.ComplexInfinity:
                return S.NaN

            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return S.ImaginaryUnit * sin(i_coeff)
            else:
                if arg.could_extract_minus_sign():
                    return -cls(-arg)

            if arg.is_Add:
                x, m = _peeloff_ipi(arg)
                if m:
                    m = m*S.Pi*S.ImaginaryUnit
                    return sinh(m)*cosh(x) + cosh(m)*sinh(x)

            if arg.is_zero:
                return S.Zero

            if arg.func == asinh:
                return arg.args[0]

            if arg.func == acosh:
                x = arg.args[0]
                return sqrt(x - 1) * sqrt(x + 1)

            if arg.func == atanh:
                x = arg.args[0]
                return x/sqrt(1 - x**2)

            if arg.func == acoth:
                x = arg.args[0]
                return 1/(sqrt(x - 1) * sqrt(x + 1))

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        """
        Returns the next term in the Taylor series expansion.
        """
        if n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)

            if len(previous_terms) > 2:
                p = previous_terms[-2]
                return p * x**2 / (n*(n - 1))
            else:
                return x**(n) / factorial(n)

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def as_real_imag(self, deep=True, **hints):
        """
        Returns this function as a complex coordinate.
        """
        from sympy.functions.elementary.trigonometric import (cos, sin)
        if self.args[0].is_extended_real:
            if deep:
                hints['complex'] = False
                return (self.expand(deep, **hints), S.Zero)
            else:
                return (self, S.Zero)
        if deep:
            re, im = self.args[0].expand(deep, **hints).as_real_imag()
        else:
            re, im = self.args[0].as_real_imag()
        return (sinh(re)*cos(im), cosh(re)*sin(im))

    def _eval_expand_complex(self, deep=True, **hints):
        re_part, im_part = self.as_real_imag(deep=deep, **hints)
        return re_part + im_part*S.ImaginaryUnit

    def _eval_expand_trig(self, deep=True, **hints):
        if deep:
            arg = self.args[0].expand(deep, **hints)
        else:
            arg = self.args[0]
        x = None
        if arg.is_Add: # TODO, implement more if deep stuff here
            x, y = arg.as_two_terms()
        else:
            coeff, terms = arg.as_coeff_Mul(rational=True)
            if coeff is not S.One and coeff.is_Integer and terms is not S.One:
                x = terms
                y = (coeff - 1)*x
        if x is not None:
            return (sinh(x)*cosh(y) + sinh(y)*cosh(x)).expand(trig=True)
        return sinh(arg)

    def _eval_rewrite_as_tractable(self, arg, limitvar=None, **kwargs):
        return (exp(arg) - exp(-arg)) / 2

    def _eval_rewrite_as_exp(self, arg, **kwargs):
        return (exp(arg) - exp(-arg)) / 2

    def _eval_rewrite_as_cosh(self, arg, **kwargs):
        return -S.ImaginaryUnit*cosh(arg + S.Pi*S.ImaginaryUnit/2)

    def _eval_rewrite_as_tanh(self, arg, **kwargs):
        tanh_half = tanh(S.Half*arg)
        return 2*tanh_half/(1 - tanh_half**2)

    def _eval_rewrite_as_coth(self, arg, **kwargs):
        coth_half = coth(S.Half*arg)
        return 2*coth_half/(coth_half**2 - 1)

    def _eval_as_leading_term(self, x, logx=None, cdir=0):
        arg = self.args[0].as_leading_term(x, logx=logx, cdir=cdir)
        arg0 = arg.subs(x, 0)

        if arg0 is S.NaN:
            arg0 = arg.limit(x, 0, dir='-' if cdir.is_negative else '+')
        if arg0.is_zero:
            return arg
        elif arg0.is_finite:
            return self.func(arg0)
        else:
            return self

    def _eval_is_real(self):
        arg = self.args[0]
        if arg.is_real:
            return True

        # if `im` is of the form n*pi
        # else, check if it is a number
        re, im = arg.as_real_imag()
        return (im%pi).is_zero

    def _eval_is_extended_real(self):
        if self.args[0].is_extended_real:
            return True

    def _eval_is_positive(self):
        if self.args[0].is_extended_real:
            return self.args[0].is_positive

    def _eval_is_negative(self):
        if self.args[0].is_extended_real:
            return self.args[0].is_negative

    def _eval_is_finite(self):
        arg = self.args[0]
        return arg.is_finite

    def _eval_is_zero(self):
        rest, ipi_mult = _peeloff_ipi(self.args[0])
        if rest.is_zero:
            return ipi_mult.is_integer


class cosh(HyperbolicFunction):
    r"""
    ``cosh(x)`` is the hyperbolic cosine of ``x``.

    The hyperbolic cosine function is $\frac{e^x + e^{-x}}{2}$.

    Examples
    ========

    >>> from sympy import cosh
    >>> from sympy.abc import x
    >>> cosh(x)
    cosh(x)

    See Also
    ========

    sinh, tanh, acosh
    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return sinh(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        from sympy.functions.elementary.trigonometric import cos
        arg = sympify(arg)

        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.Infinity
            elif arg is S.NegativeInfinity:
                return S.Infinity
            elif arg.is_zero:
                return S.One
            elif arg.is_negative:
                return cls(-arg)
        else:
            if arg is S.ComplexInfinity:
                return S.NaN

            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return cos(i_coeff)
            else:
                if arg.could_extract_minus_sign():
                    return cls(-arg)

            if arg.is_Add:
                x, m = _peeloff_ipi(arg)
                if m:
                    m = m*S.Pi*S.ImaginaryUnit
                    return cosh(m)*cosh(x) + sinh(m)*sinh(x)

            if arg.is_zero:
                return S.One

            if arg.func == asinh:
                return sqrt(1 + arg.args[0]**2)

            if arg.func == acosh:
                return arg.args[0]

            if arg.func == atanh:
                return 1/sqrt(1 - arg.args[0]**2)

            if arg.func == acoth:
                x = arg.args[0]
                return x/(sqrt(x - 1) * sqrt(x + 1))

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0 or n % 2 == 1:
            return S.Zero
        else:
            x = sympify(x)

            if len(previous_terms) > 2:
                p = previous_terms[-2]
                return p * x**2 / (n*(n - 1))
            else:
                return x**(n)/factorial(n)

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def as_real_imag(self, deep=True, **hints):
        from sympy.functions.elementary.trigonometric import (cos, sin)
        if self.args[0].is_extended_real:
            if deep:
                hints['complex'] = False
                return (self.expand(deep, **hints), S.Zero)
            else:
                return (self, S.Zero)
        if deep:
            re, im = self.args[0].expand(deep, **hints).as_real_imag()
        else:
            re, im = self.args[0].as_real_imag()

        return (cosh(re)*cos(im), sinh(re)*sin(im))

    def _eval_expand_complex(self, deep=True, **hints):
        re_part, im_part = self.as_real_imag(deep=deep, **hints)
        return re_part + im_part*S.ImaginaryUnit

    def _eval_expand_trig(self, deep=True, **hints):
        if deep:
            arg = self.args[0].expand(deep, **hints)
        else:
            arg = self.args[0]
        x = None
        if arg.is_Add: # TODO, implement more if deep stuff here
            x, y = arg.as_two_terms()
        else:
            coeff, terms = arg.as_coeff_Mul(rational=True)
            if coeff is not S.One and coeff.is_Integer and terms is not S.One:
                x = terms
                y = (coeff - 1)*x
        if x is not None:
            return (cosh(x)*cosh(y) + sinh(x)*sinh(y)).expand(trig=True)
        return cosh(arg)

    def _eval_rewrite_as_tractable(self, arg, limitvar=None, **kwargs):
        return (exp(arg) + exp(-arg)) / 2

    def _eval_rewrite_as_exp(self, arg, **kwargs):
        return (exp(arg) + exp(-arg)) / 2

    def _eval_rewrite_as_sinh(self, arg, **kwargs):
        return -S.ImaginaryUnit*sinh(arg + S.Pi*S.ImaginaryUnit/2)

    def _eval_rewrite_as_tanh(self, arg, **kwargs):
        tanh_half = tanh(S.Half*arg)**2
        return (1 + tanh_half)/(1 - tanh_half)

    def _eval_rewrite_as_coth(self, arg, **kwargs):
        coth_half = coth(S.Half*arg)**2
        return (coth_half + 1)/(coth_half - 1)

    def _eval_as_leading_term(self, x, logx=None, cdir=0):
        arg = self.args[0].as_leading_term(x, logx=logx, cdir=cdir)
        arg0 = arg.subs(x, 0)

        if arg0 is S.NaN:
            arg0 = arg.limit(x, 0, dir='-' if cdir.is_negative else '+')
        if arg0.is_zero:
            return S.One
        elif arg0.is_finite:
            return self.func(arg0)
        else:
            return self

    def _eval_is_real(self):
        arg = self.args[0]

        # `cosh(x)` is real for real OR purely imaginary `x`
        if arg.is_real or arg.is_imaginary:
            return True

        # cosh(a+ib) = cos(b)*cosh(a) + i*sin(b)*sinh(a)
        # the imaginary part can be an expression like n*pi
        # if not, check if the imaginary part is a number
        re, im = arg.as_real_imag()
        return (im%pi).is_zero

    def _eval_is_positive(self):
        # cosh(x+I*y) = cos(y)*cosh(x) + I*sin(y)*sinh(x)
        # cosh(z) is positive iff it is real and the real part is positive.
        # So we need sin(y)*sinh(x) = 0 which gives x=0 or y=n*pi
        # Case 1 (y=n*pi): cosh(z) = (-1)**n * cosh(x) -> positive for n even
        # Case 2 (x=0): cosh(z) = cos(y) -> positive when cos(y) is positive
        z = self.args[0]

        x, y = z.as_real_imag()
        ymod = y % (2*pi)

        yzero = ymod.is_zero
        # shortcut if ymod is zero
        if yzero:
            return True

        xzero = x.is_zero
        # shortcut x is not zero
        if xzero is False:
            return yzero

        return fuzzy_or([
                # Case 1:
                yzero,
                # Case 2:
                fuzzy_and([
                    xzero,
                    fuzzy_or([ymod < pi/2, ymod > 3*pi/2])
                ])
            ])


    def _eval_is_nonnegative(self):
        z = self.args[0]

        x, y = z.as_real_imag()
        ymod = y % (2*pi)

        yzero = ymod.is_zero
        # shortcut if ymod is zero
        if yzero:
            return True

        xzero = x.is_zero
        # shortcut x is not zero
        if xzero is False:
            return yzero

        return fuzzy_or([
                # Case 1:
                yzero,
                # Case 2:
                fuzzy_and([
                    xzero,
                    fuzzy_or([ymod <= pi/2, ymod >= 3*pi/2])
                ])
            ])

    def _eval_is_finite(self):
        arg = self.args[0]
        return arg.is_finite

    def _eval_is_zero(self):
        rest, ipi_mult = _peeloff_ipi(self.args[0])
        if ipi_mult and rest.is_zero:
            return (ipi_mult - S.Half).is_integer


class tanh(HyperbolicFunction):
    r"""
    ``tanh(x)`` is the hyperbolic tangent of ``x``.

    The hyperbolic tangent function is $\frac{\sinh(x)}{\cosh(x)}$.

    Examples
    ========

    >>> from sympy import tanh
    >>> from sympy.abc import x
    >>> tanh(x)
    tanh(x)

    See Also
    ========

    sinh, cosh, atanh
    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return S.One - tanh(self.args[0])**2
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        """
        Returns the inverse of this function.
        """
        return atanh

    @classmethod
    def eval(cls, arg):
        from sympy.functions.elementary.trigonometric import tan
        arg = sympify(arg)

        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.One
            elif arg is S.NegativeInfinity:
                return S.NegativeOne
            elif arg.is_zero:
                return S.Zero
            elif arg.is_negative:
                return -cls(-arg)
        else:
            if arg is S.ComplexInfinity:
                return S.NaN

            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                if i_coeff.could_extract_minus_sign():
                    return -S.ImaginaryUnit * tan(-i_coeff)
                return S.ImaginaryUnit * tan(i_coeff)
            else:
                if arg.could_extract_minus_sign():
                    return -cls(-arg)

            if arg.is_Add:
                x, m = _peeloff_ipi(arg)
                if m:
                    tanhm = tanh(m*S.Pi*S.ImaginaryUnit)
                    if tanhm is S.ComplexInfinity:
                        return coth(x)
                    else: # tanhm == 0
                        return tanh(x)

            if arg.is_zero:
                return S.Zero

            if arg.func == asinh:
                x = arg.args[0]
                return x/sqrt(1 + x**2)

            if arg.func == acosh:
                x = arg.args[0]
                return sqrt(x - 1) * sqrt(x + 1) / x

            if arg.func == atanh:
                return arg.args[0]

            if arg.func == acoth:
                return 1/arg.args[0]

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        from sympy.functions.combinatorial.numbers import bernoulli
        if n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)

            a = 2**(n + 1)

            B = bernoulli(n + 1)
            F = factorial(n + 1)

            return a*(a - 1) * B/F * x**n

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def as_real_imag(self, deep=True, **hints):
        from sympy.functions.elementary.trigonometric import (cos, sin)
        if self.args[0].is_extended_real:
            if deep:
                hints['complex'] = False
                return (self.expand(deep, **hints), S.Zero)
            else:
                return (self, S.Zero)
        if deep:
            re, im = self.args[0].expand(deep, **hints).as_real_imag()
        else:
            re, im = self.args[0].as_real_imag()
        denom = sinh(re)**2 + cos(im)**2
        return (sinh(re)*cosh(re)/denom, sin(im)*cos(im)/denom)

    def _eval_expand_trig(self, **hints):
        arg = self.args[0]
        if arg.is_Add:
            from sympy.polys.specialpolys import symmetric_poly
            n = len(arg.args)
            TX = [tanh(x, evaluate=False)._eval_expand_trig()
                for x in arg.args]
            p = [0, 0]  # [den, num]
            for i in range(n + 1):
                p[i % 2] += symmetric_poly(i, TX)
            return p[1]/p[0]
        elif arg.is_Mul:
            from sympy.functions.combinatorial.numbers import nC
            coeff, terms = arg.as_coeff_Mul()
            if coeff.is_Integer and coeff > 1:
                n = []
                d = []
                T = tanh(terms)
                for k in range(1, coeff + 1, 2):
                    n.append(nC(range(coeff), k)*T**k)
                for k in range(0, coeff + 1, 2):
                    d.append(nC(range(coeff), k)*T**k)
                return Add(*n)/Add(*d)
        return tanh(arg)

    def _eval_rewrite_as_tractable(self, arg, limitvar=None, **kwargs):
        neg_exp, pos_exp = exp(-arg), exp(arg)
        return (pos_exp - neg_exp)/(pos_exp + neg_exp)

    def _eval_rewrite_as_exp(self, arg, **kwargs):
        neg_exp, pos_exp = exp(-arg), exp(arg)
        return (pos_exp - neg_exp)/(pos_exp + neg_exp)

    def _eval_rewrite_as_sinh(self, arg, **kwargs):
        return S.ImaginaryUnit*sinh(arg)/sinh(S.Pi*S.ImaginaryUnit/2 - arg)

    def _eval_rewrite_as_cosh(self, arg, **kwargs):
        return S.ImaginaryUnit*cosh(S.Pi*S.ImaginaryUnit/2 - arg)/cosh(arg)

    def _eval_rewrite_as_coth(self, arg, **kwargs):
        return 1/coth(arg)

    def _eval_as_leading_term(self, x, logx=None, cdir=0):
        from sympy.series.order import Order
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and Order(1, x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_is_real(self):
        arg = self.args[0]
        if arg.is_real:
            return True

        re, im = arg.as_real_imag()

        # if denom = 0, tanh(arg) = zoo
        if re == 0 and im % pi == pi/2:
            return None

        # check if im is of the form n*pi/2 to make sin(2*im) = 0
        # if not, im could be a number, return False in that case
        return (im % (pi/2)).is_zero

    def _eval_is_extended_real(self):
        if self.args[0].is_extended_real:
            return True

    def _eval_is_positive(self):
        if self.args[0].is_extended_real:
            return self.args[0].is_positive

    def _eval_is_negative(self):
        if self.args[0].is_extended_real:
            return self.args[0].is_negative

    def _eval_is_finite(self):
        from sympy.functions.elementary.trigonometric import cos
        arg = self.args[0]

        re, im = arg.as_real_imag()
        denom = cos(im)**2 + sinh(re)**2
        if denom == 0:
            return False
        elif denom.is_number:
            return True
        if arg.is_extended_real:
            return True

    def _eval_is_zero(self):
        arg = self.args[0]
        if arg.is_zero:
            return True


class coth(HyperbolicFunction):
    r"""
    ``coth(x)`` is the hyperbolic cotangent of ``x``.

    The hyperbolic cotangent function is $\frac{\cosh(x)}{\sinh(x)}$.

    Examples
    ========

    >>> from sympy import coth
    >>> from sympy.abc import x
    >>> coth(x)
    coth(x)

    See Also
    ========

    sinh, cosh, acoth
    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -1/sinh(self.args[0])**2
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        """
        Returns the inverse of this function.
        """
        return acoth

    @classmethod
    def eval(cls, arg):
        from sympy.functions.elementary.trigonometric import cot
        arg = sympify(arg)

        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.One
            elif arg is S.NegativeInfinity:
                return S.NegativeOne
            elif arg.is_zero:
                return S.ComplexInfinity
            elif arg.is_negative:
                return -cls(-arg)
        else:
            if arg is S.ComplexInfinity:
                return S.NaN

            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                if i_coeff.could_extract_minus_sign():
                    return S.ImaginaryUnit * cot(-i_coeff)
                return -S.ImaginaryUnit * cot(i_coeff)
            else:
                if arg.could_extract_minus_sign():
                    return -cls(-arg)

            if arg.is_Add:
                x, m = _peeloff_ipi(arg)
                if m:
                    cothm = coth(m*S.Pi*S.ImaginaryUnit)
                    if cothm is S.ComplexInfinity:
                        return coth(x)
                    else: # cothm == 0
                        return tanh(x)

            if arg.is_zero:
                return S.ComplexInfinity

            if arg.func == asinh:
                x = arg.args[0]
                return sqrt(1 + x**2)/x

            if arg.func == acosh:
                x = arg.args[0]
                return x/(sqrt(x - 1) * sqrt(x + 1))

            if arg.func == atanh:
                return 1/arg.args[0]

            if arg.func == acoth:
                return arg.args[0]

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        from sympy.functions.combinatorial.numbers import bernoulli
        if n == 0:
            return 1 / sympify(x)
        elif n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)

            B = bernoulli(n + 1)
            F = factorial(n + 1)

            return 2**(n + 1) * B/F * x**n

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def as_real_imag(self, deep=True, **hints):
        from sympy.functions.elementary.trigonometric import (cos, sin)
        if self.args[0].is_extended_real:
            if deep:
                hints['complex'] = False
                return (self.expand(deep, **hints), S.Zero)
            else:
                return (self, S.Zero)
        if deep:
            re, im = self.args[0].expand(deep, **hints).as_real_imag()
        else:
            re, im = self.args[0].as_real_imag()
        denom = sinh(re)**2 + sin(im)**2
        return (sinh(re)*cosh(re)/denom, -sin(im)*cos(im)/denom)

    def _eval_rewrite_as_tractable(self, arg, limitvar=None, **kwargs):
        neg_exp, pos_exp = exp(-arg), exp(arg)
        return (pos_exp + neg_exp)/(pos_exp - neg_exp)

    def _eval_rewrite_as_exp(self, arg, **kwargs):
        neg_exp, pos_exp = exp(-arg), exp(arg)
        return (pos_exp + neg_exp)/(pos_exp - neg_exp)

    def _eval_rewrite_as_sinh(self, arg, **kwargs):
        return -S.ImaginaryUnit*sinh(S.Pi*S.ImaginaryUnit/2 - arg)/sinh(arg)

    def _eval_rewrite_as_cosh(self, arg, **kwargs):
        return -S.ImaginaryUnit*cosh(arg)/cosh(S.Pi*S.ImaginaryUnit/2 - arg)

    def _eval_rewrite_as_tanh(self, arg, **kwargs):
        return 1/tanh(arg)

    def _eval_is_positive(self):
        if self.args[0].is_extended_real:
            return self.args[0].is_positive

    def _eval_is_negative(self):
        if self.args[0].is_extended_real:
            return self.args[0].is_negative

    def _eval_as_leading_term(self, x, logx=None, cdir=0):
        from sympy.series.order import Order
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and Order(1, x).contains(arg):
            return 1/arg
        else:
            return self.func(arg)

    def _eval_expand_trig(self, **hints):
        arg = self.args[0]
        if arg.is_Add:
            from sympy.polys.specialpolys import symmetric_poly
            CX = [coth(x, evaluate=False)._eval_expand_trig() for x in arg.args]
            p = [[], []]
            n = len(arg.args)
            for i in range(n, -1, -1):
                p[(n - i) % 2].append(symmetric_poly(i, CX))
            return Add(*p[0])/Add(*p[1])
        elif arg.is_Mul:
            from sympy.functions.combinatorial.factorials import binomial
            coeff, x = arg.as_coeff_Mul(rational=True)
            if coeff.is_Integer and coeff > 1:
                c = coth(x, evaluate=False)
                p = [[], []]
                for i in range(coeff, -1, -1):
                    p[(coeff - i) % 2].append(binomial(coeff, i)*c**i)
                return Add(*p[0])/Add(*p[1])
        return coth(arg)


class ReciprocalHyperbolicFunction(HyperbolicFunction):
    """Base class for reciprocal functions of hyperbolic functions. """

    #To be defined in class
    _reciprocal_of = None
    _is_even = None  # type: FuzzyBool
    _is_odd = None  # type: FuzzyBool

    @classmethod
    def eval(cls, arg):
        if arg.could_extract_minus_sign():
            if cls._is_even:
                return cls(-arg)
            if cls._is_odd:
                return -cls(-arg)

        t = cls._reciprocal_of.eval(arg)
        if hasattr(arg, 'inverse') and arg.inverse() == cls:
            return arg.args[0]
        return 1/t if t is not None else t

    def _call_reciprocal(self, method_name, *args, **kwargs):
        # Calls method_name on _reciprocal_of
        o = self._reciprocal_of(self.args[0])
        return getattr(o, method_name)(*args, **kwargs)

    def _calculate_reciprocal(self, method_name, *args, **kwargs):
        # If calling method_name on _reciprocal_of returns a value != None
        # then return the reciprocal of that value
        t = self._call_reciprocal(method_name, *args, **kwargs)
        return 1/t if t is not None else t

    def _rewrite_reciprocal(self, method_name, arg):
        # Special handling for rewrite functions. If reciprocal rewrite returns
        # unmodified expression, then return None
        t = self._call_reciprocal(method_name, arg)
        if t is not None and t != self._reciprocal_of(arg):
            return 1/t

    def _eval_rewrite_as_exp(self, arg, **kwargs):
        return self._rewrite_reciprocal("_eval_rewrite_as_exp", arg)

    def _eval_rewrite_as_tractable(self, arg, limitvar=None, **kwargs):
        return self._rewrite_reciprocal("_eval_rewrite_as_tractable", arg)

    def _eval_rewrite_as_tanh(self, arg, **kwargs):
        return self._rewrite_reciprocal("_eval_rewrite_as_tanh", arg)

    def _eval_rewrite_as_coth(self, arg, **kwargs):
        return self._rewrite_reciprocal("_eval_rewrite_as_coth", arg)

    def as_real_imag(self, deep = True, **hints):
        return (1 / self._reciprocal_of(self.args[0])).as_real_imag(deep, **hints)

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def _eval_expand_complex(self, deep=True, **hints):
        re_part, im_part = self.as_real_imag(deep=True, **hints)
        return re_part + S.ImaginaryUnit*im_part

    def _eval_expand_trig(self, **hints):
        return self._calculate_reciprocal("_eval_expand_trig", **hints)

    def _eval_as_leading_term(self, x, logx=None, cdir=0):
        return (1/self._reciprocal_of(self.args[0]))._eval_as_leading_term(x)

    def _eval_is_extended_real(self):
        return self._reciprocal_of(self.args[0]).is_extended_real

    def _eval_is_finite(self):
        return (1/self._reciprocal_of(self.args[0])).is_finite


class csch(ReciprocalHyperbolicFunction):
    r"""
    ``csch(x)`` is the hyperbolic cosecant of ``x``.

    The hyperbolic cosecant function is $\frac{2}{e^x - e^{-x}}$

    Examples
    ========

    >>> from sympy import csch
    >>> from sympy.abc import x
    >>> csch(x)
    csch(x)

    See Also
    ========

    sinh, cosh, tanh, sech, asinh, acosh
    """

    _reciprocal_of = sinh
    _is_odd = True

    def fdiff(self, argindex=1):
        """
        Returns the first derivative of this function
        """
        if argindex == 1:
            return -coth(self.args[0]) * csch(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        """
        Returns the next term in the Taylor series expansion
        """
        from sympy.functions.combinatorial.numbers import bernoulli
        if n == 0:
            return 1/sympify(x)
        elif n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)

            B = bernoulli(n + 1)
            F = factorial(n + 1)

            return 2 * (1 - 2**n) * B/F * x**n

    def _eval_rewrite_as_cosh(self, arg, **kwargs):
        return S.ImaginaryUnit / cosh(arg + S.ImaginaryUnit * S.Pi / 2)

    def _eval_is_positive(self):
        if self.args[0].is_extended_real:
            return self.args[0].is_positive

    def _eval_is_negative(self):
        if self.args[0].is_extended_real:
            return self.args[0].is_negative


class sech(ReciprocalHyperbolicFunction):
    r"""
    ``sech(x)`` is the hyperbolic secant of ``x``.

    The hyperbolic secant function is $\frac{2}{e^x + e^{-x}}$

    Examples
    ========

    >>> from sympy import sech
    >>> from sympy.abc import x
    >>> sech(x)
    sech(x)

    See Also
    ========

    sinh, cosh, tanh, coth, csch, asinh, acosh
    """

    _reciprocal_of = cosh
    _is_even = True

    def fdiff(self, argindex=1):
        if argindex == 1:
            return - tanh(self.args[0])*sech(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        from sympy.functions.combinatorial.numbers import euler
        if n < 0 or n % 2 == 1:
            return S.Zero
        else:
            x = sympify(x)
            return euler(n) / factorial(n) * x**(n)

    def _eval_rewrite_as_sinh(self, arg, **kwargs):
        return S.ImaginaryUnit / sinh(arg + S.ImaginaryUnit * S.Pi /2)

    def _eval_is_positive(self):
        if self.args[0].is_extended_real:
            return True


###############################################################################
############################# HYPERBOLIC INVERSES #############################
###############################################################################

class InverseHyperbolicFunction(Function):
    """Base class for inverse hyperbolic functions."""

    pass


class asinh(InverseHyperbolicFunction):
    """
    ``asinh(x)`` is the inverse hyperbolic sine of ``x``.

    The inverse hyperbolic sine function.

    Examples
    ========

    >>> from sympy import asinh
    >>> from sympy.abc import x
    >>> asinh(x).diff(x)
    1/sqrt(x**2 + 1)
    >>> asinh(1)
    log(1 + sqrt(2))

    See Also
    ========

    acosh, atanh, sinh
    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/sqrt(self.args[0]**2 + 1)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        from sympy.functions.elementary.trigonometric import asin
        arg = sympify(arg)

        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.Infinity
            elif arg is S.NegativeInfinity:
                return S.NegativeInfinity
            elif arg.is_zero:
                return S.Zero
            elif arg is S.One:
                return log(sqrt(2) + 1)
            elif arg is S.NegativeOne:
                return log(sqrt(2) - 1)
            elif arg.is_negative:
                return -cls(-arg)
        else:
            if arg is S.ComplexInfinity:
                return S.ComplexInfinity

            if arg.is_zero:
                return S.Zero

            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return S.ImaginaryUnit * asin(i_coeff)
            else:
                if arg.could_extract_minus_sign():
                    return -cls(-arg)

        if isinstance(arg, sinh) and arg.args[0].is_number:
            z = arg.args[0]
            if z.is_real:
                return z
            r, i = match_real_imag(z)
            if r is not None and i is not None:
                f = floor((i + pi/2)/pi)
                m = z - I*pi*f
                even = f.is_even
                if even is True:
                    return m
                elif even is False:
                    return -m

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)
            if len(previous_terms) >= 2 and n > 2:
                p = previous_terms[-2]
                return -p * (n - 2)**2/(n*(n - 1)) * x**2
            else:
                k = (n - 1) // 2
                R = RisingFactorial(S.Half, k)
                F = factorial(k)
                return S.NegativeOne**k * R / F * x**n / n

    def _eval_as_leading_term(self, x, logx=None, cdir=0):
        from sympy.series.order import Order
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and Order(1, x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_rewrite_as_log(self, x, **kwargs):
        return log(x + sqrt(x**2 + 1))

    def inverse(self, argindex=1):
        """
        Returns the inverse of this function.
        """
        return sinh

    def _eval_is_zero(self):
        return self.args[0].is_zero


class acosh(InverseHyperbolicFunction):
    """
    ``acosh(x)`` is the inverse hyperbolic cosine of ``x``.

    The inverse hyperbolic cosine function.

    Examples
    ========

    >>> from sympy import acosh
    >>> from sympy.abc import x
    >>> acosh(x).diff(x)
    1/sqrt(x**2 - 1)
    >>> acosh(1)
    0

    See Also
    ========

    asinh, atanh, cosh
    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/sqrt(self.args[0]**2 - 1)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        arg = sympify(arg)

        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.Infinity
            elif arg is S.NegativeInfinity:
                return S.Infinity
            elif arg.is_zero:
                return S.Pi*S.ImaginaryUnit / 2
            elif arg is S.One:
                return S.Zero
            elif arg is S.NegativeOne:
                return S.Pi*S.ImaginaryUnit

        if arg.is_number:
            cst_table = {
                S.ImaginaryUnit: log(S.ImaginaryUnit*(1 + sqrt(2))),
                -S.ImaginaryUnit: log(-S.ImaginaryUnit*(1 + sqrt(2))),
                S.Half: S.Pi/3,
                Rational(-1, 2): S.Pi*Rational(2, 3),
                sqrt(2)/2: S.Pi/4,
                -sqrt(2)/2: S.Pi*Rational(3, 4),
                1/sqrt(2): S.Pi/4,
                -1/sqrt(2): S.Pi*Rational(3, 4),
                sqrt(3)/2: S.Pi/6,
                -sqrt(3)/2: S.Pi*Rational(5, 6),
                (sqrt(3) - 1)/sqrt(2**3): S.Pi*Rational(5, 12),
                -(sqrt(3) - 1)/sqrt(2**3): S.Pi*Rational(7, 12),
                sqrt(2 + sqrt(2))/2: S.Pi/8,
                -sqrt(2 + sqrt(2))/2: S.Pi*Rational(7, 8),
                sqrt(2 - sqrt(2))/2: S.Pi*Rational(3, 8),
                -sqrt(2 - sqrt(2))/2: S.Pi*Rational(5, 8),
                (1 + sqrt(3))/(2*sqrt(2)): S.Pi/12,
                -(1 + sqrt(3))/(2*sqrt(2)): S.Pi*Rational(11, 12),
                (sqrt(5) + 1)/4: S.Pi/5,
                -(sqrt(5) + 1)/4: S.Pi*Rational(4, 5)
            }

            if arg in cst_table:
                if arg.is_extended_real:
                    return cst_table[arg]*S.ImaginaryUnit
                return cst_table[arg]

        if arg is S.ComplexInfinity:
            return S.ComplexInfinity
        if arg == S.ImaginaryUnit*S.Infinity:
            return S.Infinity + S.ImaginaryUnit*S.Pi/2
        if arg == -S.ImaginaryUnit*S.Infinity:
            return S.Infinity - S.ImaginaryUnit*S.Pi/2

        if arg.is_zero:
            return S.Pi*S.ImaginaryUnit*S.Half

        if isinstance(arg, cosh) and arg.args[0].is_number:
            z = arg.args[0]
            if z.is_real:
                from sympy.functions.elementary.complexes import Abs
                return Abs(z)
            r, i = match_real_imag(z)
            if r is not None and i is not None:
                f = floor(i/pi)
                m = z - I*pi*f
                even = f.is_even
                if even is True:
                    if r.is_nonnegative:
                        return m
                    elif r.is_negative:
                        return -m
                elif even is False:
                    m -= I*pi
                    if r.is_nonpositive:
                        return -m
                    elif r.is_positive:
                        return m

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n == 0:
            return S.Pi*S.ImaginaryUnit / 2
        elif n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)
            if len(previous_terms) >= 2 and n > 2:
                p = previous_terms[-2]
                return p * (n - 2)**2/(n*(n - 1)) * x**2
            else:
                k = (n - 1) // 2
                R = RisingFactorial(S.Half, k)
                F = factorial(k)
                return -R / F * S.ImaginaryUnit * x**n / n

    def _eval_as_leading_term(self, x, logx=None, cdir=0):
        from sympy.series.order import Order
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and Order(1, x).contains(arg):
            return S.ImaginaryUnit*S.Pi/2
        else:
            return self.func(arg)

    def _eval_rewrite_as_log(self, x, **kwargs):
        return log(x + sqrt(x + 1) * sqrt(x - 1))

    def inverse(self, argindex=1):
        """
        Returns the inverse of this function.
        """
        return cosh

    def _eval_is_zero(self):
        if (self.args[0] - 1).is_zero:
            return True


class atanh(InverseHyperbolicFunction):
    """
    ``atanh(x)`` is the inverse hyperbolic tangent of ``x``.

    The inverse hyperbolic tangent function.

    Examples
    ========

    >>> from sympy import atanh
    >>> from sympy.abc import x
    >>> atanh(x).diff(x)
    1/(1 - x**2)

    See Also
    ========

    asinh, acosh, tanh
    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/(1 - self.args[0]**2)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        from sympy.functions.elementary.trigonometric import atan
        arg = sympify(arg)

        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg.is_zero:
                return S.Zero
            elif arg is S.One:
                return S.Infinity
            elif arg is S.NegativeOne:
                return S.NegativeInfinity
            elif arg is S.Infinity:
                return -S.ImaginaryUnit * atan(arg)
            elif arg is S.NegativeInfinity:
                return S.ImaginaryUnit * atan(-arg)
            elif arg.is_negative:
                return -cls(-arg)
        else:
            if arg is S.ComplexInfinity:
                from sympy.calculus.accumulationbounds import AccumBounds
                return S.ImaginaryUnit*AccumBounds(-S.Pi/2, S.Pi/2)

            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return S.ImaginaryUnit * atan(i_coeff)
            else:
                if arg.could_extract_minus_sign():
                    return -cls(-arg)

        if arg.is_zero:
            return S.Zero

        if isinstance(arg, tanh) and arg.args[0].is_number:
            z = arg.args[0]
            if z.is_real:
                return z
            r, i = match_real_imag(z)
            if r is not None and i is not None:
                f = floor(2*i/pi)
                even = f.is_even
                m = z - I*f*pi/2
                if even is True:
                    return m
                elif even is False:
                    return m - I*pi/2

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)
            return x**n / n

    def _eval_as_leading_term(self, x, logx=None, cdir=0):
        from sympy.series.order import Order
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and Order(1, x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_rewrite_as_log(self, x, **kwargs):
        return (log(1 + x) - log(1 - x)) / 2

    def _eval_is_zero(self):
        if self.args[0].is_zero:
            return True

    def _eval_is_imaginary(self):
        return self.args[0].is_imaginary

    def inverse(self, argindex=1):
        """
        Returns the inverse of this function.
        """
        return tanh


class acoth(InverseHyperbolicFunction):
    """
    ``acoth(x)`` is the inverse hyperbolic cotangent of ``x``.

    The inverse hyperbolic cotangent function.

    Examples
    ========

    >>> from sympy import acoth
    >>> from sympy.abc import x
    >>> acoth(x).diff(x)
    1/(1 - x**2)

    See Also
    ========

    asinh, acosh, coth
    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/(1 - self.args[0]**2)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        from sympy.functions.elementary.trigonometric import acot
        arg = sympify(arg)

        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.Zero
            elif arg is S.NegativeInfinity:
                return S.Zero
            elif arg.is_zero:
                return S.Pi*S.ImaginaryUnit / 2
            elif arg is S.One:
                return S.Infinity
            elif arg is S.NegativeOne:
                return S.NegativeInfinity
            elif arg.is_negative:
                return -cls(-arg)
        else:
            if arg is S.ComplexInfinity:
                return S.Zero

            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return -S.ImaginaryUnit * acot(i_coeff)
            else:
                if arg.could_extract_minus_sign():
                    return -cls(-arg)

        if arg.is_zero:
            return S.Pi*S.ImaginaryUnit*S.Half

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n == 0:
            return S.Pi*S.ImaginaryUnit / 2
        elif n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)
            return x**n / n

    def _eval_as_leading_term(self, x, logx=None, cdir=0):
        from sympy.series.order import Order
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and Order(1, x).contains(arg):
            return S.ImaginaryUnit*S.Pi/2
        else:
            return self.func(arg)

    def _eval_rewrite_as_log(self, x, **kwargs):
        return (log(1 + 1/x) - log(1 - 1/x)) / 2

    def inverse(self, argindex=1):
        """
        Returns the inverse of this function.
        """
        return coth


class asech(InverseHyperbolicFunction):
    """
    ``asech(x)`` is the inverse hyperbolic secant of ``x``.

    The inverse hyperbolic secant function.

    Examples
    ========

    >>> from sympy import asech, sqrt, S
    >>> from sympy.abc import x
    >>> asech(x).diff(x)
    -1/(x*sqrt(1 - x**2))
    >>> asech(1).diff(x)
    0
    >>> asech(1)
    0
    >>> asech(S(2))
    I*pi/3
    >>> asech(-sqrt(2))
    3*I*pi/4
    >>> asech((sqrt(6) - sqrt(2)))
    I*pi/12

    See Also
    ========

    asinh, atanh, cosh, acoth

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Hyperbolic_function
    .. [2] http://dlmf.nist.gov/4.37
    .. [3] http://functions.wolfram.com/ElementaryFunctions/ArcSech/

    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            z = self.args[0]
            return -1/(z*sqrt(1 - z**2))
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        arg = sympify(arg)

        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.Pi*S.ImaginaryUnit / 2
            elif arg is S.NegativeInfinity:
                return S.Pi*S.ImaginaryUnit / 2
            elif arg.is_zero:
                return S.Infinity
            elif arg is S.One:
                return S.Zero
            elif arg is S.NegativeOne:
                return S.Pi*S.ImaginaryUnit

        if arg.is_number:
            cst_table = {
                S.ImaginaryUnit: - (S.Pi*S.ImaginaryUnit / 2) + log(1 + sqrt(2)),
                -S.ImaginaryUnit: (S.Pi*S.ImaginaryUnit / 2) + log(1 + sqrt(2)),
                (sqrt(6) - sqrt(2)): S.Pi / 12,
                (sqrt(2) - sqrt(6)): 11*S.Pi / 12,
                sqrt(2 - 2/sqrt(5)): S.Pi / 10,
                -sqrt(2 - 2/sqrt(5)): 9*S.Pi / 10,
                2 / sqrt(2 + sqrt(2)): S.Pi / 8,
                -2 / sqrt(2 + sqrt(2)): 7*S.Pi / 8,
                2 / sqrt(3): S.Pi / 6,
                -2 / sqrt(3): 5*S.Pi / 6,
                (sqrt(5) - 1): S.Pi / 5,
                (1 - sqrt(5)): 4*S.Pi / 5,
                sqrt(2): S.Pi / 4,
                -sqrt(2): 3*S.Pi / 4,
                sqrt(2 + 2/sqrt(5)): 3*S.Pi / 10,
                -sqrt(2 + 2/sqrt(5)): 7*S.Pi / 10,
                S(2): S.Pi / 3,
                -S(2): 2*S.Pi / 3,
                sqrt(2*(2 + sqrt(2))): 3*S.Pi / 8,
                -sqrt(2*(2 + sqrt(2))): 5*S.Pi / 8,
                (1 + sqrt(5)): 2*S.Pi / 5,
                (-1 - sqrt(5)): 3*S.Pi / 5,
                (sqrt(6) + sqrt(2)): 5*S.Pi / 12,
                (-sqrt(6) - sqrt(2)): 7*S.Pi / 12,
            }

            if arg in cst_table:
                if arg.is_extended_real:
                    return cst_table[arg]*S.ImaginaryUnit
                return cst_table[arg]

        if arg is S.ComplexInfinity:
            from sympy.calculus.accumulationbounds import AccumBounds
            return S.ImaginaryUnit*AccumBounds(-S.Pi/2, S.Pi/2)

        if arg.is_zero:
            return S.Infinity

    @staticmethod
    @cacheit
    def expansion_term(n, x, *previous_terms):
        if n == 0:
            return log(2 / x)
        elif n < 0 or n % 2 == 1:
            return S.Zero
        else:
            x = sympify(x)
            if len(previous_terms) > 2 and n > 2:
                p = previous_terms[-2]
                return p * (n - 1)**2 // (n // 2)**2 * x**2 / 4
            else:
                k = n // 2
                R = RisingFactorial(S.Half, k) *  n
                F = factorial(k) * n // 2 * n // 2
                return -1 * R / F * x**n / 4

    def inverse(self, argindex=1):
        """
        Returns the inverse of this function.
        """
        return sech

    def _eval_rewrite_as_log(self, arg, **kwargs):
        return log(1/arg + sqrt(1/arg - 1) * sqrt(1/arg + 1))


class acsch(InverseHyperbolicFunction):
    """
    ``acsch(x)`` is the inverse hyperbolic cosecant of ``x``.

    The inverse hyperbolic cosecant function.

    Examples
    ========

    >>> from sympy import acsch, sqrt, S
    >>> from sympy.abc import x
    >>> acsch(x).diff(x)
    -1/(x**2*sqrt(1 + x**(-2)))
    >>> acsch(1).diff(x)
    0
    >>> acsch(1)
    log(1 + sqrt(2))
    >>> acsch(S.ImaginaryUnit)
    -I*pi/2
    >>> acsch(-2*S.ImaginaryUnit)
    I*pi/6
    >>> acsch(S.ImaginaryUnit*(sqrt(6) - sqrt(2)))
    -5*I*pi/12

    See Also
    ========

    asinh

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Hyperbolic_function
    .. [2] http://dlmf.nist.gov/4.37
    .. [3] http://functions.wolfram.com/ElementaryFunctions/ArcCsch/

    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            z = self.args[0]
            return -1/(z**2*sqrt(1 + 1/z**2))
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        arg = sympify(arg)

        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.Zero
            elif arg is S.NegativeInfinity:
                return S.Zero
            elif arg.is_zero:
                return S.ComplexInfinity
            elif arg is S.One:
                return log(1 + sqrt(2))
            elif arg is S.NegativeOne:
                return - log(1 + sqrt(2))

        if arg.is_number:
            cst_table = {
                S.ImaginaryUnit: -S.Pi / 2,
                S.ImaginaryUnit*(sqrt(2) + sqrt(6)): -S.Pi / 12,
                S.ImaginaryUnit*(1 + sqrt(5)): -S.Pi / 10,
                S.ImaginaryUnit*2 / sqrt(2 - sqrt(2)): -S.Pi / 8,
                S.ImaginaryUnit*2: -S.Pi / 6,
                S.ImaginaryUnit*sqrt(2 + 2/sqrt(5)): -S.Pi / 5,
                S.ImaginaryUnit*sqrt(2): -S.Pi / 4,
                S.ImaginaryUnit*(sqrt(5)-1): -3*S.Pi / 10,
                S.ImaginaryUnit*2 / sqrt(3): -S.Pi / 3,
                S.ImaginaryUnit*2 / sqrt(2 + sqrt(2)): -3*S.Pi / 8,
                S.ImaginaryUnit*sqrt(2 - 2/sqrt(5)): -2*S.Pi / 5,
                S.ImaginaryUnit*(sqrt(6) - sqrt(2)): -5*S.Pi / 12,
                S(2): -S.ImaginaryUnit*log((1+sqrt(5))/2),
            }

            if arg in cst_table:
                return cst_table[arg]*S.ImaginaryUnit

        if arg is S.ComplexInfinity:
            return S.Zero

        if arg.is_zero:
            return S.ComplexInfinity

        if arg.could_extract_minus_sign():
            return -cls(-arg)

    def inverse(self, argindex=1):
        """
        Returns the inverse of this function.
        """
        return csch

    def _eval_rewrite_as_log(self, arg, **kwargs):
        return log(1/arg + sqrt(1/arg**2 + 1))

    def _eval_is_zero(self):
        return self.args[0].is_infinite

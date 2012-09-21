from sympy.core.add import Add
from sympy.core.numbers import Rational
from sympy.core.basic import C, sympify, cacheit
from sympy.core.singleton import S
from sympy.core.function import Function, ArgumentIndexError
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.exponential import log
from sympy.functions.elementary.hyperbolic import HyperbolicFunction

###############################################################################
########################## TRIGONOMETRIC FUNCTIONS ############################
###############################################################################

class TrigonometricFunction(Function):
    """Base class for trigonometric functions. """

    unbranched = True

def _peeloff_pi(arg):
    """
    Split ARG into two parts, a "rest" and a multiple of pi/2.
    This assumes ARG to be an Add.
    The multiple of pi returned in the second position is always a Rational.

    Examples:
    >>> from sympy.functions.elementary.trigonometric import _peeloff_pi as peel
    >>> from sympy import pi
    >>> from sympy.abc import x, y
    >>> peel(x + pi/2)
    (x, pi/2)
    >>> peel(x + 2*pi/3 + pi*y)
    (x + pi*y + pi/6, pi/2)
    """
    for a in Add.make_args(arg):
        if a is S.Pi:
            K = S.One
            break
        elif a.is_Mul:
            K, p = a.as_two_terms()
            if p is S.Pi and K.is_Rational:
                break
    else:
        return arg, S.Zero

    m1 = (K % S.Half) * S.Pi
    m2 = K*S.Pi - m1
    return arg - m2, m2

def _pi_coeff(arg, cycles=1):
    """
    When arg is a Number times pi (e.g. 3*pi/2) then return the Number
    normalized to be in the range [0, 2], else None.

    When an even multiple of pi is encountered, if it is multiplying
    something with known parity then the multiple is returned as 0 otherwise
    as 2.

    Examples:
    >>> from sympy.functions.elementary.trigonometric import _pi_coeff as coeff
    >>> from sympy import pi
    >>> from sympy.abc import x, y
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

    """
    arg = sympify(arg)
    if arg is S.Pi:
        return S.One
    elif not arg:
        return S.Zero
    elif arg.is_Mul:
        cx = arg.coeff(S.Pi)
        if cx:
            c, x = cx.as_coeff_Mul() # pi is not included as coeff
            if c.is_Float:
                # recast exact binary fractions to Rationals
                f = abs(c) % 1
                if f != 0:
                    p = -round(log(f, 2).evalf())
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
                    if x.is_even is not None: # known parity
                        return S.Zero
                    return 2*x
                else:
                    return c2*x
            return cx

class sin(TrigonometricFunction):
    """
    The sine function.

    * sin(x) -> Returns the sine of x (measured in radians)

    Notes
    =====

    * sin(x) will evaluate automatically in the case x
      is a multiple of pi, pi/2, pi/3, pi/4 and pi/6.

    Examples
    ========

    >>> from sympy import sin, pi
    >>> from sympy.abc import x
    >>> sin(x**2).diff(x)
    2*x*cos(x**2)
    >>> sin(1).diff(x)
    0
    >>> sin(pi)
    0
    >>> sin(pi/2)
    1
    >>> sin(pi/6)
    1/2

    See Also
    ========

    cos, tan, asin

    References
    ==========

    U{Definitions in trigonometry<http://planetmath.org/encyclopedia/DefinitionsInTrigonometry.html>}

    """

    nargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return cos(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        """
        Returns the inverse of this function.
        """
        return asin

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Zero:
                return S.Zero
            elif arg is S.Infinity:
                return

        if arg.could_extract_minus_sign():
            return -cls(-arg)

        i_coeff = arg.as_coefficient(S.ImaginaryUnit)
        if i_coeff is not None:
            return S.ImaginaryUnit * C.sinh(i_coeff)

        pi_coeff = _pi_coeff(arg)
        if pi_coeff is not None:
            if pi_coeff.is_integer:
                return S.Zero

            if not pi_coeff.is_Rational:
                narg = pi_coeff*S.Pi
                if narg != arg:
                    return cls(narg)
                return None

            cst_table_some = {
                2 : S.One,
                3 : S.Half*sqrt(3),
                4 : S.Half*sqrt(2),
                6 : S.Half,
            }

            cst_table_more = {
                (1, 5) : sqrt((5 - sqrt(5)) / 8),
                (2, 5) : sqrt((5 + sqrt(5)) / 8)
            }

            p = pi_coeff.p
            q = pi_coeff.q

            Q, P = p // q, p % q

            try:
                result = cst_table_some[q]
            except KeyError:
                if abs(P) > q // 2:
                    P = q - P

                try:
                    result = cst_table_more[(P, q)]
                except KeyError:
                    if P != p:
                        result = cls(C.Rational(P, q)*S.Pi)
                    else:
                        newarg = pi_coeff*S.Pi
                        if newarg != arg:
                            return cls(newarg)
                        return None

            if Q % 2 == 1:
                return -result
            else:
                return result

        if arg.is_Add:
            x, m = _peeloff_pi(arg)
            if m:
                return sin(m)*cos(x)+cos(m)*sin(x)

        if arg.func is asin:
            return arg.args[0]

        if arg.func is atan:
            x = arg.args[0]
            return x / sqrt(1 + x**2)

        if arg.func is atan2:
            y, x = arg.args
            return y / sqrt(x**2 + y**2)

        if arg.func is acos:
            x = arg.args[0]
            return sqrt(1 - x**2)

        if arg.func is acot:
            x = arg.args[0]
            return 1 / (sqrt(1 + 1 / x**2) * x)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)

            if len(previous_terms) > 2:
                p = previous_terms[-2]
                return -p * x**2 / (n*(n-1))
            else:
                return (-1)**(n//2) * x**(n)/C.factorial(n)

    def _eval_rewrite_as_exp(self, arg):
        exp, I = C.exp, S.ImaginaryUnit
        if isinstance(arg, TrigonometricFunction) or isinstance(arg, HyperbolicFunction) :
            arg = arg.func(arg.args[0]).rewrite(exp)
        return (exp(arg*I) - exp(-arg*I)) / (2*I)

    def _eval_rewrite_as_Pow(self, arg):
        if arg.func is log:
            I = S.ImaginaryUnit
            x = arg.args[0]
            return I*x**-I / 2 - I*x**I /2

    def _eval_rewrite_as_cos(self, arg):
        return -cos(arg + S.Pi/2)

    def _eval_rewrite_as_tan(self, arg):
        tan_half = tan(S.Half*arg)
        return 2*tan_half/(1 + tan_half**2)

    def _eval_rewrite_as_sincos(self, arg):
        return sin(arg)*cos(arg)/cos(arg)

    def _eval_rewrite_as_cot(self, arg):
        cot_half = cot(S.Half*arg)
        return 2*cot_half/(1 + cot_half**2)

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def as_real_imag(self, deep=True, **hints):
        if self.args[0].is_real:
            if deep:
                hints['complex'] = False
                return (self.expand(deep, **hints), S.Zero)
            else:
                return (self, S.Zero)
        if deep:
            re, im = self.args[0].expand(deep, **hints).as_real_imag()
        else:
            re, im = self.args[0].as_real_imag()
        return (sin(re)*C.cosh(im), cos(re)*C.sinh(im))

    def _eval_expand_trig(self, **hints):
        from sympy import expand_mul
        arg = self.args[0]
        x = None
        if arg.is_Add: # TODO, implement more if deep stuff here
            # TODO: Do this more efficiently for more than two terms
            x, y = arg.as_two_terms()
            sx = sin(x, evaluate=False)._eval_expand_trig()
            sy = sin(y, evaluate=False)._eval_expand_trig()
            cx = cos(x, evaluate=False)._eval_expand_trig()
            cy = cos(y, evaluate=False)._eval_expand_trig()
            return sx*cy + sy*cx
        else:
            n, x = arg.as_coeff_Mul(rational=True)
            if n.is_Integer: # n will be positive because of .eval
                # canonicalization

                # See http://mathworld.wolfram.com/Multiple-AngleFormulas.html
                if n.is_odd:
                    return (-1)**((n - 1)/2)*C.chebyshevt(n, sin(x))
                else:
                    return expand_mul((-1)**(n/2 - 1)*cos(x)*C.chebyshevu(n -
                        1, sin(x)), deep=False)
        return sin(arg)

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and C.Order(1, x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_is_real(self):
        return self.args[0].is_real

    def _eval_is_bounded(self):
        arg = self.args[0]
        if arg.is_real:
            return True

    def _sage_(self):
        import sage.all as sage
        return sage.sin(self.args[0]._sage_())

class cos(TrigonometricFunction):
    """
    The cosine function.

    * cos(x) -> Returns the cosine of x (measured in radians)

    Notes
    =====

    * cos(x) will evaluate automatically in the case x
      is a multiple of pi, pi/2, pi/3, pi/4 and pi/6.

    Examples
    ========

    >>> from sympy import cos, pi
    >>> from sympy.abc import x
    >>> cos(x**2).diff(x)
    -2*x*sin(x**2)
    >>> cos(1).diff(x)
    0
    >>> cos(pi)
    -1
    >>> cos(pi/2)
    0
    >>> cos(2*pi/3)
    -1/2

    See Also
    ========

    sin, tan, acos

    References
    ==========

    U{Definitions in trigonometry<http://planetmath.org/encyclopedia/DefinitionsInTrigonometry.html>}

    """

    nargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -sin(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return acos

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Zero:
                return S.One
            elif arg is S.Infinity:
                return

        if arg.could_extract_minus_sign():
            return cls(-arg)

        i_coeff = arg.as_coefficient(S.ImaginaryUnit)
        if i_coeff is not None:
            return C.cosh(i_coeff)

        pi_coeff = _pi_coeff(arg)
        if pi_coeff is not None:
            if not pi_coeff.is_Rational:
                if pi_coeff.is_integer:
                    return (S.NegativeOne)**pi_coeff
                narg = pi_coeff*S.Pi
                if narg != arg:
                    return cls(narg)
                return None

            cst_table_some = {
                1 : S.One,
                2 : S.Zero,
                3 : S.Half,
                4 : S.Half*sqrt(2),
                6 : S.Half*sqrt(3),
            }

            cst_table_more = {
                (1, 5) : (sqrt(5) + 1)/4,
                (2, 5) : (sqrt(5) - 1)/4
            }

            p = pi_coeff.p
            q = pi_coeff.q

            Q, P = 2*p // q, p % q

            try:
                result = cst_table_some[q]
            except KeyError:
                if abs(P) > q // 2:
                    P = q - P

                try:
                    result = cst_table_more[(P, q)]
                except KeyError:
                    if P != p:
                        result = cls(C.Rational(P, q)*S.Pi)
                    else:
                        newarg = pi_coeff*S.Pi
                        if newarg != arg:
                            return cls(newarg)
                        return None

            if Q % 4 in (1, 2):
                return -result
            else:
                return result

        if arg.is_Add:
            x, m = _peeloff_pi(arg)
            if m:
                return cos(m)*cos(x)-sin(m)*sin(x)

        if arg.func is acos:
            return arg.args[0]

        if arg.func is atan:
            x = arg.args[0]
            return 1 / sqrt(1 + x**2)

        if arg.func is atan2:
            y, x = arg.args
            return x / sqrt(x**2 + y**2)

        if arg.func is asin:
            x = arg.args[0]
            return sqrt(1 - x ** 2)

        if arg.func is acot:
            x = arg.args[0]
            return 1 / sqrt(1 + 1 / x**2)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0 or n % 2 == 1:
            return S.Zero
        else:
            x = sympify(x)

            if len(previous_terms) > 2:
                p = previous_terms[-2]
                return -p * x**2 / (n*(n-1))
            else:
                return (-1)**(n//2)*x**(n)/C.factorial(n)

    def _eval_rewrite_as_exp(self, arg):
        exp, I = C.exp, S.ImaginaryUnit
        if isinstance(arg, TrigonometricFunction) or isinstance(arg, HyperbolicFunction) :
            arg = arg.func(arg.args[0]).rewrite(exp)
        return (exp(arg*I) + exp(-arg*I)) / 2

    def _eval_rewrite_as_Pow(self, arg):
        if arg.func is log:
            I = S.ImaginaryUnit
            x = arg.args[0]
            return x**I/2 + x**-I/2

    def _eval_rewrite_as_sin(self, arg):
        return sin(arg + S.Pi/2)

    def _eval_rewrite_as_tan(self, arg):
        tan_half = tan(S.Half*arg)**2
        return (1-tan_half)/(1+tan_half)

    def _eval_rewrite_as_sincos(self, arg):
        return sin(arg)*cos(arg)/sin(arg)

    def _eval_rewrite_as_cot(self, arg):
        cot_half = cot(S.Half*arg)**2
        return (cot_half-1)/(cot_half+1)

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def as_real_imag(self, deep=True, **hints):
        if self.args[0].is_real:
            if deep:
                hints['complex'] = False
                return (self.expand(deep, **hints), S.Zero)
            else:
                return (self, S.Zero)
        if deep:
            re, im = self.args[0].expand(deep, **hints).as_real_imag()
        else:
            re, im = self.args[0].as_real_imag()
        return (cos(re)*C.cosh(im), -sin(re)*C.sinh(im))

    def _eval_expand_trig(self, **hints):
        arg = self.args[0]
        x = None
        if arg.is_Add: # TODO: Do this more efficiently for more than two terms
            x, y = arg.as_two_terms()
            sx = sin(x, evaluate=False)._eval_expand_trig()
            sy = sin(y, evaluate=False)._eval_expand_trig()
            cx = cos(x, evaluate=False)._eval_expand_trig()
            cy = cos(y, evaluate=False)._eval_expand_trig()
            return cx*cy - sx*sy
        else:
            coeff, terms = arg.as_coeff_Mul(rational=True)
            if coeff.is_Integer:
                return C.chebyshevt(coeff, cos(terms))
        return cos(arg)

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and C.Order(1, x).contains(arg):
            return S.One
        else:
            return self.func(arg)

    def _eval_is_real(self):
        return self.args[0].is_real

    def _eval_is_bounded(self):
        arg = self.args[0]

        if arg.is_real:
            return True

    def _sage_(self):
        import sage.all as sage
        return sage.cos(self.args[0]._sage_())

class sec(TrigonometricFunction):    #TODO implement rest all functions for sec. see cos, sin, tan.

    def _eval_rewrite_as_cos(self,arg):
        return (1/cos(arg))

    def _eval_rewrite_as_sincos(self, arg):
        return sin(arg)/(cos(arg)*sin(arg))

class csc(TrigonometricFunction):    #TODO implement rest all functions for csc. see cos, sin, tan.

    def _eval_rewrite_as_sin(self,arg):
        return (1/sin(arg))

    def _eval_rewrite_as_sincos(self, arg):
        return cos(arg)/(sin(arg)*cos(arg))

class tan(TrigonometricFunction):
    """
    tan(x) -> Returns the tangent of x (measured in radians)

    Notes
    =====

    * tan(x) will evaluate automatically in the case x is a
      multiple of pi.

    Examples
    ========

    >>> from sympy import tan
    >>> from sympy.abc import x
    >>> tan(x**2).diff(x)
    2*x*(tan(x**2)**2 + 1)
    >>> tan(1).diff(x)
    0

    See Also
    ========

    sin, cos, atan

    References
    ==========

    U{Definitions in trigonometry<http://planetmath.org/encyclopedia/DefinitionsInTrigonometry.html>}

    """

    nargs = 1

    def fdiff(self, argindex=1):
        if argindex==1:
            return S.One + self**2
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        """
        Returns the inverse of this function.
        """
        return atan

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Zero:
                return S.Zero

        if arg.could_extract_minus_sign():
            return -cls(-arg)

        i_coeff = arg.as_coefficient(S.ImaginaryUnit)
        if i_coeff is not None:
            return S.ImaginaryUnit * C.tanh(i_coeff)

        pi_coeff = _pi_coeff(arg, 2)
        if pi_coeff is not None:
            if pi_coeff.is_integer:
                return S.Zero

            if not pi_coeff.is_Rational:
                narg = pi_coeff*S.Pi
                if narg != arg:
                    return cls(narg)
                return None

            cst_table = {
                2 : S.ComplexInfinity,
                3 : sqrt(3),
                4 : S.One,
                6 : 1 / sqrt(3),
            }

            try:
                result = cst_table[pi_coeff.q]

                if (2*pi_coeff.p // pi_coeff.q) % 4 in (1, 3):
                    return -result
                else:
                    return result
            except KeyError:
                if pi_coeff.p > pi_coeff.q:
                    p, q = pi_coeff.p % pi_coeff.q, pi_coeff.q
                    if 2 * p > q:
                        return -cls(Rational(q - p, q)*S.Pi)
                    return cls(Rational(p, q)*S.Pi)
                else:
                    newarg = pi_coeff*S.Pi
                    if newarg != arg:
                        return cls(newarg)
                    return None

        if arg.is_Add:
            x, m = _peeloff_pi(arg)
            if m:
                if (m*2/S.Pi) % 2 == 0:
                    return tan(x)
                else:
                    return -cot(x)

        if arg.func is atan:
            return arg.args[0]

        if arg.func is atan2:
            y, x = arg.args
            return y/x

        if arg.func is asin:
            x = arg.args[0]
            return x / sqrt(1 - x**2)

        if arg.func is acos:
            x = arg.args[0]
            return sqrt(1 - x**2) / x

        if arg.func is acot:
            x = arg.args[0]
            return 1 / x

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)

            a, b = ((n-1)//2), 2**(n+1)

            B = C.bernoulli(n+1)
            F = C.factorial(n+1)

            return (-1)**a * b*(b-1) * B/F * x**n

    def _eval_nseries(self, x, n, logx):
        i = self.args[0].limit(x, 0)*2/S.Pi
        if i and i.is_Integer:
            return self.rewrite(cos)._eval_nseries(x, n=n, logx=logx)
        return Function._eval_nseries(self, x, n=n, logx=logx)

    def _eval_rewrite_as_Pow(self, arg):
        if arg.func is log:
            I = S.ImaginaryUnit
            x = arg.args[0]
            return I*(x**-I - x**I)/(x**-I + x**I)

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def as_real_imag(self, deep=True, **hints):
        if self.args[0].is_real:
            if deep:
                hints['complex'] = False
                return (self.expand(deep, **hints), S.Zero)
            else:
                return (self, S.Zero)
        if deep:
            re, im = self.args[0].expand(deep, **hints).as_real_imag()
        else:
            re, im = self.args[0].as_real_imag()
        denom = cos(re)**2 + C.sinh(im)**2
        return (sin(re)*cos(re)/denom, C.sinh(im)*C.cosh(im)/denom)

    # TODO: Implement _eval_expand_trig

    def _eval_rewrite_as_exp(self, arg):
        exp, I = C.exp, S.ImaginaryUnit
        if isinstance(arg, TrigonometricFunction) or isinstance(arg, HyperbolicFunction) :
            arg = arg.func(arg.args[0]).rewrite(exp)
        neg_exp, pos_exp = exp(-arg*I), exp(arg*I)
        return I*(neg_exp-pos_exp)/(neg_exp+pos_exp)

    def _eval_rewrite_as_sin(self, x):
        return 2*sin(x)**2/sin(2*x)

    def _eval_rewrite_as_cos(self, x):
        return -cos(x + S.Pi/2)/cos(x)

    def _eval_rewrite_as_sincos(self, arg):
        return sin(arg)/cos(arg)

    def _eval_rewrite_as_cot(self, arg):
        return 1/cot(arg)

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and C.Order(1, x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_is_real(self):
        return self.args[0].is_real

    def _eval_is_bounded(self):
        arg = self.args[0]

        if arg.is_imaginary:
            return True

    def _sage_(self):
        import sage.all as sage
        return sage.tan(self.args[0]._sage_())

class cot(TrigonometricFunction):
    """
    cot(x) -> Returns the cotangent of x (measured in radians)
    """

    nargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return S.NegativeOne - self**2
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        """
        Return the inverse of this function.
        """
        return acot

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            if arg is S.Zero:
                return S.ComplexInfinity

        if arg.could_extract_minus_sign():
            return -cls(-arg)

        i_coeff = arg.as_coefficient(S.ImaginaryUnit)
        if i_coeff is not None:
            return -S.ImaginaryUnit * C.coth(i_coeff)

        pi_coeff = _pi_coeff(arg, 2)
        if pi_coeff is not None:
            if pi_coeff.is_integer:
                return S.ComplexInfinity

            if not pi_coeff.is_Rational:
                narg = pi_coeff*S.Pi
                if narg != arg:
                    return cls(narg)
                return None

            cst_table = {
                2 : S.Zero,
                3 : 1 / sqrt(3),
                4 : S.One,
                6 : sqrt(3)
            }

            try:
                result = cst_table[pi_coeff.q]

                if (2*pi_coeff.p // pi_coeff.q) % 4 in (1, 3):
                    return -result
                else:
                    return result
            except KeyError:
                if pi_coeff.p > pi_coeff.q:
                    p, q = pi_coeff.p % pi_coeff.q, pi_coeff.q
                    if 2 * p > q:
                        return -cls(Rational(q - p, q)*S.Pi)
                    return cls(Rational(p, q)*S.Pi)
                else:
                    newarg = pi_coeff*S.Pi
                    if newarg != arg:
                        return cls(newarg)
                    return None

        if arg.is_Add:
            x, m = _peeloff_pi(arg)
            if m:
                if (m*2/S.Pi) % 2 == 0:
                    return cot(x)
                else:
                    return -tan(x)

        if arg.func is acot:
            return arg.args[0]

        if arg.func is atan:
            x = arg.args[0]
            return 1 / x

        if arg.func is atan2:
            y, x = arg.args
            return x/y

        if arg.func is asin:
            x = arg.args[0]
            return sqrt(1 - x**2) / x

        if arg.func is acos:
            x = arg.args[0]
            return x / sqrt(1 - x**2)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n == 0:
            return 1 / sympify(x)
        elif n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)

            B = C.bernoulli(n+1)
            F = C.factorial(n+1)

            return (-1)**((n+1)//2) * 2**(n+1) * B/F * x**n

    def _eval_nseries(self, x, n, logx):
        i = self.args[0].limit(x, 0)/S.Pi
        if i and i.is_Integer:
            return self.rewrite(cos)._eval_nseries(x, n=n, logx=logx)
        return Function._eval_nseries(self, x, n=n, logx=logx)

    def _eval_conjugate(self):
        assert len(self.args) == 1
        return self.func(self.args[0].conjugate())

    def as_real_imag(self, deep=True, **hints):
        if self.args[0].is_real:
            if deep:
                hints['complex'] = False
                return (self.expand(deep, **hints), S.Zero)
            else:
                return (self, S.Zero)
        if deep:
            re, im = self.args[0].expand(deep, **hints).as_real_imag()
        else:
            re, im = self.args[0].as_real_imag()
        denom = sin(re)**2 + C.sinh(im)**2
        return (sin(re)*cos(re)/denom, -C.sinh(im)*C.cosh(im)/denom)

    def _eval_rewrite_as_exp(self, arg):
        exp, I = C.exp, S.ImaginaryUnit
        if isinstance(arg, TrigonometricFunction) or isinstance(arg, HyperbolicFunction) :
            arg = arg.func(arg.args[0]).rewrite(exp)
        neg_exp, pos_exp = exp(-arg*I), exp(arg*I)
        return I*(pos_exp+neg_exp)/(pos_exp-neg_exp)

    def _eval_rewrite_as_Pow(self, arg):
        if arg.func is log:
            I = S.ImaginaryUnit
            x = arg.args[0]
            return -I*(x**-I + x**I)/(x**-I - x**I)

    def _eval_rewrite_as_sin(self, x):
        return 2*sin(2*x)/sin(x)**2

    def _eval_rewrite_as_cos(self, x):
        return -cos(x)/cos(x + S.Pi/2)

    def _eval_rewrite_as_sincos(self, arg):
        return cos(arg)/sin(arg)

    def _eval_rewrite_as_tan(self, arg):
        return 1/tan(arg)

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and C.Order(1, x).contains(arg):
            return 1/arg
        else:
            return self.func(arg)

    def _eval_is_real(self):
        return self.args[0].is_real

    def _sage_(self):
        import sage.all as sage
        return sage.cot(self.args[0]._sage_())

###############################################################################
########################### TRIGONOMETRIC INVERSES ############################
###############################################################################

class asin(Function):
    """
    asin(x) -> Returns the arc sine of x (measured in radians)

    Notes
    =====

    * asin(x) will evaluate automatically in the cases
      oo, -oo, 0, 1, -1

    Examples
    ========

    >>> from sympy import asin, oo, pi
    >>> asin(1)
    pi/2
    >>> asin(-1)
    -pi/2

    See Also
    ========

    acos, atan, sin
    """

    nargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/sqrt(1 - self.args[0]**2)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.NegativeInfinity * S.ImaginaryUnit
            elif arg is S.NegativeInfinity:
                return S.Infinity * S.ImaginaryUnit
            elif arg is S.Zero:
                return S.Zero
            elif arg is S.One:
                return S.Pi / 2
            elif arg is S.NegativeOne:
                return -S.Pi / 2

        if arg.could_extract_minus_sign():
            return -cls(-arg)

        if arg.is_number:
            cst_table = {
                sqrt(3)/2  : 3,
                -sqrt(3)/2 : -3,
                sqrt(2)/2  : 4,
                -sqrt(2)/2 : -4,
                1/sqrt(2)  : 4,
                -1/sqrt(2) : -4,
                sqrt((5-sqrt(5))/8) : 5,
                -sqrt((5-sqrt(5))/8) : -5,
                S.Half     : 6,
                -S.Half    : -6,
                sqrt(2-sqrt(2))/2 : 8,
                -sqrt(2-sqrt(2))/2 : -8,
                (sqrt(5)-1)/4 : 10,
                (1-sqrt(5))/4 : -10,
                (sqrt(3)-1)/sqrt(2**3) : 12,
                (1-sqrt(3))/sqrt(2**3) : -12,
                (sqrt(5)+1)/4 : S(10)/3,
                -(sqrt(5)+1)/4 : -S(10)/3
                }

            if arg in cst_table:
                return S.Pi / cst_table[arg]

        i_coeff = arg.as_coefficient(S.ImaginaryUnit)
        if i_coeff is not None:
            return S.ImaginaryUnit * C.asinh(i_coeff)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)
            if len(previous_terms) >= 2 and n > 2:
                p = previous_terms[-2]
                return p * (n-2)**2/(n*(n-1)) * x**2
            else:
                k = (n - 1) // 2
                R = C.RisingFactorial(S.Half, k)
                F = C.factorial(k)
                return R / F * x**n / n

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and C.Order(1, x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_rewrite_as_acos(self, x):
        return S.Pi/2 - acos(x)

    def _eval_rewrite_as_atan(self, x):
        return 2*atan(x/(1 + sqrt(1 - x**2)))

    def _eval_rewrite_as_log(self, x):
        return -S.ImaginaryUnit*C.log(S.ImaginaryUnit*x + sqrt(1-x**2))

    def _eval_is_real(self):
        return self.args[0].is_real and (self.args[0]>=-1 and self.args[0]<=1)

    def _sage_(self):
        import sage.all as sage
        return sage.asin(self.args[0]._sage_())

class acos(Function):
    """
    acos(x) -> Returns the arc cosine of x (measured in radians)

    Notes
    =====

    * acos(x) will evaluate automatically in the cases
      oo, -oo, 0, 1, -1

    Examples
    ========

    >>> from sympy import acos, oo, pi
    >>> acos(1)
    0
    >>> acos(0)
    pi/2
    >>> acos(oo)
    oo*I

    See Also
    ========

    asin, atan, cos
    """

    nargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -1/sqrt(1 - self.args[0]**2)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.Infinity * S.ImaginaryUnit
            elif arg is S.NegativeInfinity:
                return S.NegativeInfinity * S.ImaginaryUnit
            elif arg is S.Zero:
                return S.Pi / 2
            elif arg is S.One:
                return S.Zero
            elif arg is S.NegativeOne:
                return S.Pi

        if arg.is_number:
            cst_table = {
                S.Half     : S.Pi/3,
                -S.Half    : 2*S.Pi/3,
                sqrt(2)/2  : S.Pi/4,
                -sqrt(2)/2 : 3*S.Pi/4,
                1/sqrt(2)  : S.Pi/4,
                -1/sqrt(2) : 3*S.Pi/4,
                sqrt(3)/2  : S.Pi/6,
                -sqrt(3)/2 : 5*S.Pi/6,
                }

            if arg in cst_table:
                return cst_table[arg]


    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n == 0:
            return S.Pi / 2
        elif n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)
            if len(previous_terms) >= 2 and n > 2:
                p = previous_terms[-2]
                return p * (n-2)**2/(n*(n-1)) * x**2
            else:
                k = (n - 1) // 2
                R = C.RisingFactorial(S.Half, k)
                F = C.factorial(k)
                return -R / F * x**n / n

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and C.Order(1, x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_is_real(self):
        return self.args[0].is_real and (self.args[0]>=-1 and self.args[0]<=1)

    def _eval_rewrite_as_log(self, x):
        return S.Pi/2 + S.ImaginaryUnit * C.log(S.ImaginaryUnit * x + sqrt(1 - x**2))

    def _eval_rewrite_as_asin(self, x):
        return S.Pi/2 - asin(x)

    def _eval_rewrite_as_atan(self, x):
        if x > -1 and x <= 1:
            return 2 * atan(sqrt(1 - x**2)/(1 + x))
        else:
            raise ValueError("The argument must be bounded in the interval (-1,1]")

    def _sage_(self):
        import sage.all as sage
        return sage.acos(self.args[0]._sage_())

class atan(Function):
    """
    atan(x) -> Returns the arc tangent of x (measured in radians)

    Notes
    =====

    * atan(x) will evaluate automatically in the cases
      oo, -oo, 0, 1, -1

    Examples
    ========

    >>> from sympy import atan, oo, pi
    >>> atan(0)
    0
    >>> atan(1)
    pi/4
    >>> atan(oo)
    pi/2

    See Also
    ========

    acos, asin, tan
    """

    nargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/(1+self.args[0]**2)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.Pi / 2
            elif arg is S.NegativeInfinity:
                return -S.Pi / 2
            elif arg is S.Zero:
                return S.Zero
            elif arg is S.One:
                return S.Pi / 4
            elif arg is S.NegativeOne:
                return -S.Pi / 4
        if arg.could_extract_minus_sign():
            return -cls(-arg)

        if arg.is_number:
            cst_table = {
                sqrt(3)/3  : 6,
                -sqrt(3)/3 : -6,
                1/sqrt(3)  : 6,
                -1/sqrt(3) : -6,
                sqrt(3)    : 3,
                -sqrt(3)   : -3,
                (1+sqrt(2)) : S(8)/3,
                -(1+sqrt(2)) : S(8)/3,
                (sqrt(2)-1) : 8,
                (1-sqrt(2)) : -8,
                sqrt((5+2*sqrt(5))) : S(5)/2,
                -sqrt((5+2*sqrt(5))) : -S(5)/2,
                (2-sqrt(3)) : 12,
                -(2-sqrt(3)) : -12
                }

            if arg in cst_table:
                return S.Pi / cst_table[arg]

        i_coeff = arg.as_coefficient(S.ImaginaryUnit)
        if i_coeff is not None:
            return S.ImaginaryUnit * C.atanh(i_coeff)



    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)
            return (-1)**((n-1)//2) * x**n / n

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and C.Order(1, x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_is_real(self):
        return self.args[0].is_real

    def _eval_rewrite_as_log(self, x):
        return S.ImaginaryUnit/2 * \
               (C.log((S(1) - S.ImaginaryUnit * x)/(S(1) + S.ImaginaryUnit * x)))

    def _eval_aseries(self, n, args0, x, logx):
        if args0[0] == S.Infinity:
            return S.Pi/2 - atan(1/self.args[0])
        elif args0[0] == S.NegativeInfinity:
            return -S.Pi/2 - atan(1/self.args[0])
        else:
            return super(atan, self)._eval_aseries(n, args0, x, logx)

    def _sage_(self):
        import sage.all as sage
        return sage.atan(self.args[0]._sage_())

class acot(Function):
    """
    acot(x) -> Returns the arc cotangent of x (measured in radians)
    """

    nargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -1 / (1+self.args[0]**2)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.Zero
            elif arg is S.NegativeInfinity:
                return S.Zero
            elif arg is S.Zero:
                return S.Pi/ 2
            elif arg is S.One:
                return S.Pi / 4
            elif arg is S.NegativeOne:
                return -S.Pi / 4

        if arg.could_extract_minus_sign():
            return -cls(-arg)

        if arg.is_number:
            cst_table = {
                sqrt(3)/3  : 3,
                -sqrt(3)/3 : -3,
                1/sqrt(3)  : 3,
                -1/sqrt(3) : -3,
                sqrt(3)    : 6,
                -sqrt(3)   : -6,
                (1+sqrt(2)) : 8,
                -(1+sqrt(2)) : -8,
                (1-sqrt(2)) : -S(8)/3,
                (sqrt(2)-1) : S(8)/3,
                sqrt(5+2*sqrt(5)) : 10,
                -sqrt(5+2*sqrt(5)) : -10,
                (2+sqrt(3)) : 12,
                -(2+sqrt(3)) : -12,
                (2-sqrt(3)) : S(12)/5,
                -(2-sqrt(3)) : -S(12)/5,
                }

            if arg in cst_table:
                return S.Pi / cst_table[arg]

        i_coeff = arg.as_coefficient(S.ImaginaryUnit)
        if i_coeff is not None:
            return -S.ImaginaryUnit * C.acoth(i_coeff)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n == 0:
            return S.Pi / 2 # FIX THIS
        elif n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)
            return (-1)**((n+1)//2) * x**n / n

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and C.Order(1, x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_is_real(self):
        return self.args[0].is_real

    def _eval_aseries(self, n, args0, x, logx):
        if args0[0] == S.Infinity:
            return S.Pi/2 - acot(1/self.args[0])
        elif args0[0] == S.NegativeInfinity:
            return 3*S.Pi/2 - acot(1/self.args[0])
        else:
            return super(atan, self)._eval_aseries(n, args0, x, logx)

    def _sage_(self):
        import sage.all as sage
        return sage.acot(self.args[0]._sage_())

    def _eval_rewrite_as_log(self, x):
        return S.ImaginaryUnit/2 * \
               (C.log((x - S.ImaginaryUnit)/(x + S.ImaginaryUnit)))

class atan2(Function):
    """
    atan2(y,x) -> Returns the atan(y/x) taking two arguments y and x.
    Signs of both y and x are considered to determine the appropriate
    quadrant of atan(y/x). The range is (-pi, pi].
    """

    nargs = 2

    @classmethod
    def eval(cls, y, x):
        sign_y = C.sign(y)

        if y.is_zero:
            if x.is_positive:
                return S.Zero
            elif x.is_zero:
                return S.NaN
            elif x.is_negative:
                return S.Pi
        elif x.is_zero:
            if sign_y.is_Number:
                return sign_y * S.Pi/2
        elif x.is_zero is False:
            abs_yx = C.Abs(y/x)
            if sign_y.is_Number and abs_yx.is_number:
                phi = C.atan(abs_yx)
                if x.is_positive:
                    return sign_y * phi
                else:
                    return sign_y * (S.Pi - phi)

    def _eval_is_real(self):
        return self.args[0].is_real and self.args[1].is_real

    def fdiff(self, argindex):
        x, y = self.args
        if argindex == 1:
            return y/(x**2 + y**2)
        elif argindex == 2:
            return -x/(x**2 + y**2)
        else:
            raise ArgumentIndexError(self, argindex)

    def _sage_(self):
        import sage.all as sage
        return sage.atan2(self.args[0]._sage_(), self.args[1]._sage_())

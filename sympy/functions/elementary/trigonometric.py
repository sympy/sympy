from __future__ import print_function, division

from sympy.core.add import Add
from sympy.core.basic import C, sympify, cacheit
from sympy.core.singleton import S
from sympy.core.numbers import igcdex
from sympy.core.function import Function, ArgumentIndexError
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.exponential import log
from sympy.functions.elementary.hyperbolic import HyperbolicFunction
from sympy.utilities.iterables import numbered_symbols
from sympy.core.compatibility import xrange

###############################################################################
########################## TRIGONOMETRIC FUNCTIONS ############################
###############################################################################


class TrigonometricFunction(Function):
    """Base class for trigonometric functions. """

    unbranched = True

    def _eval_is_rational(self):
        s = self.func(*self.args)
        if s.func == self.func:
            if s.args[0].is_rational:
                return False
        else:
            return s.is_rational


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

    Examples
    ========

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
                        c = C.Rational(i, m)
                        cx = c*x
                else:
                    c = C.Rational(int(c))
                    cx = c*x
            if x.is_integer:
                c2 = c % 2
                if c2 == 1:
                    return x
                elif not c2:
                    if x.is_even is not None:  # known parity
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

    .. [1] http://planetmath.org/encyclopedia/DefinitionsInTrigonometry.html

    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return cos(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Zero:
                return S.Zero
            elif arg is S.Infinity or arg is S.NegativeInfinity:
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

            if (2*pi_coeff).is_integer:
                return S.NegativeOne**(pi_coeff - S.Half)

            if not pi_coeff.is_Rational:
                narg = pi_coeff*S.Pi
                if narg != arg:
                    return cls(narg)
                return None

            # http://code.google.com/p/sympy/issues/detail?id=2949
            # transform a sine to a cosine, to avoid redundant code
            if pi_coeff.is_Rational:
                x = pi_coeff % 2
                if x > 1:
                    return -cls((x % 1)*S.Pi)
                if 2*x > 1:
                    return cls((1 - x)*S.Pi)
                narg = ((pi_coeff + C.Rational(3, 2)) % 2)*S.Pi
                result = cos(narg)
                if not isinstance(result, cos):
                    return result
                if pi_coeff*S.Pi != arg:
                    return cls(pi_coeff*S.Pi)
                return None

        if arg.is_Add:
            x, m = _peeloff_pi(arg)
            if m:
                return sin(m)*cos(x) + cos(m)*sin(x)

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
                return -p * x**2 / (n*(n - 1))
            else:
                return (-1)**(n//2) * x**(n)/C.factorial(n)

    def _eval_rewrite_as_exp(self, arg):
        exp, I = C.exp, S.ImaginaryUnit
        if isinstance(arg, TrigonometricFunction) or isinstance(arg, HyperbolicFunction):
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

    def _eval_rewrite_as_pow(self, arg):
        return self.rewrite(cos).rewrite(pow)

    def _eval_rewrite_as_sqrt(self, arg):
        return self.rewrite(cos).rewrite(sqrt)

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
        if arg.is_Add:  # TODO, implement more if deep stuff here
            # TODO: Do this more efficiently for more than two terms
            x, y = arg.as_two_terms()
            sx = sin(x, evaluate=False)._eval_expand_trig()
            sy = sin(y, evaluate=False)._eval_expand_trig()
            cx = cos(x, evaluate=False)._eval_expand_trig()
            cy = cos(y, evaluate=False)._eval_expand_trig()
            return sx*cy + sy*cx
        else:
            n, x = arg.as_coeff_Mul(rational=True)
            if n.is_Integer:  # n will be positive because of .eval
                # canonicalization

                # See http://mathworld.wolfram.com/Multiple-AngleFormulas.html
                if n.is_odd:
                    return (-1)**((n - 1)/2)*C.chebyshevt(n, sin(x))
                else:
                    return expand_mul((-1)**(n/2 - 1)*cos(x)*C.chebyshevu(n -
                        1, sin(x)), deep=False)
            pi_coeff = _pi_coeff(arg)
            if pi_coeff is not None:
                if pi_coeff.is_Rational:
                    return self.rewrite(sqrt)
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

    .. [1] http://planetmath.org/encyclopedia/DefinitionsInTrigonometry.html

    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -sin(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Zero:
                return S.One
            elif arg is S.Infinity or arg is S.NegativeInfinity:
                # In this cases, it is unclear if we should
                # return S.NaN or leave un-evaluated.  One
                # useful test case is how "limit(sin(x)/x,x,oo)"
                # is handled.
                # See test_sin_cos_with_infinity() an
                # Test for issue 209
                # http://code.google.com/p/sympy/issues/detail?id=2097
                # For now, we return un-evaluated.
                return

        if arg.could_extract_minus_sign():
            return cls(-arg)

        i_coeff = arg.as_coefficient(S.ImaginaryUnit)
        if i_coeff is not None:
            return C.cosh(i_coeff)

        pi_coeff = _pi_coeff(arg)
        if pi_coeff is not None:
            if pi_coeff.is_integer:
                return (S.NegativeOne)**pi_coeff

            if (2*pi_coeff).is_integer:
                return S.Zero

            if not pi_coeff.is_Rational:
                narg = pi_coeff*S.Pi
                if narg != arg:
                    return cls(narg)
                return None

            # cosine formula #####################
            # http://code.google.com/p/sympy/issues/detail?id=2949
            # explicit calculations are preformed for
            # cos(k pi / 8), cos(k pi /10), and cos(k pi / 12)
            # Some other exact values like cos(k pi/15) can be
            # calculated using a partial-fraction decomposition
            # by calling cos( X ).rewrite(sqrt)
            cst_table_some = {
                3: S.Half,
                5: (sqrt(5) + 1)/4,
            }
            if pi_coeff.is_Rational:
                q = pi_coeff.q
                p = pi_coeff.p % (2*q)
                if p > q:
                    narg = (pi_coeff - 1)*S.Pi
                    return -cls(narg)
                if 2*p > q:
                    narg = (1 - pi_coeff)*S.Pi
                    return -cls(narg)

                # If nested sqrt's are worse than un-evaluation
                # you can require q in (1, 2, 3, 4, 6)
                # q <= 12 returns expressions with 2 or fewer nestings.
                if q > 12:
                    return None

                if q in cst_table_some:
                    cts = cst_table_some[pi_coeff.q]
                    return C.chebyshevt(pi_coeff.p, cts).expand()

                if 0 == q % 2:
                    narg = (pi_coeff*2)*S.Pi
                    nval = cls(narg)
                    if None == nval:
                        return None
                    x = (2*pi_coeff + 1)/2
                    sign_cos = (-1)**((-1 if x < 0 else 1)*int(abs(x)))
                    return sign_cos*sqrt( (1 + nval)/2 )
            return None

        if arg.is_Add:
            x, m = _peeloff_pi(arg)
            if m:
                return cos(m)*cos(x) - sin(m)*sin(x)

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
                return -p * x**2 / (n*(n - 1))
            else:
                return (-1)**(n//2)*x**(n)/C.factorial(n)

    def _eval_rewrite_as_exp(self, arg):
        exp, I = C.exp, S.ImaginaryUnit
        if isinstance(arg, TrigonometricFunction) or isinstance(arg, HyperbolicFunction):
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
        return (1 - tan_half)/(1 + tan_half)

    def _eval_rewrite_as_sincos(self, arg):
        return sin(arg)*cos(arg)/sin(arg)

    def _eval_rewrite_as_cot(self, arg):
        cot_half = cot(S.Half*arg)**2
        return (cot_half - 1)/(cot_half + 1)

    def _eval_rewrite_as_pow(self, arg):
        return self._eval_rewrite_as_sqrt(arg)

    def _eval_rewrite_as_sqrt(self, arg):
        _EXPAND_INTS = False

        def migcdex(x):
            # recursive calcuation of gcd and linear combination
            # for a sequence of integers.
            # Given  (x1, x2, x3)
            # Returns (y1, y1, y3, g)
            # such that g is the gcd and x1*y1+x2*y2+x3*y3 - g = 0
            # Note, that this is only one such linear combination.
            if len(x) == 1:
                return (1, x[0])
            if len(x) == 2:
                return igcdex(x[0], x[-1])
            g = migcdex(x[1:])
            u, v, h = igcdex(x[0], g[-1])
            return tuple([u] + [v*i for i in g[0:-1] ] + [h])

        def ipartfrac(r, factors=None):
            if isinstance(r, int):
                return r
            if not isinstance(r, C.Rational):
                raise TypeError("r is not rational")
            n = r.q
            if 2 > r.q*r.q:
                return r.q

            if None == factors:
                a = [n//x**y for x, y in factorint(r.q).items()]
            else:
                a = [n//x for x in factors]
            if len(a) == 1:
                return [ r ]
            h = migcdex(a)
            ans = [ r.p*C.Rational(i*j, r.q) for i, j in zip(h[:-1], a) ]
            assert r == sum(ans)
            return ans
        pi_coeff = _pi_coeff(arg)
        if pi_coeff is None:
            return None

        assert not pi_coeff.is_integer, "should have been simplified already"

        if not pi_coeff.is_Rational:
            return None

        cst_table_some = {
            3: S.Half,
            5: (sqrt(5) + 1)/4,
            17: sqrt((15 + sqrt(17))/32 + sqrt(2)*(sqrt(17 - sqrt(17)) +
                sqrt(sqrt(2)*(-8*sqrt(17 + sqrt(17)) - (1 - sqrt(17))
                *sqrt(17 - sqrt(17))) + 6*sqrt(17) + 34))/32)
            # 65537 and 257 are the only other known Fermat primes
            # Please add if you would like them
        }

        def fermatCoords(n):
            if not isinstance(n, int):
                raise TypeError("n is not an integer")
            if n <= 0:
                raise ValueError("n has to be greater than 0")
            if n == 1 or 0 == n % 2:
                return False
            primes = dict( [(p, 0) for p in cst_table_some ] )
            assert 1 not in primes
            for p_i in primes:
                while 0 == n % p_i:
                    n = n/p_i
                    primes[p_i] += 1
            if 1 != n:
                return False
            if max(primes.values()) > 1:
                return False
            return tuple([ p for p in primes if primes[p] == 1])

        if pi_coeff.q in cst_table_some:
            return C.chebyshevt(pi_coeff.p, cst_table_some[pi_coeff.q]).expand()

        if 0 == pi_coeff.q % 2:  # recursively remove powers of 2
            narg = (pi_coeff*2)*S.Pi
            nval = cos(narg)
            if None == nval:
                return None
            nval = nval.rewrite(sqrt)
            if not _EXPAND_INTS:
                if (isinstance(nval, cos) or isinstance(-nval, cos)):
                    return None
            x = (2*pi_coeff + 1)/2
            sign_cos = (-1)**((-1 if x < 0 else 1)*int(abs(x)))
            return sign_cos*sqrt( (1 + nval)/2 )

        FC = fermatCoords(pi_coeff.q)
        if FC:
            decomp = ipartfrac(pi_coeff, FC)
            X = [(x[1], x[0]*S.Pi) for x in zip(decomp, numbered_symbols('z'))]
            pcls = cos(sum([x[0] for x in X]))._eval_expand_trig().subs(X)
            return pcls.rewrite(sqrt)
        if _EXPAND_INTS:
            decomp = ipartfrac(pi_coeff)
            X = [(x[1], x[0]*S.Pi) for x in zip(decomp, numbered_symbols('z'))]
            pcls = cos(sum([x[0] for x in X]))._eval_expand_trig().subs(X)
            return pcls
        return None

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
        if arg.is_Add:  # TODO: Do this more efficiently for more than two terms
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
            pi_coeff = _pi_coeff(arg)
            if pi_coeff is not None:
                if pi_coeff.is_Rational:
                    return self.rewrite(sqrt)
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


class ReciprocalTrigonometricFunction(TrigonometricFunction):
    """Base class for reciprocal functions of trigonometric functions. """

    _reciprocal_of = None       # mandatory, to be defined in subclass

    # _is_even and _is_odd are used for correct evaluation of csc(-x), sec(-x)
    # TODO refactor into TrigonometricFunction common parts of
    # trigonometric functions eval() like even/odd, func(x+2*k*pi), etc.
    _is_even = None  # optional, to be defined in subclass
    _is_odd = None   # optional, to be defined in subclass

    def _call_reciprocal(self, method_name, *args, **kwargs):
        # Calls method_name on _reciprocal_of
        o = self._reciprocal_of(self.args[0])
        if kwargs:
            return getattr(o, method_name)(**kwargs)
        else:
            return getattr(o, method_name)(*args)

    def _calculate_reciprocal(self, method_name, *args, **kwargs):
        # If calling method_name on _reciprocal_of returns a value != None
        # then return the reciprocal of that value
        t = self._call_reciprocal(method_name, *args, **kwargs)
        return 1/t if t != None else t

    def _rewrite_reciprocal(self, method_name, arg):
        # Special handling for rewrite functions. If reciprocal rewrite returns
        # unmodified expression, then return None
        t = self._call_reciprocal(method_name, arg)
        if t != None and t != self._reciprocal_of(arg):
            return 1/t
        else:
            return

    def fdiff(self, argindex=1):
        return self._calculate_reciprocal("fdiff", argindex)

    def _eval_rewrite_as_exp(self, arg):
        return self._rewrite_reciprocal("_eval_rewrite_as_exp", arg)

    def _eval_rewrite_as_Pow(self, arg):
        return self._rewrite_reciprocal("_eval_rewrite_as_Pow", arg)

    def _eval_rewrite_as_sin(self, arg):
        return self._rewrite_reciprocal("_eval_rewrite_as_sin", arg)

    def _eval_rewrite_as_cos(self, arg):
        return self._rewrite_reciprocal("_eval_rewrite_as_cos", arg)

    def _eval_rewrite_as_tan(self, arg):
        return self._rewrite_reciprocal("_eval_rewrite_as_tan", arg)

    def _eval_rewrite_as_pow(self, arg):
        return self._rewrite_reciprocal("_eval_rewrite_as_pow", arg)

    def _eval_rewrite_as_sqrt(self, arg):
        return self._rewrite_reciprocal("_eval_rewrite_as_sqrt", arg)

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def as_real_imag(self, deep=True, **hints):
        return (1/self._reciprocal_of(self.args[0])).as_real_imag(deep,
                                                                  **hints)

    def _eval_expand_trig(self, **hints):
        return self._calculate_reciprocal("_eval_expand_trig", **hints)

    def _eval_is_real(self):
        return self._reciprocal_of(self.args[0])._eval_is_real()

    def _eval_as_leading_term(self, x):
        return (1/self._reciprocal_of(self.args[0]))._eval_as_leading_term(x)

    def _eval_is_bounded(self):
        return (1/self._reciprocal_of(self.args[0])).is_bounded

    def _eval_nseries(self, x, n, logx):
        return (1/self._reciprocal_of(self.args[0]))._eval_nseries(x, n, logx)

    @classmethod
    def eval(cls, arg):
        if arg.could_extract_minus_sign():
            if cls._is_even:
                return cls(-arg)
            if cls._is_odd:
                return -cls(-arg)

        pi_coeff = _pi_coeff(arg)
        if (pi_coeff is not None
            and not (2*pi_coeff).is_integer
            and pi_coeff.is_Rational):
                q = pi_coeff.q
                p = pi_coeff.p % (2*q)
                if p > q:
                    narg = (pi_coeff - 1)*S.Pi
                    return -cls(narg)
                if 2*p > q:
                    narg = (1 - pi_coeff)*S.Pi
                    return -cls(narg)
        t = cls._reciprocal_of.eval(arg)
        return 1/t if t != None else t


class sec(ReciprocalTrigonometricFunction):
    _reciprocal_of = cos
    _is_even = True


    def _eval_rewrite_as_cos(self, arg):
        return (1/cos(arg))

    def _eval_rewrite_as_sincos(self, arg):
        return sin(arg)/(cos(arg)*sin(arg))

    def fdiff(self, argindex=1):
        if argindex == 1:
            return tan(self.args[0])*sec(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    # TODO def taylor_term(n, x, *previous_terms):

    def _sage_(self):
        import sage.all as sage
        return sage.sec(self.args[0]._sage_())


class csc(ReciprocalTrigonometricFunction):
    _reciprocal_of = sin
    _is_odd = True


    def _eval_rewrite_as_sin(self, arg):
        return (1/sin(arg))

    def _eval_rewrite_as_sincos(self, arg):
        return cos(arg)/(sin(arg)*cos(arg))

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -cot(self.args[0])*csc(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    # TODO def taylor_term(n, x, *previous_terms):

    def _sage_(self):
        import sage.all as sage
        return sage.csc(self.args[0]._sage_())


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

    .. [1] http://planetmath.org/encyclopedia/DefinitionsInTrigonometry.html

    """

    def fdiff(self, argindex=1):
        if argindex == 1:
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

            if pi_coeff.is_Rational:
                narg = ((pi_coeff + S.Half) % 1 - S.Half)*S.Pi
                # see cos() to specify which expressions should  be
                # expanded automatically in terms of radicals
                cresult, sresult = cos(narg), cos(narg - S.Pi/2)
                if not isinstance(cresult, cos) \
                        and not isinstance(sresult, cos):
                    if cresult == 0:
                        return S.ComplexInfinity
                    return (sresult/cresult)
                if narg != arg:
                    return cls(narg)

        if arg.is_Add:
            x, m = _peeloff_pi(arg)
            if m:
                tanm = tan(m)
                tanx = tan(x)
                if tanm is S.ComplexInfinity:
                    return -cot(x)
                return (tanm + tanx)/(1 - tanm*tanx)

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

            a, b = ((n - 1)//2), 2**(n + 1)

            B = C.bernoulli(n + 1)
            F = C.factorial(n + 1)

            return (-1)**a * b*(b - 1) * B/F * x**n

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

    def _eval_expand_trig(self, **hints):
        arg = self.args[0]
        x = None
        if arg.is_Add:
            from sympy import symmetric_poly
            n = len(arg.args)
            TX = []
            for x in arg.args:
                tx = tan(x, evaluate=False)._eval_expand_trig()
                TX.append(tx)

            Yg = numbered_symbols('Y')
            Y = [ next(Yg) for i in xrange(n) ]

            p = [0, 0]
            for i in xrange(n + 1):
                p[1 - i % 2] += symmetric_poly(i, Y)*(-1)**((i % 4)//2)
            return (p[0]/p[1]).subs(list(zip(Y, TX)))

        else:
            coeff, terms = arg.as_coeff_Mul(rational=True)
            if coeff.is_Integer and coeff > 1:
                I = S.ImaginaryUnit
                z = C.Symbol('dummy', real=True)
                P = ((1 + I*z)**coeff).expand()
                return (C.im(P)/C.re(P)).subs([(z, tan(terms))])
        return tan(arg)

    def _eval_rewrite_as_exp(self, arg):
        exp, I = C.exp, S.ImaginaryUnit
        if isinstance(arg, TrigonometricFunction) or isinstance(arg, HyperbolicFunction):
            arg = arg.func(arg.args[0]).rewrite(exp)
        neg_exp, pos_exp = exp(-arg*I), exp(arg*I)
        return I*(neg_exp - pos_exp)/(neg_exp + pos_exp)

    def _eval_rewrite_as_sin(self, x):
        return 2*sin(x)**2/sin(2*x)

    def _eval_rewrite_as_cos(self, x):
        return -cos(x + S.Pi/2)/cos(x)

    def _eval_rewrite_as_sincos(self, arg):
        return sin(arg)/cos(arg)

    def _eval_rewrite_as_cot(self, arg):
        return 1/cot(arg)

    def _eval_rewrite_as_pow(self, arg):
        y = self.rewrite(cos).rewrite(pow)
        if y.has(cos):
            return None
        return y

    def _eval_rewrite_as_sqrt(self, arg):
        y = self.rewrite(cos).rewrite(sqrt)
        if y.has(cos):
            return None
        return y

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

    def fdiff(self, argindex=1):
        if argindex == 1:
            return S.NegativeOne - self**2
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        """
        Returns the inverse of this function.
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

            if pi_coeff.is_Rational:
                narg = (((pi_coeff + S.Half) % 1) - S.Half)*S.Pi
                # see cos() to specify which expressions should be
                # expanded automatically in terms of radicals
                cresult, sresult = cos(narg), cos(narg - S.Pi/2)
                if not isinstance(cresult, cos) \
                        and not isinstance(sresult, cos):
                    if sresult == 0:
                        return S.ComplexInfinity
                    return cresult / sresult
                if narg != arg:
                    return cls(narg)

        if arg.is_Add:
            x, m = _peeloff_pi(arg)
            if m:
                cotm = cot(m)
                if cotm == 0:
                    return -tan(x)
                cotx = cot(x)
                if cotm is S.ComplexInfinity:
                    return cotx
                if cotm.is_Rational:
                    return (cotm*cotx - 1) / (cotm + cotx)
            return None

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

            B = C.bernoulli(n + 1)
            F = C.factorial(n + 1)

            return (-1)**((n + 1)//2) * 2**(n + 1) * B/F * x**n

    def _eval_nseries(self, x, n, logx):
        i = self.args[0].limit(x, 0)/S.Pi
        if i and i.is_Integer:
            return self.rewrite(cos)._eval_nseries(x, n=n, logx=logx)
        return self.rewrite(tan)._eval_nseries(x, n=n, logx=logx)

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
        if isinstance(arg, TrigonometricFunction) or isinstance(arg, HyperbolicFunction):
            arg = arg.func(arg.args[0]).rewrite(exp)
        neg_exp, pos_exp = exp(-arg*I), exp(arg*I)
        return I*(pos_exp + neg_exp)/(pos_exp - neg_exp)

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

    def _eval_rewrite_as_pow(self, arg):
        y = self.rewrite(cos).rewrite(pow)
        if y.has(cos):
            return None
        return y

    def _eval_rewrite_as_sqrt(self, arg):
        y = self.rewrite(cos).rewrite(sqrt)
        if y.has(cos):
            return None
        return y

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and C.Order(1, x).contains(arg):
            return 1/arg
        else:
            return self.func(arg)

    def _eval_is_real(self):
        return self.args[0].is_real

    def _eval_expand_trig(self, **hints):
        arg = self.args[0]
        x = None
        if arg.is_Add:
            from sympy import symmetric_poly
            n = len(arg.args)
            CX = []
            for x in arg.args:
                cx = cot(x, evaluate=False)._eval_expand_trig()
                CX.append(cx)

            Yg = numbered_symbols('Y')
            Y = [ next(Yg) for i in xrange(n) ]

            p = [0, 0]
            for i in xrange(n, -1, -1):
                p[(n - i) % 2] += symmetric_poly(i, Y)*(-1)**(((n - i) % 4)//2)
            return (p[0]/p[1]).subs(list(zip(Y, CX)))
        else:
            coeff, terms = arg.as_coeff_Mul(rational=True)
            if coeff.is_Integer and coeff > 1:
                I = S.ImaginaryUnit
                z = C.Symbol('dummy', real=True)
                P = ((z + I)**coeff).expand()
                return (C.re(P)/C.im(P)).subs([(z, cot(terms))])
        return cot(arg)

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

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/sqrt(1 - self.args[0]**2)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_is_rational(self):
        s = self.func(*self.args)
        if s.func == self.func:
            if s.args[0].is_rational:
                return False
        else:
            return s.is_rational

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
                sqrt(3)/2: 3,
                -sqrt(3)/2: -3,
                sqrt(2)/2: 4,
                -sqrt(2)/2: -4,
                1/sqrt(2): 4,
                -1/sqrt(2): -4,
                sqrt((5 - sqrt(5))/8): 5,
                -sqrt((5 - sqrt(5))/8): -5,
                S.Half: 6,
                -S.Half: -6,
                sqrt(2 - sqrt(2))/2: 8,
                -sqrt(2 - sqrt(2))/2: -8,
                (sqrt(5) - 1)/4: 10,
                (1 - sqrt(5))/4: -10,
                (sqrt(3) - 1)/sqrt(2**3): 12,
                (1 - sqrt(3))/sqrt(2**3): -12,
                (sqrt(5) + 1)/4: S(10)/3,
                -(sqrt(5) + 1)/4: -S(10)/3
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
                return p * (n - 2)**2/(n*(n - 1)) * x**2
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
        return -S.ImaginaryUnit*C.log(S.ImaginaryUnit*x + sqrt(1 - x**2))

    def _eval_is_real(self):
        return self.args[0].is_real and (self.args[0] >= -1 and self.args[0] <= 1)

    def inverse(self, argindex=1):
        """
        Returns the inverse of this function.
        """
        return sin

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

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -1/sqrt(1 - self.args[0]**2)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_is_rational(self):
        s = self.func(*self.args)
        if s.func == self.func:
            if s.args[0].is_rational:
                return False
        else:
            return s.is_rational

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
                S.Half: S.Pi/3,
                -S.Half: 2*S.Pi/3,
                sqrt(2)/2: S.Pi/4,
                -sqrt(2)/2: 3*S.Pi/4,
                1/sqrt(2): S.Pi/4,
                -1/sqrt(2): 3*S.Pi/4,
                sqrt(3)/2: S.Pi/6,
                -sqrt(3)/2: 5*S.Pi/6,
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
                return p * (n - 2)**2/(n*(n - 1)) * x**2
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
        return self.args[0].is_real and (self.args[0] >= -1 and self.args[0] <= 1)

    def _eval_rewrite_as_log(self, x):
        return S.Pi/2 + S.ImaginaryUnit * C.log(S.ImaginaryUnit * x + sqrt(1 - x**2))

    def _eval_rewrite_as_asin(self, x):
        return S.Pi/2 - asin(x)

    def _eval_rewrite_as_atan(self, x):
        return atan(sqrt(1 - x**2)/x) + (S.Pi/2)*(1 - x*sqrt(1/x**2))

    def inverse(self, argindex=1):
        """
        Returns the inverse of this function.
        """
        return cos

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

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/(1 + self.args[0]**2)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_is_rational(self):
        s = self.func(*self.args)
        if s.func == self.func:
            if s.args[0].is_rational:
                return False
        else:
            return s.is_rational

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
                sqrt(3)/3: 6,
                -sqrt(3)/3: -6,
                1/sqrt(3): 6,
                -1/sqrt(3): -6,
                sqrt(3): 3,
                -sqrt(3): -3,
                (1 + sqrt(2)): S(8)/3,
                -(1 + sqrt(2)): S(8)/3,
                (sqrt(2) - 1): 8,
                (1 - sqrt(2)): -8,
                sqrt((5 + 2*sqrt(5))): S(5)/2,
                -sqrt((5 + 2*sqrt(5))): -S(5)/2,
                (2 - sqrt(3)): 12,
                -(2 - sqrt(3)): -12
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
            return (-1)**((n - 1)//2) * x**n / n

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and C.Order(1, x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_is_real(self):
        return self.args[0].is_real

    def _eval_rewrite_as_log(self, x):
        return S.ImaginaryUnit/2 * (C.log(
            (S(1) - S.ImaginaryUnit * x)/(S(1) + S.ImaginaryUnit * x)))

    def _eval_aseries(self, n, args0, x, logx):
        if args0[0] == S.Infinity:
            return S.Pi/2 - atan(1/self.args[0])
        elif args0[0] == S.NegativeInfinity:
            return -S.Pi/2 - atan(1/self.args[0])
        else:
            return super(atan, self)._eval_aseries(n, args0, x, logx)

    def inverse(self, argindex=1):
        """
        Returns the inverse of this function.
        """
        return tan

    def _sage_(self):
        import sage.all as sage
        return sage.atan(self.args[0]._sage_())


class acot(Function):
    """
    acot(x) -> Returns the arc cotangent of x (measured in radians)
    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -1 / (1 + self.args[0]**2)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_is_rational(self):
        s = self.func(*self.args)
        if s.func == self.func:
            if s.args[0].is_rational:
                return False
        else:
            return s.is_rational

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
                sqrt(3)/3: 3,
                -sqrt(3)/3: -3,
                1/sqrt(3): 3,
                -1/sqrt(3): -3,
                sqrt(3): 6,
                -sqrt(3): -6,
                (1 + sqrt(2)): 8,
                -(1 + sqrt(2)): -8,
                (1 - sqrt(2)): -S(8)/3,
                (sqrt(2) - 1): S(8)/3,
                sqrt(5 + 2*sqrt(5)): 10,
                -sqrt(5 + 2*sqrt(5)): -10,
                (2 + sqrt(3)): 12,
                -(2 + sqrt(3)): -12,
                (2 - sqrt(3)): S(12)/5,
                -(2 - sqrt(3)): -S(12)/5,
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
            return S.Pi / 2  # FIX THIS
        elif n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)
            return (-1)**((n + 1)//2) * x**n / n

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

    def _eval_rewrite_as_log(self, x):
        return S.ImaginaryUnit/2 * \
            (C.log((x - S.ImaginaryUnit)/(x + S.ImaginaryUnit)))

    def inverse(self, argindex=1):
        """
        Returns the inverse of this function.
        """
        return cot

    def _sage_(self):
        import sage.all as sage
        return sage.acot(self.args[0]._sage_())


class atan2(Function):
    r"""
    The function ``atan2(y, x)`` computes `\operatorname{atan}(y/x)` taking
    two arguments `y` and `x`.  Signs of both `y` and `x` are considered to
    determine the appropriate quadrant of `\operatorname{atan}(y/x)`.
    The range is `(-\pi, \pi]`. The complete definition reads as follows:

    .. math::

        \operatorname{atan2}(y, x) =
        \begin{cases}
          \arctan\left(\frac y x\right) & \qquad x > 0 \\
          \arctan\left(\frac y x\right) + \pi& \qquad y \ge 0 , x < 0 \\
          \arctan\left(\frac y x\right) - \pi& \qquad y < 0 , x < 0 \\
          +\frac{\pi}{2} & \qquad y > 0 , x = 0 \\
          -\frac{\pi}{2} & \qquad y < 0 , x = 0 \\
          \text{undefined} & \qquad y = 0, x = 0
        \end{cases}

    Attention: Note the role reversal of both arguments. The `y`-coordinate
    is the first argument and the `x`-coordinate the second.

    Examples
    ========

    Going counter-clock wise around the origin we find the
    following angles:

    >>> from sympy import atan2
    >>> atan2(0, 1)
    0
    >>> atan2(1, 1)
    pi/4
    >>> atan2(1, 0)
    pi/2
    >>> atan2(1, -1)
    3*pi/4
    >>> atan2(0, -1)
    pi
    >>> atan2(-1, -1)
    -3*pi/4
    >>> atan2(-1, 0)
    -pi/2
    >>> atan2(-1, 1)
    -pi/4

    which are all correct. Compare this to the results of the ordinary
    `\operatorname{atan}` function for the point `(x, y) = (-1, 1)`

    >>> from sympy import atan, S
    >>> atan(S(1) / -1)
    -pi/4
    >>> atan2(1, -1)
    3*pi/4

    where only the `\operatorname{atan2}` function reurns what we expect.
    We can differentiate the function with respect to both arguments:

    >>> from sympy import diff
    >>> from sympy.abc import x, y
    >>> diff(atan2(y, x), x)
    -y/(x**2 + y**2)

    >>> diff(atan2(y, x), y)
    x/(x**2 + y**2)

    We can express the `\operatorname{atan2}` function in terms of
    complex logarithms:

    >>> from sympy import log
    >>> atan2(y, x).rewrite(log)
    -I*log((x + I*y)/sqrt(x**2 + y**2))

    and in terms of `\operatorname(atan)`:

    >>> from sympy import atan
    >>> atan2(y, x).rewrite(atan)
    2*atan(y/(x + sqrt(x**2 + y**2)))

    but note that this form is undefined on the negative real axis.

    See Also
    ========

    sin, cos, sec, csc, tan, cot
    asin, acos, atan

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Atan2
    .. [2] http://functions.wolfram.com/ElementaryFunctions/ArcTan2/
    """

    @classmethod
    def eval(cls, y, x):
        if x is S.NegativeInfinity:
            if y.is_zero:
                # Special case y = 0 because we define Heaviside(0) = 1/2
                return S.Pi
            return 2*S.Pi*(C.Heaviside(C.re(y))) - S.Pi
        elif x is S.Infinity:
            return S.Zero

        if x.is_real and y.is_real:
            if x.is_positive:
                return atan(y / x)
            elif x.is_negative:
                if y.is_negative:
                    return atan(y / x) - S.Pi
                else:
                    return atan(y / x) + S.Pi
            elif x.is_zero:
                if y.is_positive:
                    return S.Pi/2
                elif y.is_negative:
                    return -S.Pi/2
                elif y.is_zero:
                    return S.NaN

        if y.is_zero and x.is_real and x.is_nonzero:
            return S.Pi * (S.One - C.Heaviside(x))

    def _eval_rewrite_as_log(self, y, x):
        return -S.ImaginaryUnit*C.log((x + S.ImaginaryUnit*y) / sqrt(x**2 + y**2))

    def _eval_rewrite_as_atan(self, y, x):
        return 2*atan(y / (sqrt(x**2 + y**2) + x))

    def _eval_rewrite_as_arg(self, y, x):
        if (x.is_real or x.is_imaginary) and (y.is_real or y.is_imaginary):
            return C.arg(x + y*S.ImaginaryUnit)

    def _eval_is_real(self):
        return self.args[0].is_real and self.args[1].is_real

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate(), self.args[1].conjugate())

    def fdiff(self, argindex):
        y, x = self.args
        if argindex == 1:
            # Diff wrt y
            return x/(x**2 + y**2)
        elif argindex == 2:
            # Diff wrt x
            return -y/(x**2 + y**2)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_evalf(self, prec):
        y, x = self.args
        if x.is_real and y.is_real:
            super(atan2, self)._eval_evalf(prec)

    def _sage_(self):
        import sage.all as sage
        return sage.atan2(self.args[0]._sage_(), self.args[1]._sage_())

from sympy.core import S, C
from sympy.core.function import Function, Derivative, ArgumentIndexError
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.piecewise import Piecewise
from sympy.core import Add, Mul
from sympy.core.relational import Eq

###############################################################################
######################### REAL and IMAGINARY PARTS ############################
###############################################################################

class re(Function):
    """Returns real part of expression. This function performs only
       elementary analysis and so it will fail to decompose properly
       more complicated expressions. If completely simplified result
       is needed then use Basic.as_real_imag() or perform complex
       expansion on instance of this function.

       >>> from sympy import re, im, I, E
       >>> from sympy.abc import x, y

       >>> re(2*E)
       2*E

       >>> re(2*I + 17)
       17

       >>> re(2*I)
       0

       >>> re(im(x) + x*I + 2)
       2

       See Also
       ========

       im
    """
    nargs = 1

    is_real = True
    unbranched = True # implicitely works on the projection to C

    @classmethod
    def eval(cls, arg):
        if arg is S.NaN:
            return S.NaN
        elif arg.is_real:
            return arg
        elif arg.is_Function and arg.func == conjugate:
            return re(arg.args[0])
        else:

            included, reverted, excluded = [], [], []
            arg = Add.make_args(arg)
            for term in arg:
                coeff = term.as_coefficient(S.ImaginaryUnit)

                if coeff is not None:
                    if not coeff.is_real:
                        reverted.append(coeff)
                elif not term.has(S.ImaginaryUnit) and term.is_real:
                    excluded.append(term)
                else:
                    included.append(term)

            if len(arg) != len(included):
                a, b, c = map(lambda xs: Add(*xs),
                    [included, reverted, excluded])

                return cls(a) - im(b) + c

    def _eval_conjugate(self):
        return self

    def as_real_imag(self, deep=True):
        """
        Returns the real number with a zero complex part.
        """
        return (self, S.Zero)

    def _eval_expand_complex(self, deep=True, **hints):
#        if deep:
#            return self.args[0].expand(deep, **hints).as_real_imag()[0]
#        else:
        return self.args[0].as_real_imag()[0]

    def _eval_derivative(self, x):
        return re(Derivative(self.args[0], x, **{'evaluate': True}))

class im(Function):
    """
    Returns imaginary part of expression. This function performs only
    elementary analysis and so it will fail to decompose properly more
    complicated expressions. If completely simplified result is needed then
    use Basic.as_real_imag() or perform complex expansion on instance of
    this function.

    Examples
    ========

    >>> from sympy import re, im, E, I
    >>> from sympy.abc import x, y

    >>> im(2*E)
    0

    >>> re(2*I + 17)
    17

    >>> im(x*I)
    re(x)

    >>> im(re(x) + y)
    im(y)

    See Also
    ========

    re
    """

    nargs = 1

    is_real = True
    unbranched = True # implicitely works on the projection to C

    @classmethod
    def eval(cls, arg):
        if arg is S.NaN:
            return S.NaN
        elif arg.is_real:
            return S.Zero
        elif arg.is_Function and arg.func == conjugate:
            return -im(arg.args[0])
        else:
            included, reverted, excluded = [], [], []
            arg = Add.make_args(arg)
            for term in arg:
                coeff = term.as_coefficient(S.ImaginaryUnit)

                if coeff is not None:
                    if not coeff.is_real:
                        reverted.append(coeff)
                    else:
                        excluded.append(coeff)
                elif term.has(S.ImaginaryUnit) or not term.is_real:
                    included.append(term)

            if len(arg) != len(included):
                a, b, c = map(lambda xs: Add(*xs),
                    [included, reverted, excluded])

                return cls(a) + re(b) + c

    def _eval_conjugate(self):
        return self

    def as_real_imag(self, deep=True):
        """
        Return the imaginary part with a zero real part.

        Examples
        ========

        >>> from sympy.functions import im
        >>> from sympy import I
        >>> im(2 + 3*I).as_real_imag()
        (3, 0)
        """
        return (self, S.Zero)

    def _eval_expand_complex(self, deep=True, **hints):
#        if deep:
#            return self.args[0].expand(deep, **hints).as_real_imag()[1]
        return self.args[0].as_real_imag()[1]

    def _eval_derivative(self, x):
        return im(Derivative(self.args[0], x, **{'evaluate': True}))

###############################################################################
############### SIGN, ABSOLUTE VALUE, ARGUMENT and CONJUGATION ################
###############################################################################

class sign(Function):
    """
    Returns the sign of an expression, that is:

    * 1 if expression is positive
    * 0 if expression is equal to zero
    * -1 if expression is negative

    Examples
    ========

    >>> from sympy.functions import sign
    >>> sign(-1)
    -1
    >>> sign(0)
    0

    See Also
    ========

    Abs, conjugate
    """

    nargs = 1

    @classmethod
    def eval(cls, arg):
        if arg is S.NaN:
            return S.NaN
        if arg is S.Zero: return S.Zero
        if arg.is_positive: return S.One
        if arg.is_negative: return S.NegativeOne
        if arg.is_Function:
            if arg.func is sign: return arg
        if arg.is_Mul:
            c, args = arg.as_coeff_mul()
            unk = []
            is_neg = c.is_negative
            for ai in args:
                if ai.is_negative is None:
                    unk.append(ai)
                elif ai.is_negative:
                    is_neg = not is_neg
            if c is S.One and len(unk) == len(args):
                return None
            return (S.NegativeOne if is_neg else S.One) * cls(arg._new_rawargs(*unk))

    is_bounded = True

    def _eval_derivative(self, x):
        return S.Zero

    def _eval_conjugate(self):
        return self

    def _eval_is_zero(self):
        return (self.args[0] is S.Zero)

    def _sage_(self):
        import sage.all as sage
        return sage.sgn(self.args[0]._sage_())

class Abs(Function):
    """
    Return the absolute value of the argument.

    This is an extension of the built-in function abs() to accept symbolic
    values.  If you pass a SymPy expression to the built-in abs(), it will
    pass it automatically to Abs().

    Examples
    ========

    >>> from sympy import Abs, Symbol, S
    >>> Abs(-1)
    1
    >>> x = Symbol('x', real=True)
    >>> Abs(-x)
    Abs(x)
    >>> Abs(x**2)
    x**2
    >>> abs(-x) # The Python built-in
    Abs(x)

    Note that the Python built-in will return either an Expr or int depending on
    the argument::

        >>> type(abs(-1))
        <... 'int'>
        >>> type(abs(S.NegativeOne))
        <class 'sympy.core.numbers.One'>

    Abs will always return a sympy object.

    See Also
    ========

    sign, conjugate
    """

    nargs = 1

    is_real = True
    is_negative = False

    def fdiff(self, argindex=1):
        """
        Get the first derivative of the argument to Abs().

        Examples
        ========

        >>> from sympy.abc import x
        >>> from sympy.functions import Abs
        >>> Abs(-x).fdiff()
        sign(x)
        """
        if argindex == 1:
            return sign(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        if arg is S.NaN:
            return S.NaN
        if arg.is_zero:     return arg
        if arg.is_positive: return arg
        if arg.is_negative: return -arg
        coeff, terms = arg.as_coeff_mul()
        if coeff is not S.One:
            return cls(coeff) * cls(Mul(*terms))
        if arg.is_real is False:
            from sympy import expand_mul
            return sqrt( expand_mul(arg * arg.conjugate()) )
        if arg.is_Pow:
            base, exponent = arg.as_base_exp()
            if exponent.is_even and base.is_real:
                return arg
        return

    def _eval_is_nonzero(self):
        return self._args[0].is_nonzero

    def _eval_is_positive(self):
        return self.is_nonzero

    def _eval_conjugate(self):
        return self

    def _eval_power(self,other):
        if self.args[0].is_real and other.is_integer:
            if other.is_even:
                return self.args[0]**other
            elif other is not S.NegativeOne and other.is_Integer:
                e = other - sign(other)
                return self.args[0]**e*self
        return

    def _eval_nseries(self, x, n, logx):
        direction = self.args[0].leadterm(x)[0]
        s = self.args[0]._eval_nseries(x, n=n, logx=logx)
        when = Eq(direction, 0)
        return Piecewise(
                         ((s.subs(direction, 0)), when),
                         (sign(direction)*s, True),
                         )

    def _sage_(self):
        import sage.all as sage
        return sage.abs_symbolic(self.args[0]._sage_())

    def _eval_derivative(self, x):
        if self.args[0].is_real:
            return Derivative(self.args[0], x, **{'evaluate': True}) * sign(self.args[0])
        return (re(self.args[0]) * re(Derivative(self.args[0], x,
            **{'evaluate': True})) + im(self.args[0]) * im(Derivative(self.args[0],
                x, **{'evaluate': True}))) / Abs(self.args[0])

class arg(Function):
    """Returns the argument (in radians) of a complex number"""

    nargs = 1

    is_real = True
    is_bounded = True

    @classmethod
    def eval(cls, arg):
        x, y = re(arg), im(arg)
        arg = C.atan2(y, x)
        if arg.is_number:
            return arg

    def _eval_conjugate(self):
        return self

    def _eval_derivative(self, t):
        x, y = re(self.args[0]), im(self.args[0])
        return (x * Derivative(y, t, **{'evaluate': True}) - y *
                Derivative(x, t, **{'evaluate': True})) / (x**2 + y**2)

class conjugate(Function):
    """
    Changes the sign of the imaginary part of a complex number.

    Examples
    ========

    >>> from sympy import conjugate, I

    >>> conjugate(1 + I)
    1 - I

    See Also
    ========

    sign, Abs
    """

    nargs = 1

    @classmethod
    def eval(cls, arg):
        obj = arg._eval_conjugate()
        if obj is not None:
            return obj

    def _eval_conjugate(self):
        return self.args[0]

    def _eval_derivative(self, x):
        return conjugate(Derivative(self.args[0], x, **{'evaluate': True}))

###############################################################################
############### HANDLING OF POLAR NUMBERS #####################################
###############################################################################

class polar_lift(Function):
    """
    Lift argument to the riemann surface of the logarithm, using the
    standard branch.

    >>> from sympy import Symbol, polar_lift, I
    >>> p = Symbol('p', polar=True)
    >>> x = Symbol('x')
    >>> polar_lift(4)
    4*exp_polar(0)
    >>> polar_lift(-4)
    4*exp_polar(I*pi)
    >>> polar_lift(-I)
    exp_polar(-I*pi/2)
    >>> polar_lift(I + 2)
    polar_lift(2 + I)

    >>> polar_lift(4*x)
    4*polar_lift(x)
    >>> polar_lift(4*p)
    4*p

    See Also
    ========

    sympy.functions.elementary.exponential.exp_polar
    periodic_argument
    """

    nargs = 1

    is_polar = True
    is_comparable = False # Cannot be evalf'd.

    @classmethod
    def eval(cls, arg):
        from sympy import exp_polar, pi, I, arg as argument
        if arg.is_number:
            ar = argument(arg)
            #if not ar.has(argument) and not ar.has(atan):
            if ar in (0, pi/2, -pi/2, pi):
                return exp_polar(I*ar)*abs(arg)

        if arg.is_Mul:
            args = arg.args
        else:
            args = [arg]
        included = []
        excluded = []
        positive = []
        for arg in args:
            if arg.is_polar:
                included += [arg]
            elif arg.is_positive:
                positive += [arg]
            else:
                excluded += [arg]
        if len(excluded) < len(args):
            if excluded:
                return Mul(*(included + positive))*polar_lift(Mul(*excluded))
            elif included:
                return Mul(*(included + positive))
            else:
                return Mul(*positive)*exp_polar(0)

    def _eval_evalf(self, prec):
        """ Careful! any evalf of polar numbers is flaky """
        return self.args[0]._eval_evalf(prec)

class periodic_argument(Function):
    """
    Represent the argument on a quotient of the riemann surface of the
    logarithm. That is, given a period P, always return a value in
    (-P/2, P/2], by using exp(P*I) == 1.

    >>> from sympy import exp, exp_polar, periodic_argument, unbranched_argument
    >>> from sympy import I, pi
    >>> unbranched_argument(exp(5*I*pi))
    pi
    >>> unbranched_argument(exp_polar(5*I*pi))
    5*pi
    >>> periodic_argument(exp_polar(5*I*pi), 2*pi)
    pi
    >>> periodic_argument(exp_polar(5*I*pi), 3*pi)
    -pi
    >>> periodic_argument(exp_polar(5*I*pi), pi)
    0

    See Also
    ========

    sympy.functions.elementary.exponential.exp_polar
    polar_lift : Lift argument to the riemann surface of the logarithm
    principal_branch
    """

    nargs = 2

    @classmethod
    def _getunbranched(cls, ar):
        from sympy import exp_polar, log
        if ar.is_Mul:
            args = ar.args
        else:
            args = [ar]
        unbranched = 0
        for a in args:
            if not a.is_polar:
                unbranched += arg(a)
            elif a.func is exp_polar:
                unbranched += a.args[0].as_real_imag()[1]
            elif a.is_Pow:
                re, im = a.exp.as_real_imag()
                unbranched += re*unbranched_argument(a.base) + im*log(abs(a.base))
            else:
                return None
        return unbranched

    @classmethod
    def eval(cls, ar, period):
        # Our strategy is to evaluate the argument on the riemann surface of the
        # logarithm, and then reduce.
        # NOTE evidently this means it is a rather bad idea to use this with
        # period != 2*pi and non-polar numbers.
        from sympy import ceiling, oo, atan2, atan, polar_lift, pi
        if not period.is_positive:
            return None
        if period == oo and isinstance(ar, principal_branch):
            return periodic_argument(*ar.args)
        if ar.func is polar_lift and period >= 2*pi:
            return periodic_argument(ar.args[0], period)
        unbranched = cls._getunbranched(ar)
        if unbranched is None:
            return None
        if unbranched.has(periodic_argument, atan2, arg, atan):
            return None
        if period == oo:
            return unbranched
        if period != oo:
            n = ceiling(unbranched/period - S(1)/2)*period
            if not n.has(ceiling):
                return unbranched - n

    def _eval_evalf(self, prec):
        from sympy import ceiling, oo
        z, period = self.args
        if period == oo:
            unbranched = periodic_argument._getunbranched(z)
            if unbranched is None:
                return self
            return unbranched._eval_evalf(prec)
        ub = periodic_argument(z, oo)._eval_evalf(prec)
        return (ub - ceiling(ub/period - S(1)/2)*period)._eval_evalf(prec)

def unbranched_argument(arg):
    from sympy import oo
    return periodic_argument(arg, oo)

class principal_branch(Function):
    """
    Represent a polar number reduced to its principal branch on a quotient
    of the riemann surface of the logarithm.

    This is a function of two arguments. The first argument is a polar
    number `z`, and the second one a positive real number of infinity, `p`.
    The result is "z mod exp_polar(I*p)".

    >>> from sympy import exp_polar, principal_branch, oo, I, pi
    >>> from sympy.abc import z
    >>> principal_branch(z, oo)
    z
    >>> principal_branch(exp_polar(2*pi*I)*3, 2*pi)
    3*exp_polar(0)
    >>> principal_branch(exp_polar(2*pi*I)*3*z, 2*pi)
    3*principal_branch(z, 2*pi)

    See Also
    ========

    sympy.functions.elementary.exponential.exp_polar
    polar_lift : Lift argument to the riemann surface of the logarithm
    periodic_argument
    """

    nargs = 2
    is_polar = True
    is_comparable = False # cannot always be evalf'd

    @classmethod
    def eval(self, x, period):
        from sympy import oo, exp_polar, I, Mul, polar_lift, Symbol
        if isinstance(x, polar_lift):
            return principal_branch(x.args[0], period)
        if period == oo:
            return x
        ub = periodic_argument(x, oo)
        barg = periodic_argument(x, period)
        if ub != barg and not ub.has(periodic_argument) \
           and not barg.has(periodic_argument):
            pl = polar_lift(x)
            def mr(expr):
                if not isinstance(expr, Symbol):
                    return polar_lift(expr)
                return expr
            pl = pl.replace(polar_lift, mr)
            if not pl.has(polar_lift):
                res = exp_polar(I*(barg - ub))*pl
                if not res.is_polar and not res.has(exp_polar):
                    res *= exp_polar(0)
                return res

        if not x.free_symbols:
            c, m = x, ()
        else:
            c, m = x.as_coeff_mul(*x.free_symbols)
        others = []
        for y in m:
            if y.is_positive:
                c *= y
            else:
                others += [y]
        m = tuple(others)
        arg = periodic_argument(c, period)
        if arg.has(periodic_argument):
            return None
        if arg.is_number and (unbranched_argument(c) != arg or \
                              (arg == 0 and m != () and c != 1)):
            if arg == 0:
                return abs(c)*principal_branch(Mul(*m), period)
            return principal_branch(exp_polar(I*arg)*Mul(*m), period)*abs(c)
        if arg.is_number and ((abs(arg) < period/2) is True or arg == period/2) \
           and m == ():
            return exp_polar(arg*I)*abs(c)

    def _eval_evalf(self, prec):
        from sympy import exp, pi, I
        z, period = self.args
        p = periodic_argument(z, period)._eval_evalf(prec)
        if abs(p) > pi or p == -pi:
            return self # Cannot evalf for this argument.
        return (abs(z)*exp(I*p))._eval_evalf(prec)

# /cyclic/
from sympy.core import basic as _
_.abs_ = Abs
del _

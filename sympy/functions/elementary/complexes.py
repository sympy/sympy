from __future__ import print_function, division

from sympy.core import S, C
from sympy.core.compatibility import u
from sympy.core.exprtools import factor_terms
from sympy.core.function import (Function, Derivative, ArgumentIndexError,
    AppliedUndef)
from sympy.core.logic import fuzzy_not
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.piecewise import Piecewise
from sympy.core import Add, Mul
from sympy.core.relational import Eq
from sympy.functions.elementary.trigonometric import atan, atan2

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

    is_real = True
    unbranched = True  # implicitely works on the projection to C

    @classmethod
    def eval(cls, arg):
        if arg is S.NaN:
            return S.NaN
        elif arg.is_real:
            return arg
        elif arg.is_imaginary or (S.ImaginaryUnit*arg).is_real:
            return S.Zero
        elif arg.is_Function and arg.func is conjugate:
            return re(arg.args[0])
        else:

            included, reverted, excluded = [], [], []
            args = Add.make_args(arg)
            for term in args:
                coeff = term.as_coefficient(S.ImaginaryUnit)

                if coeff is not None:
                    if not coeff.is_real:
                        reverted.append(coeff)
                elif not term.has(S.ImaginaryUnit) and term.is_real:
                    excluded.append(term)
                else:
                    # Try to do some advanced expansion.  If
                    # impossible, don't try to do re(arg) again
                    # (because this is what we are trying to do now).
                    real_imag = term.as_real_imag(ignore=arg)
                    if real_imag:
                        excluded.append(real_imag[0])
                    else:
                        included.append(term)

            if len(args) != len(included):
                a, b, c = map(lambda xs: Add(*xs),
                    [included, reverted, excluded])

                return cls(a) - im(b) + c

    def as_real_imag(self, deep=True, **hints):
        """
        Returns the real number with a zero complex part.
        """
        return (self, S.Zero)

    def _eval_derivative(self, x):
        if x.is_real or self.args[0].is_real:
            return re(Derivative(self.args[0], x, evaluate=True))
        if x.is_imaginary or self.args[0].is_imaginary:
            return -S.ImaginaryUnit \
                * im(Derivative(self.args[0], x, evaluate=True))

    def _sage_(self):
        import sage.all as sage
        return sage.real_part(self.args[0]._sage_())


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

    is_real = True
    unbranched = True  # implicitely works on the projection to C

    @classmethod
    def eval(cls, arg):
        if arg is S.NaN:
            return S.NaN
        elif arg.is_real:
            return S.Zero
        elif arg.is_imaginary or (S.ImaginaryUnit*arg).is_real:
            return -S.ImaginaryUnit * arg
        elif arg.is_Function and arg.func is conjugate:
            return -im(arg.args[0])
        else:
            included, reverted, excluded = [], [], []
            args = Add.make_args(arg)
            for term in args:
                coeff = term.as_coefficient(S.ImaginaryUnit)

                if coeff is not None:
                    if not coeff.is_real:
                        reverted.append(coeff)
                    else:
                        excluded.append(coeff)
                elif term.has(S.ImaginaryUnit) or not term.is_real:
                    # Try to do some advanced expansion.  If
                    # impossible, don't try to do im(arg) again
                    # (because this is what we are trying to do now).
                    real_imag = term.as_real_imag(ignore=arg)
                    if real_imag:
                        excluded.append(real_imag[1])
                    else:
                        included.append(term)

            if len(args) != len(included):
                a, b, c = map(lambda xs: Add(*xs),
                    [included, reverted, excluded])

                return cls(a) + re(b) + c

    def as_real_imag(self, deep=True, **hints):
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

    def _eval_derivative(self, x):
        if x.is_real or self.args[0].is_real:
            return im(Derivative(self.args[0], x, evaluate=True))
        if x.is_imaginary or self.args[0].is_imaginary:
            return -S.ImaginaryUnit \
                * re(Derivative(self.args[0], x, evaluate=True))

    def _sage_(self):
        import sage.all as sage
        return sage.imag_part(self.args[0]._sage_())


###############################################################################
############### SIGN, ABSOLUTE VALUE, ARGUMENT and CONJUGATION ################
###############################################################################

class sign(Function):
    """
    Returns the complex sign of an expression:

    If the expresssion is real the sign will be:

        * 1 if expression is positive
        * 0 if expression is equal to zero
        * -1 if expression is negative

    If the expresssion is imaginary the sign will be:

        * I if im(expression) is positive
        * -I if im(expression) is negative

    Otherwise an unevaluated expression will be returned. When evaluated, the
    result (in general) will be ``cos(arg(expr)) + I*sin(arg(expr))``.

    Examples
    ========

    >>> from sympy.functions import sign
    >>> from sympy.core.numbers import I

    >>> sign(-1)
    -1
    >>> sign(0)
    0
    >>> sign(-3*I)
    -I
    >>> sign(1 + I)
    sign(1 + I)
    >>> _.evalf()
    0.707106781186548 + 0.707106781186548*I

    See Also
    ========

    Abs, conjugate
    """

    is_finite = True
    is_complex = True

    def doit(self):
        if self.args[0].is_nonzero:
            return self.args[0] / Abs(self.args[0])
        return self

    @classmethod
    def eval(cls, arg):
        # handle what we can
        if arg.is_Mul:
            c, args = arg.as_coeff_mul()
            unk = []
            s = sign(c)
            for a in args:
                if a.is_negative:
                    s = -s
                elif a.is_positive:
                    pass
                else:
                    ai = im(a)
                    if a.is_imaginary and ai.is_comparable:  # i.e. a = I*real
                        s *= S.ImaginaryUnit
                        if ai.is_negative:
                            # can't use sign(ai) here since ai might not be
                            # a Number
                            s = -s
                    else:
                        unk.append(a)
            if c is S.One and len(unk) == len(args):
                return None
            return s * cls(arg._new_rawargs(*unk))
        if arg is S.NaN:
            return S.NaN
        if arg.is_zero:  # it may be an Expr that is zero
            return S.Zero
        if arg.is_positive:
            return S.One
        if arg.is_negative:
            return S.NegativeOne
        if arg.is_Function:
            if arg.func is sign:
                return arg
        if arg.is_imaginary:
            if arg.is_Pow and arg.exp is S.Half:
                # we catch this because non-trivial sqrt args are not expanded
                # e.g. sqrt(1-sqrt(2)) --x-->  to I*sqrt(sqrt(2) - 1)
                return S.ImaginaryUnit
            arg2 = -S.ImaginaryUnit * arg
            if arg2.is_positive:
                return S.ImaginaryUnit
            if arg2.is_negative:
                return -S.ImaginaryUnit

    def _eval_Abs(self):
        if self.args[0].is_nonzero:
            return S.One

    def _eval_conjugate(self):
        return sign(conjugate(self.args[0]))

    def _eval_derivative(self, x):
        if self.args[0].is_real:
            from sympy.functions.special.delta_functions import DiracDelta
            return 2 * Derivative(self.args[0], x, evaluate=True) \
                * DiracDelta(self.args[0])
        elif self.args[0].is_imaginary:
            from sympy.functions.special.delta_functions import DiracDelta
            return 2 * Derivative(self.args[0], x, evaluate=True) \
                * DiracDelta(-S.ImaginaryUnit * self.args[0])

    def _eval_is_nonnegative(self):
        if self.args[0].is_nonnegative:
            return True

    def _eval_is_nonpositive(self):
        if self.args[0].is_nonpositive:
            return True

    def _eval_is_imaginary(self):
        return self.args[0].is_imaginary

    def _eval_is_integer(self):
        return self.args[0].is_real

    def _eval_is_zero(self):
        return self.args[0].is_zero

    def _eval_power(self, other):
        if (
            self.args[0].is_real and
            self.args[0].is_nonzero and
            other.is_integer and
            other.is_even
        ):
            return S.One

    def _sage_(self):
        import sage.all as sage
        return sage.sgn(self.args[0]._sage_())

    def _eval_rewrite_as_Piecewise(self, arg):
        if arg.is_real:
            return Piecewise((1, arg > 0), (-1, arg < 0), (0, True))

    def _eval_rewrite_as_Heaviside(self, arg):
        if arg.is_real:
            return C.Heaviside(arg)*2-1

    def _eval_simplify(self, ratio, measure):
        return self.func(self.args[0].factor())


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

    is_real = True
    is_negative = False
    unbranched = True

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
        from sympy.simplify.simplify import signsimp
        if hasattr(arg, '_eval_Abs'):
            obj = arg._eval_Abs()
            if obj is not None:
                return obj
        # handle what we can
        arg = signsimp(arg, evaluate=False)
        if arg.is_Mul:
            known = []
            unk = []
            for t in arg.args:
                tnew = cls(t)
                if tnew.func is cls:
                    unk.append(tnew.args[0])
                else:
                    known.append(tnew)
            known = Mul(*known)
            unk = cls(Mul(*unk), evaluate=False) if unk else S.One
            return known*unk
        if arg is S.NaN:
            return S.NaN
        if arg.is_zero:  # it may be an Expr that is zero
            return S.Zero
        if arg.is_nonnegative:
            return arg
        if arg.is_nonpositive:
            return -arg
        if arg.is_imaginary:
            arg2 = -S.ImaginaryUnit * arg
            if arg2.is_nonnegative:
                return arg2
        if arg.is_real is False and arg.is_imaginary is False:
            from sympy import expand_mul
            return sqrt( expand_mul(arg * arg.conjugate()) )
        if arg.is_real is None and arg.is_imaginary is None and arg.is_Add:
            if all(a.is_real or a.is_imaginary or (S.ImaginaryUnit*a).is_real for a in arg.args):
                from sympy import expand_mul
                return sqrt(expand_mul(arg * arg.conjugate()))
        if arg.is_Pow:
            base, exponent = arg.as_base_exp()
            if exponent.is_even and base.is_real:
                return arg
            if exponent.is_integer and base is S.NegativeOne:
                return S.One

    def _eval_is_nonzero(self):
        return self._args[0].is_nonzero

    def _eval_is_positive(self):
        return self.is_nonzero

    def _eval_power(self, other):
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
        if self.args[0].is_real or self.args[0].is_imaginary:
            return Derivative(self.args[0], x, evaluate=True) \
                * sign(conjugate(self.args[0]))
        return (re(self.args[0]) * Derivative(re(self.args[0]), x,
            evaluate=True) + im(self.args[0]) * Derivative(im(self.args[0]),
                x, evaluate=True)) / Abs(self.args[0])

    def _eval_rewrite_as_Heaviside(self, arg):
        # Note this only holds for real arg (since Heaviside is not defined
        # for complex arguments).
        if arg.is_real:
            return arg*(C.Heaviside(arg) - C.Heaviside(-arg))

    def _eval_rewrite_as_Piecewise(self, arg):
        if arg.is_real:
            return Piecewise((arg, arg >= 0), (-arg, True))

    def _eval_rewrite_as_sign(self, arg):
        return arg/C.sign(arg)


class arg(Function):
    """Returns the argument (in radians) of a complex number"""

    is_real = True
    is_finite = True

    @classmethod
    def eval(cls, arg):
        if not arg.is_Atom:
            c, arg_ = factor_terms(arg).as_coeff_Mul()
            if arg_.is_Mul:
                arg_ = Mul(*[a if (sign(a) not in (-1, 1)) else
                    sign(a) for a in arg_.args])
            arg_ = sign(c)*arg_
        else:
            arg_ = arg
        x, y = re(arg_), im(arg_)
        rv = C.atan2(y, x)
        if rv.is_number and not rv.atoms(AppliedUndef):
            return rv
        if arg_ != arg:
            return cls(arg_, evaluate=False)

    def _eval_derivative(self, t):
        x, y = re(self.args[0]), im(self.args[0])
        return (x * Derivative(y, t, evaluate=True) - y *
                    Derivative(x, t, evaluate=True)) / (x**2 + y**2)

    def _eval_rewrite_as_atan2(self, arg):
        x, y = re(self.args[0]), im(self.args[0])
        return atan2(y, x)

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

    @classmethod
    def eval(cls, arg):
        obj = arg._eval_conjugate()
        if obj is not None:
            return obj

    def _eval_Abs(self):
        return Abs(self.args[0], evaluate=True)

    def _eval_adjoint(self):
        return transpose(self.args[0])

    def _eval_conjugate(self):
        return self.args[0]

    def _eval_derivative(self, x):
        if x.is_real:
            return conjugate(Derivative(self.args[0], x, evaluate=True))
        elif x.is_imaginary:
            return -conjugate(Derivative(self.args[0], x, evaluate=True))

    def _eval_transpose(self):
        return adjoint(self.args[0])


class transpose(Function):
    """
    Linear map transposition.
    """

    @classmethod
    def eval(cls, arg):
        obj = arg._eval_transpose()
        if obj is not None:
            return obj

    def _eval_adjoint(self):
        return conjugate(self.args[0])

    def _eval_conjugate(self):
        return adjoint(self.args[0])

    def _eval_transpose(self):
        return self.args[0]


class adjoint(Function):
    """
    Conjugate transpose or Hermite conjugation.
    """

    @classmethod
    def eval(cls, arg):
        obj = arg._eval_adjoint()
        if obj is not None:
            return obj
        obj = arg._eval_transpose()
        if obj is not None:
            return conjugate(obj)

    def _eval_adjoint(self):
        return self.args[0]

    def _eval_conjugate(self):
        return transpose(self.args[0])

    def _eval_transpose(self):
        return conjugate(self.args[0])

    def _latex(self, printer, exp=None, *args):
        arg = printer._print(self.args[0])
        tex = r'%s^{\dag}' % arg
        if exp:
            tex = r'\left(%s\right)^{%s}' % (tex, printer._print(exp))
        return tex

    def _pretty(self, printer, *args):
        from sympy.printing.pretty.stringpict import prettyForm
        pform = printer._print(self.args[0], *args)
        if printer._use_unicode:
            pform = pform**prettyForm(u('\u2020'))
        else:
            pform = pform**prettyForm('+')
        return pform

###############################################################################
############### HANDLING OF POLAR NUMBERS #####################################
###############################################################################


class polar_lift(Function):
    """
    Lift argument to the Riemann surface of the logarithm, using the
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

    is_polar = True
    is_comparable = False  # Cannot be evalf'd.

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

    def _eval_Abs(self):
        return Abs(self.args[0], evaluate=True)


class periodic_argument(Function):
    """
    Represent the argument on a quotient of the Riemann surface of the
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
    polar_lift : Lift argument to the Riemann surface of the logarithm
    principal_branch
    """

    @classmethod
    def _getunbranched(cls, ar):
        from sympy import exp_polar, log, polar_lift
        if ar.is_Mul:
            args = ar.args
        else:
            args = [ar]
        unbranched = 0
        for a in args:
            if not a.is_polar:
                unbranched += arg(a)
            elif a.func is exp_polar:
                unbranched += a.exp.as_real_imag()[1]
            elif a.is_Pow:
                re, im = a.exp.as_real_imag()
                unbranched += re*unbranched_argument(
                    a.base) + im*log(abs(a.base))
            elif a.func is polar_lift:
                unbranched += arg(a.args[0])
            else:
                return None
        return unbranched

    @classmethod
    def eval(cls, ar, period):
        # Our strategy is to evaluate the argument on the Riemann surface of the
        # logarithm, and then reduce.
        # NOTE evidently this means it is a rather bad idea to use this with
        # period != 2*pi and non-polar numbers.
        from sympy import ceiling, oo, atan2, atan, polar_lift, pi, Mul
        if not period.is_positive:
            return None
        if period == oo and isinstance(ar, principal_branch):
            return periodic_argument(*ar.args)
        if ar.func is polar_lift and period >= 2*pi:
            return periodic_argument(ar.args[0], period)
        if ar.is_Mul:
            newargs = [x for x in ar.args if not x.is_positive]
            if len(newargs) != len(ar.args):
                return periodic_argument(Mul(*newargs), period)
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
    of the Riemann surface of the logarithm.

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
    polar_lift : Lift argument to the Riemann surface of the logarithm
    periodic_argument
    """

    is_polar = True
    is_comparable = False  # cannot always be evalf'd

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
        if arg.is_number and (unbranched_argument(c) != arg or
                              (arg == 0 and m != () and c != 1)):
            if arg == 0:
                return abs(c)*principal_branch(Mul(*m), period)
            return principal_branch(exp_polar(I*arg)*Mul(*m), period)*abs(c)
        if arg.is_number and ((abs(arg) < period/2) == True or arg == period/2) \
                and m == ():
            return exp_polar(arg*I)*abs(c)

    def _eval_evalf(self, prec):
        from sympy import exp, pi, I
        z, period = self.args
        p = periodic_argument(z, period)._eval_evalf(prec)
        if abs(p) > pi or p == -pi:
            return self  # Cannot evalf for this argument.
        return (abs(z)*exp(I*p))._eval_evalf(prec)

# /cyclic/
from sympy.core import basic as _
_.abs_ = Abs
del _

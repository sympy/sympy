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

    """
    nargs = 1

    is_real = True

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
        return (self, S.Zero)

    def _eval_expand_complex(self, deep=True, **hints):
#        if deep:
#            return self.args[0].expand(deep, **hints).as_real_imag()[0]
#        else:
        return self.args[0].as_real_imag()[0]

    def _eval_derivative(self, x):
        return re(Derivative(self.args[0], x, **{'evaluate': True}))

class im(Function):
    """Returns imaginary part of expression. This function performs
       only elementary analysis and so it will fail to decompose
       properly more complicated expressions. If completely simplified
       result is needed then use Basic.as_real_imag() or perform complex
       expansion on instance of this function.

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

    """

    nargs = 1

    is_real = True

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
    """Return the sign of an expression, that is:
        -1 if expr <  0
         0 if expr == 0
         1 if expr >  0
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
    """Return the absolute value of the argument.

    This is an extension of the built-in function abs() to accept symbolic
    values.  If you pass a SymPy expression to the built-in abs(), it will
    pass it automatically to Abs().

    Examples

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

    """

    nargs = 1

    is_real = True
    is_negative = False

    def fdiff(self, argindex=1):
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
            return sqrt( (arg * arg.conjugate()).expand() )
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
    """Changes the sign of the imaginary part of a complex number.

        >>> from sympy import conjugate, I

        >>> conjugate(1 + I)
        1 - I

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

# /cyclic/
from sympy.core import basic as _
_.abs_ = Abs
del _

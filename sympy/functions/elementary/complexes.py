
from sympy.core.basic import Basic, S, C, sympify
from sympy.core.function import Function
from sympy.functions.elementary.miscellaneous import sqrt

from sympy.utilities.decorator import deprecated

###############################################################################
######################### REAL and IMAGINARY PARTS ############################
###############################################################################

class re(Function):
    """Returns real part of expression. This function performs only
       elementary analysis and so it will fail to decompose properly
       more complicated expressions. If completely simplified result
       is needed then use Basic.as_real_imag() or perform complex
       expansion on instance of this function.

       >>> from sympy import *

       >>> x, y = symbols('x', 'y')

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
    def _eval_apply_subs(self, *args):
        return

    @classmethod
    @deprecated
    def canonize(cls, arg):
        return cls.eval(arg)

    @classmethod
    def eval(cls, arg):
        if arg is S.NaN:
            return S.NaN
        elif arg.is_real:
            return arg
        else:
            if not arg.is_Add:
                arg = [arg]

            included, reverted, excluded = [], [], []

            if isinstance(arg, Basic):
                arg = arg.args
            for term in arg:
                coeff = term.as_coefficient(S.ImaginaryUnit)

                if coeff is not None:
                    if not coeff.is_real:
                        reverted.append(coeff)
                elif not term.has(S.ImaginaryUnit) and term.is_real:
                    excluded.append(term)
                else:
                    included.append(term)

            if len(arg[:]) != len(included):
                a, b, c = map(lambda xs: C.Add(*xs),
                    [included, reverted, excluded])

                return cls(a) - im(b) + c

    def _eval_conjugate(self):
        return self

    def _eval_expand_complex(self, deep=True, **hints):
#        if deep:
#            return self.args[0].expand(deep, **hints).as_real_imag()[0]
#        else:
        return self.args[0].as_real_imag()[0]

class im(Function):
    """Returns imaginary part of expression. This function performs
       only elementary analysis and so it will fail to decompose
       properly more complicated expressions. If completely simplified
       result is needed then use Basic.as_real_imag() or perform complex
       expansion on instance of this function.

       >>> from sympy import *

       >>> x, y = symbols('x', 'y')

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
    def _eval_apply_subs(self, *args):
        return

    @classmethod
    @deprecated
    def canonize(cls, arg):
        return cls.eval(arg)

    @classmethod
    def eval(cls, arg):
        if arg is S.NaN:
            return S.NaN
        elif arg.is_real:
            return S.Zero
        else:
            if not arg.is_Add:
                arg = [arg]

            included, reverted, excluded = [], [], []

            if isinstance(arg, Basic):
                arg = arg.args
            for term in arg:
                coeff = term.as_coefficient(S.ImaginaryUnit)

                if coeff is not None:
                    if not coeff.is_real:
                        reverted.append(coeff)
                    else:
                        excluded.append(coeff)
                elif term.has(S.ImaginaryUnit) or not term.is_real:
                    included.append(term)

            if len(arg[:]) != len(included):
                a, b, c = map(lambda xs: C.Add(*xs),
                    [included, reverted, excluded])

                return cls(a) + re(b) + c

    def _eval_conjugate(self):
        return self

    def _eval_expand_complex(self, deep=True, **hints):
#        if deep:
#            return self.args[0].expand(deep, **hints).as_real_imag()[1]
        return self.args[0].as_real_imag()[1]

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
    @deprecated
    def canonize(cls, arg):
        return cls.eval(arg)

    @classmethod
    def eval(cls, arg):
        if arg is S.NaN:
            return S.NaN
        if arg is S.Zero: return S.Zero
        if arg.is_positive: return S.One
        if arg.is_negative: return S.NegativeOne
        if arg.is_Mul:
            coeff, terms = arg.as_coeff_terms()
            if coeff is not S.One:
                return cls(coeff) * cls(C.Mul(*terms))

    is_bounded = True

    def _eval_conjugate(self):
        return self

    def _eval_is_zero(self):
        return (self[0] is S.Zero)

class abs(Function):
    """Return the absolute value of the argument. This is an extension of the built-in
    function abs to accept symbolic values

    Examples

        >>> from sympy import abs, Symbol
        >>> abs(-1)
        1
        >>> x = Symbol('x', real=True)
        >>> abs(-x)
        abs(x)
        >>> abs(x**2)
        x**2
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
    def _eval_apply_subs(self, *args):
        return

    @classmethod
    @deprecated
    def canonize(cls, arg):
        return cls.eval(arg)

    @classmethod
    def eval(cls, arg):
        if arg is S.NaN:
            return S.NaN
        if arg.is_zero:     return arg
        if arg.is_positive: return arg
        if arg.is_negative: return -arg
        coeff, terms = arg.as_coeff_terms()
        if coeff is not S.One:
            return cls(coeff) * cls(C.Mul(*terms))
        if arg.is_real is False:
            return sqrt( (arg * arg.conjugate()).expand() )
        if arg.is_Pow:
            base, exponent = arg.as_base_exp()
            if exponent.is_even and base.is_real:
                return arg
        return

    @classmethod
    def _eval_apply_evalf(cls, arg):
        # XXX this is weird!!!
        # XXX remove me when 'abs -> abs_' is done
        arg = arg.evalf()

        if arg.is_Number:
            import operator
            return operator.abs(float(arg))

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

    def nseries(self, x, x0, n):
        direction = self.args[0].leadterm(x)[0]
        return sign(direction)*self.args[0].nseries(x, x0, n)

    def _sage_(self):
        import sage.all as sage
        return sage.abs_symbolic(self.args[0]._sage_())

class arg(Function):
    """Returns the argument (in radians) of a complex number"""

    nargs = 1

    is_real = True
    is_bounded = True

    @classmethod
    @deprecated
    def canonize(cls, arg):
        return cls.eval(arg)

    @classmethod
    def eval(cls, arg):
        x, y = re(arg), im(arg)
        arg = C.atan2(y, x)
        if arg.is_number:
            return arg

    def _eval_conjugate(self):
        return self

class conjugate(Function):
    """Changes the sign of the imaginary part of a complex number.

        >>> from sympy import *

        >>> conjugate(1 + I)
        1 - I

    """

    nargs = 1

    @classmethod
    @deprecated
    def canonize(cls, arg):
        return cls.eval(arg)

    @classmethod
    def eval(cls, arg):
        obj = arg._eval_conjugate()
        if obj is not None:
            return obj

    def _eval_conjugate(self):
        return self.args[0]


# /cyclic/
from sympy.core import basic as _
_.abs_ = abs
del _

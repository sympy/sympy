
from sympy.core.basic import Basic, S, C, sympify
from sympy.core.function import Function, Lambda
from sympy.core.cache import cacheit

from sympy.functions.elementary.miscellaneous import sqrt

from sympy.utilities.decorator import deprecated

###############################################################################
########################### HYPERBOLIC FUNCTIONS ##############################
###############################################################################

class sinh(Function):
    """
    Usage
    =====
      sinh(x) -> Returns the hyperbolic sine of x
    """
    nargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return cosh(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return asinh

    @classmethod
    def _eval_apply_subs(self, *args):
        return

    @classmethod
    @deprecated
    def canonize(cls, arg):
        return eval(cls, arg)

    @classmethod
    def eval(cls, arg):
        arg = sympify(arg)

        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.Infinity
            elif arg is S.NegativeInfinity:
                return S.NegativeInfinity
            elif arg is S.Zero:
                return S.Zero
            elif arg.is_negative:
                return -cls(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return S.ImaginaryUnit * C.sin(i_coeff)
            else:
                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return -cls(-arg)

            if isinstance(arg, asinh):
                return arg.args[0]

            if isinstance(arg, acosh):
                x = arg.args[0]
                return sqrt(x-1) * sqrt(x+1)

            if isinstance(arg, atanh):
                x = arg.args[0]
                return x/sqrt(1-x**2)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)

            if len(previous_terms) > 2:
                p = previous_terms[-2]
                return p * x**2 / (n*(n-1))
            else:
                return x**(n) / C.Factorial(n)

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def _eval_expand_complex(self, deep=True, **hints):
        if self.args[0].is_real:
            if deep:
                hints['complex'] = False
                return self.expand(deep, **hints)
            else:
                return self
        if deep:
            re, im = self.args[0].expand(deep, **hints).as_real_imag()
        else:
            re, im = self.args[0].as_real_imag()
        return sinh(re)*C.cos(im) + cosh(re)*C.sin(im)*S.ImaginaryUnit

    def _eval_rewrite_as_exp(self, arg):
        return (S.Exp(arg) - S.Exp(-arg)) / 2

    def _eval_rewrite_as_cosh(self, arg):
        return -S.ImaginaryUnit*cosh(arg + S.Pi*S.ImaginaryUnit/2)

    def _eval_rewrite_as_tanh(self, arg):
        tanh_half = tanh(S.Half*arg)
        return 2*tanh_half/(1 - tanh_half**2)

    def _eval_rewrite_as_coth(self, arg):
        coth_half = coth(S.Half*arg)
        return 2*coth_half/(coth_half**2 - 1)

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if C.Order(1,x).contains(arg):
            return S.One
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
        return sage.sinh(self.args[0]._sage_())

class cosh(Function):
    """
    Usage
    =====
      cosh(x) -> Returns the hyperbolic cosine of x
    """
    nargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return sinh(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return acosh

    @classmethod
    def _eval_apply_subs(self, *args):
        return

    @classmethod
    @deprecated
    def canonize(cls, arg):
        return cls.eval(arg)

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
            elif arg is S.Zero:
                return S.One
            elif arg.is_negative:
                return cls(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return C.cos(i_coeff)
            else:
                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return cls(-arg)

            if isinstance(arg, asinh):
                return sqrt(1+arg.args[0]**2)

            if isinstance(arg, acosh):
                return arg.args[0]

            if isinstance(arg, atanh):
                return 1/sqrt(1-arg.args[0]**2)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0 or n % 2 == 1:
            return S.Zero
        else:
            x = sympify(x)

            if len(previous_terms) > 2:
                p = previous_terms[-2]
                return p * x**2 / (n*(n-1))
            else:
                return x**(n)/C.Factorial(n)

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def _eval_expand_complex(self, deep=True, **hints):
        if self.args[0].is_real:
            if deep:
                hints['complex'] = False
                return self.expand(deep, **hints)
            else:
                return self
        if deep:
            re, im = self.args[0].expand(deep, **hints).as_real_imag()
        else:
            re, im = self.args[0].as_real_imag()
        return cosh(re)*C.cos(im) + sinh(re)*C.sin(im)*S.ImaginaryUnit

    def _eval_rewrite_as_exp(self, arg):
        return (S.Exp(arg) + S.Exp(-arg)) / 2

    def _eval_rewrite_as_sinh(self, arg):
        return -S.ImaginaryUnit*sinh(arg + S.Pi*S.ImaginaryUnit/2)

    def _eval_rewrite_as_tanh(self, arg):
        tanh_half = tanh(S.Half*arg)**2
        return (1+tanh_half)/(1-tanh_half)

    def _eval_rewrite_as_coth(self, arg):
        coth_half = coth(S.Half*arg)**2
        return (coth_half+1)/(coth_half-1)

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if C.Order(1,x).contains(arg):
            return S.One
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
        return sage.cosh(self.args[0]._sage_())

class tanh(Function):
    """
    Usage
    =====
      tanh(x) -> Returns the hyperbolic tangent of x
    """
    nargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return S.One - tanh(self.args[0])**2
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return atanh

    @classmethod
    def _eval_apply_subs(self, *args):
        return

    @classmethod
    @deprecated
    def canonize(cls, arg):
        return cls.eval(arg)

    @classmethod
    def eval(cls, arg):
        arg = sympify(arg)

        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.One
            elif arg is S.NegativeInfinity:
                return S.NegativeOne
            elif arg is S.Zero:
                return S.Zero
            elif arg.is_negative:
                return -cls(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return S.ImaginaryUnit * C.tan(i_coeff)
            else:
                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return -cls(-arg)

            if isinstance(arg, asinh):
                x = arg.args[0]
                return x/sqrt(1+x**2)

            if isinstance(arg, acosh):
                x = arg.args[0]
                return sqrt(x-1) * sqrt(x+1) / x

            if isinstance(arg, atanh):
                return arg.args[0]

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)

            a = 2**(n+1)

            B = C.bernoulli(n+1)
            F = C.Factorial(n+1)

            return a*(a-1) * B/F * x**n

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def _eval_expand_complex(self, deep=True, **hints):
        if self.args[0].is_real:
            if deep:
                hints['complex'] = False
                return self.expand(deep, **hints)
            else:
                return self
        if deep:
            re, im = self.args[0].expand(deep, **hints).as_real_imag()
        else:
            re, im = self.args[0].as_real_imag()
        denom = sinh(re)**2 + C.cos(im)**2
        return (sinh(re)*cosh(re) + \
            S.ImaginaryUnit*C.sin(im)*C.cos(im))/denom

    def _eval_rewrite_as_exp(self, arg):
        neg_exp, pos_exp = S.Exp(-arg), S.Exp(arg)
        return (pos_exp-neg_exp)/(pos_exp+neg_exp)

    def _eval_rewrite_as_sinh(self, arg):
        return S.ImaginaryUnit*sinh(arg)/sinh(S.Pi*S.ImaginaryUnit/2 - arg)

    def _eval_rewrite_as_cosh(self, arg):
        return S.ImaginaryUnit*cosh(S.Pi*S.ImaginaryUnit/2 - arg)/cosh(arg)

    def _eval_rewrite_as_coth(self, arg):
        return 1/coth(arg)

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if C.Order(1,x).contains(arg):
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
        return sage.tanh(self.args[0]._sage_())

class coth(Function):
    """
    Usage
    =====
      coth(x) -> Returns the hyperbolic cotangent of x
    """
    nargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/sinh(self.args[0])**2
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return acoth

    @classmethod
    def _eval_apply_subs(self, *args):
        return

    @classmethod
    @deprecated
    def canonize(cls, arg):
        return cls.eval(arg)

    @classmethod
    def eval(cls, arg):
        arg = sympify(arg)

        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.One
            elif arg is S.NegativeInfinity:
                return S.NegativeOne
            elif arg is S.Zero:
                return S.Zero
            elif arg.is_negative:
                return -cls(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return -S.ImaginaryUnit * C.cot(i_coeff)
            else:
                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return -cls(-arg)

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
            F = C.Factorial(n+1)

            return 2**(n+1) * B/F * x**n

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def _eval_expand_complex(self, deep=True, **hints):
        if self.args[0].is_real:
            if deep:
                return self.expand(deep, **hints)
            else:
                return self
        if deep:
            re, im = self.args[0].expand(deep, **hints).as_real_imag()
        else:
            re, im = self.args[0].as_real_imag()
        denom = sinh(re)**2 + C.sin(im)**2
        return (sinh(re)*cosh(re) - \
            S.ImaginaryUnit*C.sin(im)*C.cos(im))/denom

    def _eval_rewrite_as_exp(self, arg):
        neg_exp, pos_exp = S.Exp(-arg), S.Exp(arg)
        return (pos_exp+neg_exp)/(pos_exp-neg_exp)

    def _eval_rewrite_as_sinh(self, arg):
        return -S.ImaginaryUnit*sinh(S.Pi*S.ImaginaryUnit/2 - arg)/sinh(arg)

    def _eval_rewrite_as_cosh(self, arg):
        return -S.ImaginaryUnit*cosh(arg)/cosh(S.Pi*S.ImaginaryUnit/2 - arg)

    def _eval_rewrite_as_tanh(self, arg):
        return 1/tanh(arg)

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if C.Order(1,x).contains(arg):
            return S.One
        else:
            return self.func(arg)

    def _sage_(self):
        import sage.all as sage
        return sage.coth(self.args[0]._sage_())


###############################################################################
############################# HYPERBOLIC INVERSES #############################
###############################################################################

class asinh(Function):
    """
    Usage
    =====
      asinh(x) -> Returns the inverse hyperbolic sine of x
    """
    nargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return (self.args[0]**2 + 1)**(-S.Half)
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
        arg = sympify(arg)

        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.Infinity
            elif arg is S.NegativeInfinity:
                return S.NegativeInfinity
            elif arg is S.Zero:
                return S.Zero
            elif arg is S.One:
                return C.log(2**S.Half + 1)
            elif arg is S.NegativeOne:
                return C.log(2**S.Half - 1)
            elif arg.is_negative:
                return -cls(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return S.ImaginaryUnit * C.asin(i_coeff)
            else:
                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return -cls(-arg)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)

            if len(previous_terms) > 2:
                p = previous_terms[-2]
                return -p * (n-2)**2/(k*(k-1)) * x**2
            else:
                k = (n - 1) // 2

                R = C.RisingFactorial(S.Half, k)
                F = C.Factorial(k)

                return (-1)**k * R / F * x**n / n

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if C.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _sage_(self):
        import sage.all as sage
        return sage.asinh(self.args[0]._sage_())

class acosh(Function):
    """
    Usage
    =====
      acosh(x) -> Returns the inverse hyperbolic cosine of x
    """
    nargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return (self.args[0]**2 - 1)**(-S.Half)
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
        arg = sympify(arg)

        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.Infinity * S.ImaginaryUnit
            elif arg is S.NegativeInfinity:
                return S.NegativeInfinity * S.ImaginaryUnit
            elif arg is S.Zero:
                return S.Pi*S.ImaginaryUnit / 2
            elif arg is S.One:
                return S.Zero
            elif arg is S.NegativeOne:
                return S.Pi*S.ImaginaryUnit
            else:
                cst_table = {
                    S.Half       : S.Pi/3,
                    -S.Half      : 2*S.Pi/3,
                    sqrt(2)/2    : S.Pi/4,
                    -sqrt(2)/2   : 3*S.Pi/4,
                    1/sqrt(2)    : S.Pi/4,
                    -1/sqrt(2)   : 3*S.Pi/4,
                    sqrt(3)/2    : S.Pi/6,
                    -sqrt(3)/2   : 5*S.Pi/6,
                }

                if arg in cst_table:
                    return cst_table[arg]*S.ImaginaryUnit

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n == 0:
            return S.Pi*S.ImaginaryUnit / 2
        elif n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)

            if len(previous_terms) > 2:
                p = previous_terms[-2]
                return p * (n-2)**2/(k*(k-1)) * x**2
            else:
                k = (n - 1) // 2

                R = C.RisingFactorial(S.Half, k)
                F = C.Factorial(k)

                return -R / F * S.ImaginaryUnit * x**n / n

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if C.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _sage_(self):
        import sage.all as sage
        return sage.acosh(self.args[0]._sage_())

class atanh(Function):
    """
    Usage
    =====
      atanh(x) -> Returns the inverse hyperbolic tangent of x
    """
    nargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/(1-self.args[0]**2)
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
        arg = sympify(arg)

        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Zero:
                return S.Zero
            elif arg is S.One:
                return S.Infinity
            elif arg is S.NegativeOne:
                return S.NegativeInfinity
            elif arg.is_negative:
                return -cls(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return S.ImaginaryUnit * C.atan(i_coeff)
            else:
                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return -cls(-arg)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = sympify(x)
            return x**n / n

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if C.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _sage_(self):
        import sage.all as sage
        return sage.atanh(self.args[0]._sage_())

class acoth(Function):
    """
    Usage
    =====
      acoth(x) -> Returns the inverse hyperbolic cotangent of x
    """
    nargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/(1-self.args[0]**2)
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
        arg = sympify(arg)

        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Infinity:
                return S.Zero
            elif arg is S.NegativeInfinity:
                return S.Zero
            elif arg is S.Zero:
                return S.Pi*S.ImaginaryUnit / 2
            elif arg is S.One:
                return S.Infinity
            elif arg is S.NegativeOne:
                return S.NegativeInfinity
            elif arg.is_negative:
                return -cls(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return -S.ImaginaryUnit * C.acot(i_coeff)
            else:
                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return -cls(-arg)

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

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if C.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _sage_(self):
        import sage.all as sage
        return sage.acoth(self.args[0]._sage_())

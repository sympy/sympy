
from sympy.core.basic import Basic, S, cache_it, cache_it_immutable
from sympy.core.function import DefinedFunction, Apply, Lambda

###############################################################################
########################### HYPERBOLIC FUNCTIONS ##############################
###############################################################################

class Sinh(DefinedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return S.Cosh
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return S.ASinh

    def _eval_apply(self, arg):
        arg = Basic.sympify(arg)

        if isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.NaN):
                return S.NaN
            elif isinstance(arg, Basic.Infinity):
                return S.Infinity
            elif isinstance(arg, Basic.NegativeInfinity):
                return S.NegativeInfinity
            elif isinstance(arg, Basic.Zero):
                return S.Zero
            elif arg.is_negative:
                return -self(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return S.ImaginaryUnit * S.Sin(i_coeff)
            else:
                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return -self(-arg)

    def _eval_apply_evalf(self, arg):
        arg = arg.evalf()

        if isinstance(arg, Basic.Number):
            return arg.sinh()

    @cache_it_immutable
    def taylor_term(self, n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = Basic.sympify(x)

            if len(previous_terms) > 2:
                p = previous_terms[-2]
                return p * x**2 / (n*(n-1))
            else:
                return x**(n) / S.Factorial(n)

class ApplySinh(Apply):

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def _eval_expand_complex(self, *args):
        if self.args[0].is_real:
            return self
        re, im = self.args[0].as_real_imag()
        return S.Sinh(re)*S.Cos(im) + S.Cosh(re)*S.Sin(im)*S.ImaginaryUnit

    def _eval_rewrite_as_exp(self, arg):
        return (S.Exp(arg) - S.Exp(-arg)) / 2

    def _eval_rewrite_as_cosh(self, arg):
        return -S.ImaginaryUnit*S.Cosh(arg + S.Pi*S.ImaginaryUnit/2)

    def _eval_rewrite_as_tanh(self, arg):
        tanh_half = S.Tanh(S.Half*arg)
        return 2*tanh_half/(1 - tanh_half**2)

    def _eval_rewrite_as_coth(self, arg):
        coth_half = S.Coth(S.Half*arg)
        return 2*coth_half/(coth_half**2 - 1)

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return S.One
        else:
            return self.func(arg)

    def _eval_is_real(self):
        return self.args[0].is_real
    def _eval_is_bounded(self):
        arg = self.args[0]
        if arg.is_imaginary:
            return True

class Cosh(DefinedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return S.Sinh
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return S.ACosh

    def _eval_apply(self, arg):
        arg = Basic.sympify(arg)

        if isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.NaN):
                return S.NaN
            elif isinstance(arg, Basic.Infinity):
                return S.Infinity
            elif isinstance(arg, Basic.NegativeInfinity):
                return S.Infinity
            elif isinstance(arg, Basic.Zero):
                return S.One
            elif arg.is_negative:
                return self(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return S.Cos(i_coeff)
            else:
                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return self(-arg)

    def _eval_apply_evalf(self, arg):
        arg = arg.evalf()

        if isinstance(arg, Basic.Number):
            return arg.cosh()

    @cache_it_immutable
    def taylor_term(self, n, x, *previous_terms):
        if n < 0 or n % 2 == 1:
            return S.Zero
        else:
            x = Basic.sympify(x)

            if len(previous_terms) > 2:
                p = previous_terms[-2]
                return p * x**2 / (n*(n-1))
            else:
                return x**(n)/S.Factorial(n)

class ApplyCosh(Apply):
    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def _eval_expand_complex(self, *args):
        if self.args[0].is_real:
            return self
        re, im = self.args[0].as_real_imag()
        return S.Cosh(re)*S.Cos(im) + S.Sinh(re)*S.Sin(im)*S.ImaginaryUnit

    def _eval_rewrite_as_exp(self, arg):
        return (S.Exp(arg) + S.Exp(-arg)) / 2

    def _eval_rewrite_as_sinh(self, arg):
        return -S.ImaginaryUnit*S.Sinh(arg + S.Pi*S.ImaginaryUnit/2)

    def _eval_rewrite_as_tanh(self, arg):
        tanh_half = S.Tanh(S.Half*arg)**2
        return (1+tanh_half)/(1-tanh_half)

    def _eval_rewrite_as_coth(self, arg):
        coth_half = S.Coth(S.Half*arg)**2
        return (coth_half+1)/(coth_half-1)

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return S.One
        else:
            return self.func(arg)

    def _eval_is_real(self):
        return self.args[0].is_real
    def _eval_is_bounded(self):
        arg = self.args[0]
        if arg.is_imaginary:
            return True

class Tanh(DefinedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return S.One - S.Tanh**2
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return S.ATanh

    def _eval_apply(self, arg):
        arg = Basic.sympify(arg)

        if isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.NaN):
                return S.NaN
            elif isinstance(arg, Basic.Infinity):
                return S.One
            elif isinstance(arg, Basic.NegativeInfinity):
                return S.NegativeOne
            elif isinstance(arg, Basic.Zero):
                return S.Zero
            elif arg.is_negative:
                return -self(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return S.ImaginaryUnit * S.Tan(i_coeff)
            else:
                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return -self(-arg)

    def _eval_apply_evalf(self, arg):
        arg = arg.evalf()

        if isinstance(arg, Basic.Number):
            return arg.tanh()

    @cache_it_immutable
    def taylor_term(self, n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = Basic.sympify(x)

            a = 2**(n+1)

            B = S.Bernoulli(n+1)
            F = S.Factorial(n+1)

            return a*(a-1) * B/F * x**n

class ApplyTanh(Apply):

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def _eval_expand_complex(self, *args):
        if self.args[0].is_real:
            return self
        re, im = self.args[0].as_real_imag()
        denom = S.Sinh(re)**2 + S.Cos(im)**2
        return (S.Sinh(re)*S.Cosh(re) + \
            S.ImaginaryUnit*S.Sin(im)*S.Cos(im))/denom

    def _eval_rewrite_as_exp(self, arg):
        neg_exp, pos_exp = S.Exp(-arg), S.Exp(arg)
        return (pos_exp-neg_exp)/(pos_exp+neg_exp)

    def _eval_rewrite_as_sinh(self, arg):
        return S.ImaginaryUnit*S.Sinh(arg)/S.Sinh(S.Pi*S.ImaginaryUnit/2 - arg)

    def _eval_rewrite_as_cosh(self, arg):
        return S.ImaginaryUnit*S.Cosh(S.Pi*S.ImaginaryUnit/2 - arg)/S.Cosh(arg)

    def _eval_rewrite_as_coth(self, arg):
        return 1/S.Coth(arg)

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return S.One
        else:
            return self.func(arg)

    def _eval_is_real(self):
        return self.args[0].is_real
    def _eval_is_bounded(self):
        arg = self.args[0]
        if arg.is_real:
            return True

class Coth(DefinedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/S.Sinh**2
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return S.ACoth

    def _eval_apply(self, arg):
        arg = Basic.sympify(arg)

        if isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.NaN):
                return S.NaN
            elif isinstance(arg, Basic.Infinity):
                return S.One
            elif isinstance(arg, Basic.NegativeInfinity):
                return S.NegativeOne
            elif isinstance(arg, Basic.Zero):
                return S.Zero
            elif arg.is_negative:
                return -self(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return -S.ImaginaryUnit * Basic.cot(i_coeff)
            else:
                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return -self(-arg)

    def _eval_apply_evalf(self, arg):
        arg = arg.evalf()

        if isinstance(arg, Basic.Number):
            return arg.coth()

    @cache_it_immutable
    def taylor_term(self, n, x, *previous_terms):
        if n == 0:
            return 1 / Basic.sympify(x)
        elif n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = Basic.sympify(x)

            B = S.Bernoulli(n+1)
            F = S.Factorial(n+1)

            return 2**(n+1) * B/F * x**n

class ApplyCoth(Apply):

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def _eval_expand_complex(self, *args):
        if self.args[0].is_real:
            return self
        re, im = self.args[0].as_real_imag()
        denom = S.Sinh(re)**2 + S.Sin(im)**2
        return (S.Sinh(re)*S.Cosh(re) - \
            S.ImaginaryUnit*S.Sin(im)*S.Cos(im))/denom

    def _eval_rewrite_as_exp(self, arg):
        neg_exp, pos_exp = S.Exp(-arg), S.Exp(arg)
        return (pos_exp+neg_exp)/(pos_exp-neg_exp)

    def _eval_rewrite_as_sinh(self, arg):
        return -S.ImaginaryUnit*S.Sinh(S.Pi*S.ImaginaryUnit/2 - arg)/S.Sinh(arg)

    def _eval_rewrite_as_cosh(self, arg):
        return -S.ImaginaryUnit*S.Cosh(arg)/S.Cosh(S.Pi*S.ImaginaryUnit/2 - arg)

    def _eval_rewrite_as_tanh(self, arg):
        return 1/S.Tanh(arg)

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return S.One
        else:
            return self.func(arg)

Basic.singleton['sinh'] = Sinh
Basic.singleton['cosh'] = Cosh
Basic.singleton['tanh'] = Tanh
Basic.singleton['coth'] = Coth

###############################################################################
############################# HYPERBOLIC INVERSES #############################
###############################################################################

class ASinh(DefinedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            z = Basic.Symbol('z', dummy=True)
            return Lambda((1 + z**2)**(-S.Half), z)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_apply(self, arg):
        arg = Basic.sympify(arg)

        if isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.NaN):
                return S.NaN
            elif isinstance(arg, Basic.Infinity):
                return S.Infinity
            elif isinstance(arg, Basic.NegativeInfinity):
                return S.NegativeInfinity
            elif isinstance(arg, Basic.Zero):
                return S.Zero
            elif isinstance(arg, Basic.One):
                return S.Log(S.Sqrt(2) + 2)
            elif isinstance(arg, Basic.NegativeOne):
                return S.Log(S.Sqrt(2) - 2)
            elif arg.is_negative:
                return -self(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return S.ImaginaryUnit * S.ASin(i_coeff)
            else:
                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return -self(-arg)

    def _eval_apply_evalf(self, arg):
        arg = arg.evalf()

        if isinstance(arg, Basic.Number):
            return arg.asinh()

    @cache_it_immutable
    def taylor_term(self, n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = Basic.sympify(x)

            if len(previous_terms) > 2:
                p = previous_terms[-2]
                return -p * (n-2)**2/(k*(k-1)) * x**2
            else:
                k = (n - 1) // 2

                R = S.RisingFactorial(S.Half, k)
                F = S.Factorial(k)

                return (-1)**k * R / F * x**n / n

class ApplyASinh(Apply):

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

class ACosh(DefinedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            z = Basic.Symbol('z', dummy=True)
            return Lambda(((z-1)*(z+1))**(-1), z)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_apply(self, arg):
        arg = Basic.sympify(arg)

        if isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.NaN):
                return S.NaN
            elif isinstance(arg, Basic.Infinity):
                return S.Infinity * S.ImaginaryUnit
            elif isinstance(arg, Basic.NegativeInfinity):
                return S.NegativeInfinity * S.ImaginaryUnit
            elif isinstance(arg, Basic.Zero):
                return S.Pi*S.ImaginaryUnit / 2
            elif isinstance(arg, Basic.One):
                return S.Zero
            elif isinstance(arg, Basic.NegativeOne):
                return S.Pi*S.ImaginaryUnit
            else:
                cst_table = {
                    S.Half       : S.Pi/3,
                    -S.Half      : 2*S.Pi/3,
                    S.Sqrt(2)/2  : S.Pi/4,
                    -S.Sqrt(2)/2 : 3*S.Pi/4,
                    1/S.Sqrt(2)  : S.Pi/4,
                    -1/S.Sqrt(2) : 3*S.Pi/4,
                    S.Sqrt(3)/2  : S.Pi/6,
                    -S.Sqrt(3)/2 : 5*S.Pi/6,
                }

                if arg in cst_table:
                    return cst_table[arg]*S.ImaginaryUnit

    def _eval_apply_evalf(self, arg):
        arg = arg.evalf()

        if isinstance(arg, Basic.Number):
            return arg.acosh()

    @cache_it_immutable
    def taylor_term(self, n, x, *previous_terms):
        if n == 0:
            return S.Pi*S.ImaginaryUnit / 2
        elif n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = Basic.sympify(x)

            if len(previous_terms) > 2:
                p = previous_terms[-2]
                return p * (n-2)**2/(k*(k-1)) * x**2
            else:
                k = (n - 1) // 2

                R = S.RisingFactorial(S.Half, k)
                F = S.Factorial(k)

                return -R / F * S.ImaginaryUnit * x**n / n

class ApplyACosh(Apply):

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

class ATanh(DefinedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            z = Basic.Symbol('z', dummy=True)
            return Lambda((z-1)**2, z)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_apply(self, arg):
        arg = Basic.sympify(arg)

        if isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.NaN):
                return S.NaN
            elif isinstance(arg, Basic.Zero):
                return S.Zero
            elif isinstance(arg, Basic.One):
                return S.Infinity
            elif isinstance(arg, Basic.NegativeOne):
                return S.NegativeInfinity
            elif arg.is_negative:
                return -self(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return S.ImaginaryUnit * S.ATan(i_coeff)
            else:
                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return -self(-arg)

    def _eval_apply_evalf(self, arg):
        arg = arg.evalf()

        if isinstance(arg, Basic.Number):
            return arg.atanh()

    @cache_it_immutable
    def taylor_term(self, n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = Basic.sympify(x)
            return x**n / n

class ApplyATanh(Apply):

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

class ACoth(DefinedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            z = Basic.Symbol('z', dummy=True)
            return Lambda((z-1)**2, z)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_apply(self, arg):
        arg = Basic.sympify(arg)

        if isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.NaN):
                return S.NaN
            elif isinstance(arg, Basic.Infinity):
                return S.Zero
            elif isinstance(arg, Basic.NegativeInfinity):
                return S.Zero
            elif isinstance(arg, Basic.Zero):
                return S.Pi*S.ImaginaryUnit / 2
            elif isinstance(arg, Basic.One):
                return S.Infinity
            elif isinstance(arg, Basic.NegativeOne):
                return S.NegativeInfinity
            elif arg.is_negative:
                return -self(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return -S.ImaginaryUnit * S.ACot(i_coeff)
            else:
                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return -self(-arg)

    def _eval_apply_evalf(self, arg):
        arg = arg.evalf()

        if isinstance(arg, Basic.Number):
            return arg.acoth()

    @cache_it_immutable
    def taylor_term(self, n, x, *previous_terms):
        if n == 0:
            return S.Pi*S.ImaginaryUnit / 2
        elif n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = Basic.sympify(x)
            return x**n / n

class ApplyACoth(Apply):

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

Basic.singleton['asinh'] = ASinh
Basic.singleton['acosh'] = ACosh
Basic.singleton['atanh'] = ATanh
Basic.singleton['acoth'] = ACoth


from sympy.core.basic import Basic, S, cache_it, cache_it_immutable
from sympy.core.function import SingleValuedFunction, Lambda

###############################################################################
########################### HYPERBOLIC FUNCTIONS ##############################
###############################################################################

class sinh(SingleValuedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return cosh(self[0])
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return asinh

    @classmethod
    def _eval_apply_subs(self, *args):
        return

    @classmethod
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
                return S.ImaginaryUnit * Basic.sin(i_coeff)
            else:
                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return -self(-arg)

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

    def _eval_conjugate(self):
        return self.func(self[0].conjugate())

    def _eval_expand_complex(self, *args):
        if self[0].is_real:
            return self
        re, im = self[0].as_real_imag()
        return sinh(re)*Basic.cos(im) + cosh(re)*Basic.sin(im)*S.ImaginaryUnit

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
        arg = self[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return S.One
        else:
            return self.func(arg)

    def _eval_is_real(self):
        return self[0].is_real
    def _eval_is_bounded(self):
        arg = self[0]
        if arg.is_imaginary:
            return True

class cosh(SingleValuedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return sinh(self[0])
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return acosh

    @classmethod
    def _eval_apply_subs(self, *args):
        return

    @classmethod
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
                return Basic.cos(i_coeff)
            else:
                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return self(-arg)

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

    def _eval_conjugate(self):
        return self.func(self[0].conjugate())

    def _eval_expand_complex(self, *args):
        if self[0].is_real:
            return self
        re, im = self[0].as_real_imag()
        return cosh(re)*Basic.cos(im) + sinh(re)*Basic.sin(im)*S.ImaginaryUnit

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
        arg = self[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return S.One
        else:
            return self.func(arg)

    def _eval_is_real(self):
        return self[0].is_real
    def _eval_is_bounded(self):
        arg = self[0]
        if arg.is_imaginary:
            return True

class tanh(SingleValuedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return S.One - tanh(self[0])**2
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return atanh

    @classmethod
    def _eval_apply_subs(self, *args):
        return

    @classmethod
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
                return S.ImaginaryUnit * Basic.tan(i_coeff)
            else:
                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return -self(-arg)

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

    def _eval_conjugate(self):
        return self.func(self[0].conjugate())

    def _eval_expand_complex(self, *args):
        if self[0].is_real:
            return self
        re, im = self[0].as_real_imag()
        denom = sinh(re)**2 + Basic.cos(im)**2
        return (sinh(re)*cosh(re) + \
            S.ImaginaryUnit*Basic.sin(im)*Basic.cos(im))/denom

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
        arg = self[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return S.One
        else:
            return self.func(arg)

    def _eval_is_real(self):
        return self[0].is_real
    def _eval_is_bounded(self):
        arg = self[0]
        if arg.is_real:
            return True

class coth(SingleValuedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/sinh(self[0])**2
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return acoth

    @classmethod
    def _eval_apply_subs(self, *args):
        return

    @classmethod
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

    def _eval_conjugate(self):
        return self.func(self[0].conjugate())

    def _eval_expand_complex(self, *args):
        if self[0].is_real:
            return self
        re, im = self[0].as_real_imag()
        denom = sinh(re)**2 + Basic.sin(im)**2
        return (sinh(re)*cosh(re) - \
            S.ImaginaryUnit*Basic.sin(im)*Basic.cos(im))/denom

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
        arg = self[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return S.One
        else:
            return self.func(arg)


###############################################################################
############################# HYPERBOLIC INVERSES #############################
###############################################################################

class asinh(SingleValuedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return (self[0]**2 + 1)**(-S.Half)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def _eval_apply_subs(self, *args):
        return

    @classmethod
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

    def _eval_as_leading_term(self, x):
        arg = self[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

class acosh(SingleValuedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return (self[0]**2 - 1)**(-S.Half)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def _eval_apply_subs(self, *args):
        return

    @classmethod
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

    def _eval_as_leading_term(self, x):
        arg = self[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

class atanh(SingleValuedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/(1-self[0]**2)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def _eval_apply_subs(self, *args):
        return

    @classmethod
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

    @cache_it_immutable
    def taylor_term(self, n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = Basic.sympify(x)
            return x**n / n

    def _eval_as_leading_term(self, x):
        arg = self[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

class acoth(SingleValuedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/(1-self[0]**2)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def _eval_apply_subs(self, *args):
        return

    @classmethod
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

    @cache_it_immutable
    def taylor_term(self, n, x, *previous_terms):
        if n == 0:
            return S.Pi*S.ImaginaryUnit / 2
        elif n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = Basic.sympify(x)
            return x**n / n

    def _eval_as_leading_term(self, x):
        arg = self[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)


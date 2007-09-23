
from sympy.core.basic import Basic, S, cache_it, cache_it_immutable
from sympy.core.function import DefinedFunction, Apply, Lambda, SingleValuedFunction

###############################################################################
########################## TRIGONOMETRIC FUNCTIONS ############################
###############################################################################

class Sin(DefinedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return S.Cos
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return S.ASin

    def _eval_apply(self, arg):
        arg = Basic.sympify(arg)

        if isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.NaN):
                return S.NaN
            elif isinstance(arg, Basic.Zero):
                return S.Zero
            elif arg.is_negative:
                return -self(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return S.ImaginaryUnit * S.Sinh(i_coeff)
            else:
                pi_coeff = arg.as_coefficient(S.Pi)

                if pi_coeff is not None:
                    if pi_coeff.is_integer:
                        return S.Zero
                    elif isinstance(pi_coeff, Basic.Rational):
                        cst_table = {
                            2 : S.One,
                            3 : S.Half*S.Sqrt(3),
                            4 : S.Half*S.Sqrt(2),
                            6 : S.Half,
                        }

                        try:
                            result = cst_table[pi_coeff.q]

                            if (pi_coeff.p // pi_coeff.q) % 2 == 1:
                                return -result
                            else:
                                return result
                        except KeyError:
                            pass

                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return -self(-arg)

    def _eval_apply_evalf(self, arg):
        arg = arg.evalf()

        if isinstance(arg, Basic.Number):
            return arg.sin()

    @cache_it_immutable
    def taylor_term(self, n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = Basic.sympify(x)

            if len(previous_terms) > 2:
                p = previous_terms[-2]
                return -p * x**2 / (n*(n-1))
            else:
                return (-1)**(n//2) * x**(n)/S.Factorial(n)

class ApplySin(Apply):

    def _eval_rewrite_as_exp(self, arg):
        exp, I = S.Exp, S.ImaginaryUnit
        return (exp(arg*I) - exp(-arg*I)) / (2*I)

    def _eval_rewrite_as_cos(self, arg):
        return -S.Cos(arg + S.Pi/2)

    def _eval_rewrite_as_tan(self, arg):
        tan_half = S.Tan(S.Half*arg)
        return 2*tan_half/(1 + tan_half**2)

    def _eval_rewrite_as_cot(self, arg):
        cot_half = S.Cot(S.Half*arg)
        return 2*cot_half/(1 + cot_half**2)

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def _eval_expand_complex(self, *args):
        if self.args[0].is_real:
            return self
        re, im = self.args[0].as_real_imag()
        return S.Sin(re)*S.Cosh(im) + \
            S.ImaginaryUnit*S.Cos(re)*S.Sinh(im)

    def _eval_expand_trig(self, *args):
        arg = self.args[0].expand()
        cos, sin = S.Cos, S.Sin
        x = None
        if isinstance(arg, Basic.Add):
            x = arg[0]
            y = Basic.Add(*arg[1:])
        else:
            coeff, terms = arg.as_coeff_terms()
            if not isinstance(coeff, Basic.One) and isinstance(coeff, Basic.Integer) and terms:
                x = Basic.Mul(*terms)
                y = (coeff-1)*x
        if x is not None:
            return (sin(x)*cos(y) + sin(y)*cos(x)).expand(trig=True)
        return sin(arg)

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_is_real(self):
        return self.args[0].is_real

    def _eval_is_bounded(self):
        arg = self.args[0]
        if arg.is_real:
            return True

class Cos(DefinedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -S.Sin
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return S.ACos

    def _eval_apply(self, arg):
        arg = Basic.sympify(arg)

        if isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.NaN):
                return S.NaN
            elif isinstance(arg, Basic.Zero):
                return S.One
            elif arg.is_negative:
                return self(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return S.Cosh(i_coeff)
            else:
                pi_coeff = arg.as_coefficient(S.Pi)

                if pi_coeff is not None:
                    if isinstance(pi_coeff, Basic.Rational):
                        cst_table = {
                            1 : S.One,
                            2 : S.Zero,
                            3 : S.Half,
                            4 : S.Half*S.Sqrt(2),
                            6 : S.Half*S.Sqrt(3),
                        }

                        try:
                            result = cst_table[pi_coeff.q]

                            if (2*pi_coeff.p // pi_coeff.q) % 4 in (1, 2):
                                return -result
                            else:
                                return result
                        except KeyError:
                            pass

                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return self(-arg)

    def _eval_apply_evalf(self, arg):
        arg = arg.evalf()

        if isinstance(arg, Basic.Number):
            return arg.cos()

    @cache_it_immutable
    def taylor_term(self, n, x, *previous_terms):
        if n < 0 or n % 2 == 1:
            return S.Zero
        else:
            x = Basic.sympify(x)

            if len(previous_terms) > 2:
                p = previous_terms[-2]
                return -p * x**2 / (n*(n-1))
            else:
                return (-1)**(n//2)*x**(n)/S.Factorial(n)

class ApplyCos(Apply):

    def _eval_rewrite_as_exp(self, arg):
        exp, I = S.Exp, S.ImaginaryUnit
        return (exp(arg*I) + exp(-arg*I)) / 2

    def _eval_rewrite_as_sin(self, arg):
        return S.Sin(arg + S.Pi/2)

    def _eval_rewrite_as_tan(self, arg):
        tan_half = S.Tan(S.Half*arg)**2
        return (1-tan_half)/(1+tan_half)

    def _eval_rewrite_as_cot(self, arg):
        cot_half = S.Cot(S.Half*arg)**2
        return (cot_half-1)/(cot_half+1)

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def _eval_expand_complex(self, *args):
        if self.args[0].is_real:
            return self
        re, im = self.args[0].as_real_imag()
        return S.Cos(re)*S.Cosh(im) - \
            S.ImaginaryUnit*S.Sin(re)*S.Sinh(im)

    def _eval_expand_trig(self, *args):
        arg = self.args[0].expand()
        cos = S.Cos
        sin = S.Sin
        x = None
        if isinstance(arg, Basic.Add):
            x = arg[0]
            y = Basic.Add(*arg[1:])
            return (cos(x)*cos(y) - sin(y)*sin(x)).expand(trig=True)
        else:
            coeff, terms = arg.as_coeff_terms()
            if not isinstance(coeff, Basic.One) and isinstance(coeff, Basic.Integer) and terms:
                x = Basic.Mul(*terms)
                return Basic.ChebyshevT()(coeff, cos(x))
        return cos(arg)

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return S.One
        else:
            return self.func(arg)

    def _eval_is_bounded(self):
        arg = self.args[0]

        if arg.is_real:
            return True

class Tan(DefinedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex==1:
            return S.One + S.Tan**2
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return S.ATan

    def _eval_apply(self, arg):
        arg = Basic.sympify(arg)

        if isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.NaN):
                return S.NaN
            elif isinstance(arg, Basic.Zero):
                return S.Zero
            elif arg.is_negative:
                return -self(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return S.ImaginaryUnit * S.Tanh(i_coeff)
            else:
                pi_coeff = arg.as_coefficient(S.Pi)

                if pi_coeff is not None:
                    if pi_coeff.is_integer:
                        return S.Zero
                    elif isinstance(pi_coeff, Basic.Rational):
                        cst_table = {
                           #2 : S.ComplexInfinity,
                            3 : S.Sqrt(3),
                            4 : S.One,
                            6 : 1 / S.Sqrt(3),
                        }

                        try:
                            result = cst_table[pi_coeff.q]

                            if (2*pi_coeff.p // pi_coeff.q) % 4 in (1, 3):
                                return -result
                            else:
                                return result
                        except KeyError:
                            pass

                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return -self(-arg)

    def _eval_apply_evalf(self, arg):
        arg = arg.evalf()

        if isinstance(arg, Basic.Number):
            return arg.tan()

    @cache_it_immutable
    def taylor_term(self, n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = Basic.sympify(x)

            a, b = ((n-1)//2), 2**(n+1)

            B = S.Bernoulli(n+1)
            F = S.Factorial(n+1)

            return (-1)**a * b*(b-1) * B/F * x**n

class ApplyTan(Apply):

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def _eval_expand_complex(self, *args):
        if self.args[0].is_real:
            return self
        re, im = self.args[0].as_real_imag()
        denom = S.Cos(re)**2 + S.Sinh(im)**2
        return (S.Sin(re)*S.Cos(re) + \
            S.ImaginaryUnit*S.Sinh(im)*S.Cosh(im))/denom

    def _eval_expand_trig(self, *args):
        return self

    def _eval_rewrite_as_exp(self, arg):
        exp, I = S.Exp, S.ImaginaryUnit
        neg_exp, pos_exp = exp(-arg*I), exp(arg*I)
        return I*(neg_exp-pos_exp)/(neg_exp+pos_exp)

    def _eval_rewrite_as_sin(self, arg):
        return 2*S.Sin(x)**2/S.Sin(2*x)

    def _eval_rewrite_as_cos(self, arg):
        return -S.Cos(x + S.Pi/2)/S.Cos(x)

    def _eval_rewrite_as_cot(self, arg):
        return 1/S.Cot(arg)

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return S.One
        else:
            return self.func(arg)

    def _eval_is_bounded(self):
        arg = self.args[0]

        if arg.is_imaginary:
            return True

class cot(SingleValuedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -S.One - S.Cot**2
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return S.ACot

    @classmethod
    def _eval_apply_subs(self, *args):
        return

    @classmethod
    def _eval_apply(self, arg):
        arg = Basic.sympify(arg)

        if isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.NaN):
                return S.NaN
            #elif isinstance(arg, Basic.Zero):
            #    return S.ComplexInfinity
            elif arg.is_negative:
                return -self(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return -S.ImaginaryUnit * S.Coth(i_coeff)
            else:
                pi_coeff = arg.as_coefficient(S.Pi)

                if pi_coeff is not None:
                    #if pi_coeff.is_integer:
                    #    return S.ComplexInfinity
                    if isinstance(pi_coeff, Basic.Rational):
                        cst_table = {
                            2 : S.Zero,
                            3 : 1 / S.Sqrt(3),
                            4 : S.One,
                            6 : S.Sqrt(3)
                        }

                        try:
                            result = cst_table[pi_coeff.q]

                            if (2*pi_coeff.p // pi_coeff.q) % 4 in (1, 3):
                                return -result
                            else:
                                return result
                        except KeyError:
                            pass

                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return -self(-arg)

    def _eval_apply_evalf(self, arg):
        arg = arg.evalf()

        if isinstance(arg, Basic.Number):
            return arg.cot()

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

            return (-1)**((n+1)//2) * 2**(n+1) * B/F * x**n

    def _eval_conjugate(self):
        args = self[1:] #empty!?
        args = self._args
        print type(self), self.func, self.func(Basic.Symbol("x"))
        assert len(args) == 1
        return self.func(args[0].conjugate())

    def _eval_expand_complex(self, *args):
        if self.args[0].is_real:
            return self
        re, im = self.args[0].as_real_imag()
        denom = S.Sin(re)**2 + S.Sinh(im)**2
        return (S.Sin(re)*S.Cos(re) - \
            S.ImaginaryUnit*S.Sinh(im)*S.Cosh(im))/denom

    def _eval_rewrite_as_exp(self, arg):
        exp, I = S.Exp, S.ImaginaryUnit
        neg_exp, pos_exp = exp(-arg*I), exp(arg*I)
        return I*(pos_exp+neg_exp)/(pos_exp-neg_exp)

    def _eval_rewrite_as_sin(self, arg):
        return 2*S.Sin(2*x)/S.Sin(x)**2

    def _eval_rewrite_as_cos(self, arg):
        return -S.Cos(x)/S.Cos(x + S.Pi/2)

    def _eval_rewrite_as_tan(self, arg):
        return 1/S.Tan(arg)

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return S.One
        else:
            return self.func(arg)

Basic.singleton['sin'] = Sin
Basic.singleton['cos'] = Cos
Basic.singleton['tan'] = Tan

###############################################################################
########################### TRIGONOMETRIC INVERSES ############################
###############################################################################

class ASin(DefinedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            s = Basic.Symbol('x', dummy=True)
            return Lambda((1 - s**2)**(-S.Half), s)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_apply(self, arg):
        arg = Basic.sympify(arg)

        if isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.NaN):
                return S.NaN
            elif isinstance(arg, Basic.Infinity):
                return S.NegativeInfinity * S.ImaginaryUnit
            elif isinstance(arg, Basic.NegativeInfinity):
                return S.Infinity * S.ImaginaryUnit
            elif isinstance(arg, Basic.Zero):
                return S.Zero
            elif isinstance(arg, Basic.One):
                return S.Pi / 2
            elif isinstance(arg, Basic.NegativeOne):
                return -S.Pi / 2
            else:
                cst_table = {
                    S.Half       : 6,
                    -S.Half      : -6,
                    S.Sqrt(2)/2  : 4,
                    -S.Sqrt(2)/2 : -4,
                    1/S.Sqrt(2)  : 4,
                    -1/S.Sqrt(2) : -4,
                    S.Sqrt(3)/2  : 3,
                    -S.Sqrt(3)/2 : -3,
                }

                if arg in cst_table:
                    return S.Pi / cst_table[arg]
                elif arg.is_negative:
                    return -self(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return S.ImaginaryUnit * S.ASinh(i_coeff)
            else:
                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return -self(-arg)

    def _eval_apply_evalf(self, arg):
        arg = arg.evalf()

        if isinstance(arg, Basic.Number):
            return arg.asin()

    @cache_it_immutable
    def taylor_term(self, n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
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

                return R / F * x**n / n

class ApplyASin(Apply):

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

class ACos(DefinedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            s = Basic.Symbol('x', dummy=True)
            return Lambda(-(1 - s**2)**(-S.Half), s)
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
                return S.Pi / 2
            elif isinstance(arg, Basic.One):
                return S.Zero
            elif isinstance(arg, Basic.NegativeOne):
                return S.Pi
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
                    return cst_table[arg]

    def _eval_apply_evalf(self, arg):
        arg = arg.evalf()

        if isinstance(arg, Basic.Number):
            return arg.acos()

    @cache_it_immutable
    def taylor_term(self, n, x, *previous_terms):
        if n == 0:
            return S.Pi / 2
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

                return -R / F * x**n / n

class ApplyACos(Apply):

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

class ATan(DefinedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            s = Basic.Symbol('x', dummy=True)
            return Lambda(1 / (1 + s**2), s)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_apply(self, arg):
        arg = Basic.sympify(arg)

        if isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.NaN):
                return S.NaN
            elif isinstance(arg, Basic.Infinity):
                return S.Pi / 2
            elif isinstance(arg, Basic.NegativeInfinity):
                return -S.Pi / 2
            elif isinstance(arg, Basic.Zero):
                return S.Zero
            elif isinstance(arg, Basic.One):
                return S.Pi / 4
            elif isinstance(arg, Basic.NegativeOne):
                return -S.Pi / 4
            else:
                cst_table = {
                    S.Sqrt(3)/3  : 6,
                    -S.Sqrt(3)/3 : -6,
                    1/S.Sqrt(3)  : 6,
                    -1/S.Sqrt(3) : -6,
                    S.Sqrt(3)    : 3,
                    -S.Sqrt(3)   : -3,
                }

                if arg in cst_table:
                    return S.Pi / cst_table[arg]
                elif arg.is_negative:
                    return -self(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return S.ImaginaryUnit * S.ATanh(i_coeff)
            else:
                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return -self(-arg)

    def _eval_apply_evalf(self, arg):
        arg = arg.evalf()

        if isinstance(arg, Basic.Number):
            return arg.atan()

    @cache_it_immutable
    def taylor_term(self, n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = Basic.sympify(x)
            return (-1)**((n-1)//2) * x**n / n

class ApplyATan(Apply):

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

class acot(SingleValuedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            s = Basic.Symbol('x', dummy=True)
            return Lambda(-1 / (1 + s**2), s)
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
                return S.Pi/ 2
            elif isinstance(arg, Basic.One):
                return S.Pi / 4
            elif isinstance(arg, Basic.NegativeOne):
                return -S.Pi / 4
            else:
                cst_table = {
                    S.Sqrt(3)/3  : 3,
                    -S.Sqrt(3)/3 : -3,
                    1/S.Sqrt(3)  : 3,
                    -1/S.Sqrt(3) : -3,
                    S.Sqrt(3)    : 6,
                    -S.Sqrt(3)   : -6,
                }

                if arg in cst_table:
                    return S.Pi / cst_table[arg]
                elif arg.is_negative:
                    return -self(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return -S.ImaginaryUnit * S.ACoth(i_coeff)
            else:
                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return -self(-arg)

    def _eval_apply_evalf(self, arg):
        arg = arg.evalf()

        if isinstance(arg, Basic.Number):
            return arg.acot()

    @cache_it_immutable
    def taylor_term(self, n, x, *previous_terms):
        if n == 0:
            return S.Pi / 2 # FIX THIS
        elif n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = Basic.sympify(x)
            return (-1)**((n+1)//2) * x**n / n

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

Basic.singleton['asin'] = ASin
Basic.singleton['acos'] = ACos
Basic.singleton['atan'] = ATan

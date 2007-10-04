
from sympy.core.basic import Basic, S, cache_it, cache_it_immutable
from sympy.core.function import Lambda, SingleValuedFunction

###############################################################################
########################## TRIGONOMETRIC FUNCTIONS ############################
###############################################################################

class sin(SingleValuedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return cos(self[0])
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return asin

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
            elif arg.is_negative:
                return -self(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return S.ImaginaryUnit * Basic.sinh(i_coeff)
            else:
                pi_coeff = arg.as_coefficient(S.Pi)

                if pi_coeff is not None:
                    if pi_coeff.is_integer:
                        return S.Zero
                    elif isinstance(pi_coeff, Basic.Rational):
                        cst_table = {
                            2 : S.One,
                            3 : S.Half*Basic.sqrt(3),
                            4 : S.Half*Basic.sqrt(2),
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


    @classmethod
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
                return (-1)**(n//2) * x**(n)/Basic.Factorial(n)

    def _eval_rewrite_as_exp(self, arg):
        exp, I = Basic.exp, S.ImaginaryUnit
        return (exp(arg*I) - exp(-arg*I)) / (2*I)

    def _eval_rewrite_as_cos(self, arg):
        return -cos(arg + S.Pi/2)

    def _eval_rewrite_as_tan(self, arg):
        tan_half = tan(S.Half*arg)
        return 2*tan_half/(1 + tan_half**2)

    def _eval_rewrite_as_cot(self, arg):
        cot_half = S.Cot(S.Half*arg)
        return 2*cot_half/(1 + cot_half**2)

    def _eval_conjugate(self):
        return self.func(self[0].conjugate())

    def _eval_expand_complex(self, *args):
        if self[0].is_real:
            return self
        re, im = self[0].as_real_imag()
        return sin(re)*Basic.cosh(im) + S.ImaginaryUnit*cos(re)*Basic.sinh(im)

    def _eval_expand_trig(self, *args):
        arg = self[0].expand()
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
        arg = self[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_is_real(self):
        return self[0].is_real

    def _eval_is_bounded(self):
        arg = self[0]
        if arg.is_real:
            return True

class cos(SingleValuedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -sin(self[0])
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return acos

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
                return S.One
            elif arg.is_negative:
                return self(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return Basic.cosh(i_coeff)
            else:
                pi_coeff = arg.as_coefficient(S.Pi)

                if pi_coeff is not None:
                    if isinstance(pi_coeff, Basic.Rational):
                        cst_table = {
                            1 : S.One,
                            2 : S.Zero,
                            3 : S.Half,
                            4 : S.Half*Basic.sqrt(2),
                            6 : S.Half*Basic.sqrt(3),
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


    @classmethod
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
                return (-1)**(n//2)*x**(n)/Basic.Factorial(n)

    def _eval_rewrite_as_exp(self, arg):
        exp, I = Basic.exp, S.ImaginaryUnit
        return (exp(arg*I) + exp(-arg*I)) / 2

    def _eval_rewrite_as_sin(self, arg):
        return sin(arg + S.Pi/2)

    def _eval_rewrite_as_tan(self, arg):
        tan_half = tan(S.Half*arg)**2
        return (1-tan_half)/(1+tan_half)

    def _eval_rewrite_as_cot(self, arg):
        cot_half = S.Cot(S.Half*arg)**2
        return (cot_half-1)/(cot_half+1)

    def _eval_conjugate(self):
        return self.func(self[0].conjugate())

    def _eval_expand_complex(self, *args):
        if self[0].is_real:
            return self
        re, im = self[0].as_real_imag()
        return cos(re)*Basic.cosh(im) - \
            S.ImaginaryUnit*sin(re)*Basic.sinh(im)

    def _eval_expand_trig(self, *args):
        arg = self[0].expand()
        x = None
        if isinstance(arg, Basic.Add):
            x = arg[0]
            y = Basic.Add(*arg[1:])
            return (cos(x)*cos(y) - sin(y)*sin(x)).expand(trig=True)
        else:
            coeff, terms = arg.as_coeff_terms()
            if not isinstance(coeff, Basic.One) and isinstance(coeff, Basic.Integer) and terms:
                x = Basic.Mul(*terms)
                return Basic.chebyshevt(coeff, cos(x))
        return cos(arg)

    def _eval_as_leading_term(self, x):
        arg = self[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return S.One
        else:
            return self.func(arg)

    def _eval_is_bounded(self):
        arg = self[0]

        if arg.is_real:
            return True

class tan(SingleValuedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex==1:
            return S.One + self**2
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return atan

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
            elif arg.is_negative:
                return -self(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return S.ImaginaryUnit * Basic.tanh(i_coeff)
            else:
                pi_coeff = arg.as_coefficient(S.Pi)

                if pi_coeff is not None:
                    if pi_coeff.is_integer:
                        return S.Zero
                    elif isinstance(pi_coeff, Basic.Rational):
                        cst_table = {
                           #2 : S.ComplexInfinity,
                            3 : Basic.sqrt(3),
                            4 : S.One,
                            6 : 1 / Basic.sqrt(3),
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


    @classmethod
    @cache_it_immutable
    def taylor_term(self, n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = Basic.sympify(x)

            a, b = ((n-1)//2), 2**(n+1)

            B = Basic.bernoulli(n+1)
            F = Basic.Factorial(n+1)

            return (-1)**a * b*(b-1) * B/F * x**n

    def _eval_conjugate(self):
        return self.func(self[0].conjugate())

    def _eval_expand_complex(self, *args):
        if self[0].is_real:
            return self
        re, im = self[0].as_real_imag()
        denom = cos(re)**2 + Basic.sinh(im)**2
        return (sin(re)*cos(re) + \
            S.ImaginaryUnit*Basic.sinh(im)*Basic.cosh(im))/denom

    def _eval_expand_trig(self, *args):
        return self

    def _eval_rewrite_as_exp(self, arg):
        exp, I = Basic.exp, S.ImaginaryUnit
        neg_exp, pos_exp = exp(-arg*I), exp(arg*I)
        return I*(neg_exp-pos_exp)/(neg_exp+pos_exp)

    def _eval_rewrite_as_sin(self, arg):
        return 2*sin(x)**2/sin(2*x)

    def _eval_rewrite_as_cos(self, arg):
        return -cos(x + S.Pi/2)/cos(x)

    def _eval_rewrite_as_cot(self, arg):
        return 1/S.Cot(arg)

    def _eval_as_leading_term(self, x):
        arg = self[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return S.One
        else:
            return self.func(arg)

    def _eval_is_bounded(self):
        arg = self[0]

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
                return -S.ImaginaryUnit * Basic.coth(i_coeff)
            else:
                pi_coeff = arg.as_coefficient(S.Pi)

                if pi_coeff is not None:
                    #if pi_coeff.is_integer:
                    #    return S.ComplexInfinity
                    if isinstance(pi_coeff, Basic.Rational):
                        cst_table = {
                            2 : S.Zero,
                            3 : 1 / Basic.sqrt(3),
                            4 : S.One,
                            6 : Basic.sqrt(3)
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


    @cache_it_immutable
    def taylor_term(self, n, x, *previous_terms):
        if n == 0:
            return 1 / Basic.sympify(x)
        elif n < 0 or n % 2 == 0:
            return S.Zero
        else:
            x = Basic.sympify(x)

            B = S.Bernoulli(n+1)
            F = Basic.Factorial(n+1)

            return (-1)**((n+1)//2) * 2**(n+1) * B/F * x**n

    def _eval_conjugate(self):
        assert len(self) == 1
        return self.func(self[0].conjugate())

    def _eval_expand_complex(self, *args):
        if self[0].is_real:
            return self
        re, im = self[0].as_real_imag()
        denom = sin(re)**2 + Basic.sinh(im)**2
        return (sin(re)*cos(re) - \
            S.ImaginaryUnit*Basic.sinh(im)*Basic.cosh(im))/denom

    def _eval_rewrite_as_exp(self, arg):
        exp, I = Basic.exp, S.ImaginaryUnit
        neg_exp, pos_exp = exp(-arg*I), exp(arg*I)
        return I*(pos_exp+neg_exp)/(pos_exp-neg_exp)

    def _eval_rewrite_as_sin(self, arg):
        return 2*sin(2*x)/sin(x)**2

    def _eval_rewrite_as_cos(self, arg):
        return -cos(x)/cos(x + S.Pi/2)

    def _eval_rewrite_as_tan(self, arg):
        return 1/tan(arg)

    def _eval_as_leading_term(self, x):
        arg = self[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return S.One
        else:
            return self.func(arg)

###############################################################################
########################### TRIGONOMETRIC INVERSES ############################
###############################################################################

class asin(SingleValuedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return (1 - self[0]**2)**(-S.Half)
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
                    Basic.sqrt(2)/2  : 4,
                    -Basic.sqrt(2)/2 : -4,
                    1/Basic.sqrt(2)  : 4,
                    -1/Basic.sqrt(2) : -4,
                    Basic.sqrt(3)/2  : 3,
                    -Basic.sqrt(3)/2 : -3,
                }

                if arg in cst_table:
                    return S.Pi / cst_table[arg]
                elif arg.is_negative:
                    return -self(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return S.ImaginaryUnit * Basic.asinh(i_coeff)
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
                return p * (n-2)**2/(k*(k-1)) * x**2
            else:
                k = (n - 1) // 2

                R = Basic.RisingFactorial(S.Half, k)
                F = Basic.Factorial(k)

                return R / F * x**n / n

    def _eval_as_leading_term(self, x):
        arg = self[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

class acos(SingleValuedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -(1 - self[0]**2)**(-S.Half)
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
                return S.Pi / 2
            elif isinstance(arg, Basic.One):
                return S.Zero
            elif isinstance(arg, Basic.NegativeOne):
                return S.Pi
            else:
                cst_table = {
                    S.Half       : S.Pi/3,
                    -S.Half      : 2*S.Pi/3,
                    Basic.sqrt(2)/2  : S.Pi/4,
                    -Basic.sqrt(2)/2 : 3*S.Pi/4,
                    1/Basic.sqrt(2)  : S.Pi/4,
                    -1/Basic.sqrt(2) : 3*S.Pi/4,
                    Basic.sqrt(3)/2  : S.Pi/6,
                    -Basic.sqrt(3)/2 : 5*S.Pi/6,
                }

                if arg in cst_table:
                    return cst_table[arg]


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

                R = Basic.RisingFactorial(S.Half, k)
                F = Basic.Factorial(k)

                return -R / F * x**n / n

    def _eval_as_leading_term(self, x):
        arg = self[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

class atan(SingleValuedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/(1+self[0]**2)
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
                    Basic.sqrt(3)/3  : 6,
                    -Basic.sqrt(3)/3 : -6,
                    1/Basic.sqrt(3)  : 6,
                    -1/Basic.sqrt(3) : -6,
                    Basic.sqrt(3)    : 3,
                    -Basic.sqrt(3)   : -3,
                }

                if arg in cst_table:
                    return S.Pi / cst_table[arg]
                elif arg.is_negative:
                    return -self(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return S.ImaginaryUnit * Basic.atanh(i_coeff)
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
            return (-1)**((n-1)//2) * x**n / n

    def _eval_as_leading_term(self, x):
        arg = self[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

class acot(SingleValuedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -1 / (1+self[0]**2)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def _eval_apply(cls, arg):
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
                    Basic.sqrt(3)/3  : 3,
                    -Basic.sqrt(3)/3 : -3,
                    1/Basic.sqrt(3)  : 3,
                    -1/Basic.sqrt(3) : -3,
                    Basic.sqrt(3)    : 6,
                    -Basic.sqrt(3)   : -6,
                }

                if arg in cst_table:
                    return S.Pi / cst_table[arg]
                elif arg.is_negative:
                    return -cls(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return -S.ImaginaryUnit * Basic.acoth(i_coeff)
            else:
                coeff, terms = arg.as_coeff_terms()

                if coeff.is_negative:
                    return -cls(-arg)


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
        arg = self[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

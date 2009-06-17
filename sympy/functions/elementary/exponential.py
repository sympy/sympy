
from sympy.core.basic import Basic, S, C, sympify, Wild
from sympy.core.function import Lambda, Function, Function, expand_log
from sympy.core.cache import cacheit

from sympy.utilities.decorator import deprecated

class exp(Function):

    nargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return self
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return log

    @classmethod
    def _eval_apply_subs(self, *args):
        return

    @classmethod
    @deprecated
    def canonize(cls, arg):
        return cls.eval(arg)

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg is S.NaN:
                return S.NaN
            elif arg is S.Zero:
                return S.One
            elif arg is S.One:
                return S.Exp1
            elif arg is S.Infinity:
                return S.Infinity
            elif arg is S.NegativeInfinity:
                return S.Zero
        elif arg.func is log:
            return arg.args[0]
        elif arg.is_Mul:
            coeff = arg.as_coefficient(S.Pi*S.ImaginaryUnit)

            if coeff is not None:
                if (2*coeff).is_integer:
                    if coeff.is_even:
                        return S.One
                    elif coeff.is_odd:
                        return S.NegativeOne
                    elif (coeff + S.Half).is_even:
                        return -S.ImaginaryUnit
                    elif (coeff + S.Half).is_odd:
                        return S.ImaginaryUnit
            I = S.ImaginaryUnit
            oo = S.Infinity
            a = Wild("a", exclude=[I, oo])
            r = arg.match(I*a*oo)
            if r and r[a] != 0:
                return S.NaN

        if arg.is_Add:
            args = arg.args[:]
        else:
            args = [arg]

        included, excluded = [], []

        for arg in args:
            coeff, terms = arg.as_coeff_terms()

            if coeff is S.Infinity:
                excluded.append(coeff**C.Mul(*terms))
            else:
                coeffs, log_term = [coeff], None

                for term in terms:
                    if term.func is log:
                        if log_term is None:
                            log_term = term.args[0]
                        else:
                            log_term = None
                            break
                    elif term.is_comparable:
                        coeffs.append(term)
                    else:
                        log_term = None
                        break

                if log_term is not None:
                    excluded.append(log_term**C.Mul(*coeffs))
                else:
                    included.append(arg)

        if excluded:
            return C.Mul(*(excluded+[cls(C.Add(*included))]))

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n<0: return S.Zero
        if n==0: return S.One
        x = sympify(x)
        if previous_terms:
            p = previous_terms[-1]
            if p is not None:
                return p * x / n
        return x**n/C.Factorial()(n)

    def _eval_expand_complex(self, deep=True, **hints):
        re, im = self.args[0].as_real_imag()
        if deep:
            re = re.expand(deep, **hints)
            im = im.expand(deep, **hints)
        cos, sin = C.cos(im), C.sin(im)
        return exp(re) * cos + S.ImaginaryUnit * exp(re) * sin

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def as_base_exp(self):
        return S.Exp1, C.Mul(*self.args)

    def as_coeff_terms(self, x=None):
        arg = self.args[0]
        if x is not None:
            c,f = arg.as_coeff_factors(x)
            return self.func(c), tuple( self.func(a) for a in f )
        if arg.is_Add:
            return S.One, tuple( self.func(a) for a in arg.args )
        return S.One,(self,)

    def _eval_subs(self, old, new):
        if self==old: return new
        arg = self.args[0]
        o = old
        if old.is_Pow: # handle (exp(3*log(x))).subs(x**2, z) -> z**(3/2)
            old = exp(old.exp * log(old.base))
        if old.func is exp:
            # exp(a*expr) .subs( exp(b*expr), y )  ->  y ** (a/b)
            a, expr_terms = self.args[0].as_coeff_terms()
            b, expr_terms_= old .args[0].as_coeff_terms()

            if expr_terms == expr_terms_:
                return new ** (a/b)


            if arg.is_Add: # exp(2*x+a).subs(exp(3*x),y) -> y**(2/3) * exp(a)
                # exp(exp(x) + exp(x**2)).subs(exp(exp(x)), w) -> w * exp(exp(x**2))
                oarg = old.args[0]
                new_l = []
                old_al = []
                coeff2,terms2 = oarg.as_coeff_terms()
                for a in arg.args:
                    a = a._eval_subs(old, new)
                    coeff1,terms1 = a.as_coeff_terms()
                    if terms1==terms2:
                        new_l.append(new**(coeff1/coeff2))
                    else:
                        old_al.append(a._eval_subs(old, new))
                if new_l:
                    new_l.append(self.func(C.Add(*old_al)))
                    r = C.Mul(*new_l)
                    return r
        old = o
        return Function._eval_subs(self, old, new)

    def _eval_is_real(self):
        return self.args[0].is_real

    def _eval_is_positive(self):
        if self.args[0].is_real:
            return True

    def _eval_is_bounded(self):
        arg = self.args[0]
        if arg.is_unbounded:
            if arg.is_negative: return True
            if arg.is_positive: return False
        if arg.is_bounded:
            return True
        if arg.is_real:
            return False
    def _eval_is_zero(self):
        return (self.args[0] is S.NegativeInfinity)

    def _eval_power(b, e):
        """exp(b[0])**e -> exp(b[0]*e)"""
        return exp(b.args[0] * e)

    def _eval_lseries(self, x, x0):
        s = self.args[0]
        yield exp(s.subs(x, x0))
        from sympy import Symbol, Integral, Derivative
        t = Symbol("t", dummy=True)
        f = s.subs(x, t)
        g = Integral(exp(f) * Derivative(f, t), (t, x0, x)).lseries(x, x0)
        for term in g:
            yield term

    def _eval_nseries(self, x, x0, n):
        from sympy import limit, Symbol, oo, powsimp
        arg = self.args[0]
        arg_series = arg.nseries(x, x0, n)
        if arg_series.is_Order:
            return 1+arg_series
        arg0 = limit(arg_series, x, x0)
        if arg0 in [-oo, oo]:
            return self
        s = Symbol("s", dummy=True)
        exp_series = exp(s)._taylor(s, x0, n)
        r = exp(arg0)*exp_series.subs(s, arg_series-arg0)
        r = r.expand()
        return powsimp(r, deep=True, combine='exp')

    def _taylor(self, x, x0, n):
        l = []
        g = None
        for i in xrange(n):
            g = self.taylor_term(i, self.args[0], g)
            g = g.nseries(x, x0, n)
            l.append(g)
        return C.Add(*l) + C.Order(x**n, x)

    def _eval_as_leading_term(self, x):
        arg = self.args[0]
        if arg.is_Add:
            return C.Mul(*[exp(f).as_leading_term(x) for f in arg.args])
        arg = self.args[0].as_leading_term(x)
        if C.Order(1,x).contains(arg):
            return S.One
        return exp(arg)

    def _eval_expand_power_exp(self, deep=True, **hints):
        if deep:
            arg = self.args[0].expand(deep=deep, **hints)
        else:
            arg = self.args[0]
        if arg.is_Add:
            expr = 1
            for x in arg.args:
                if deep:
                    x = x.expand(deep=deep, **hints)
                expr *= self.func(x)
            return expr
        return self.func(arg)

    def _eval_rewrite_as_sin(self, arg):
        I = S.ImaginaryUnit
        return C.sin(I*arg+S.Pi/2) - I*C.sin(I*arg)

    def _eval_rewrite_as_cos(self, arg):
        I = S.ImaginaryUnit
        return C.cos(I*arg) + I*C.cos(I*arg+S.Pi/2)

    def _sage_(self):
        import sage.all as sage
        return sage.exp(self.args[0]._sage_())

class log(Function):

    nargs = (1,2)

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/self.args[0]
            s = C.Symbol('x', dummy=True)
            return Lambda(s**(-1), s)
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return exp

    @classmethod
    def _eval_apply_subs(self, *args):
        return

    @classmethod
    @deprecated
    def canonize(cls, arg, base=None):
        return cls.eval(arg, base)

    @classmethod
    def eval(cls, arg, base=None):
        if base is not None:
            base = sympify(base)

            if base is not S.Exp1:
                return cls(arg)/cls(base)
            else:
                return cls(arg)

        arg = sympify(arg)

        if arg.is_Number:
            if arg is S.Zero:
                return S.NegativeInfinity
            elif arg is S.One:
                return S.Zero
            elif arg is S.Infinity:
                return S.Infinity
            elif arg is S.NegativeInfinity:
                return S.Infinity
            elif arg is S.NaN:
                return S.NaN
            elif arg.is_negative:
                return S.Pi * S.ImaginaryUnit + cls(-arg)
        elif arg is S.Exp1:
            return S.One
        #this doesn't work due to caching: :(
        #elif arg.func is exp and arg[0].is_real:
        #using this one instead:
        elif arg.func is exp:
            return arg.args[0]
        #this shouldn't happen automatically (see the issue 252):
        #elif arg.is_Pow:
        #    if arg.exp.is_Number or arg.exp.is_NumberSymbol or \
        #        arg.exp.is_number:
        #        return arg.exp * self(arg.base)
        #elif arg.is_Mul and arg.is_real:
        #    return C.Add(*[self(a) for a in arg])
        elif not arg.is_Add:
            coeff = arg.as_coefficient(S.ImaginaryUnit)

            if coeff is not None:
                if coeff is S.Infinity:
                    return S.Infinity
                elif coeff is S.NegativeInfinity:
                    return S.Infinity
                elif coeff.is_Rational:
                    if coeff.is_nonnegative:
                        return S.Pi * S.ImaginaryUnit * S.Half + cls(coeff)
                    else:
                        return -S.Pi * S.ImaginaryUnit * S.Half + cls(-coeff)

    def as_base_exp(self):
        return self, S.One
        #why is this here:?
        return exp, S.NegativeOne

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms): # of log(1+x)
        from sympy import powsimp
        if n<0: return S.Zero
        x = sympify(x)
        if n==0: return x
        if previous_terms:
            p = previous_terms[-1]
            if p is not None:
                return powsimp((-n) * p * x / (n+1), deep=True, combine='exp')
        return (1-2*(n%2)) * x**(n+1)/(n+1)

    def _eval_expand_log(self, deep=True, **hints):
        if deep:
            arg = self.args[0].expand(deep=deep, **hints)
        else:
            arg = self.args[0]
        if arg.is_Mul:
            expr = sympify(0)
            nonpos = sympify(1)
            for x in arg.args:
                if deep:
                    x = x.expand(deep=deep, **hints)
                if x.is_positive:
                    expr += self.func(x)._eval_expand_log(deep=deep, **hints)
                else:
                    nonpos *= x
            return expr + log(nonpos)
        elif arg.is_Pow:
            if arg.exp.is_real:# and arg.base.is_positive:
                # This should only run when base.is_positive, but it breaks
                # nseries, so it will have to wait for the new assumptions system.
                # See the variable obj2 in log._eval_nseries.
                if deep:
                    b = arg.base.expand(deep=deep, **hints)
                    e = arg.exp.expand(deep=deep, **hints)
                else:
                    b = arg.base
                    e = arg.exp
                return e * self.func(b)._eval_expand_log(deep=deep,\
                **hints)
        return self.func(arg)

    def _eval_expand_complex(self, deep=True, **hints):
        if deep:
            abs = C.abs(self.args[0].expand(deep, **hints))
            arg = C.arg(self.args[0].expand(deep, **hints))
        else:
            abs = C.abs(self.args[0])
            arg = C.arg(self.args[0])
        if hints['log']: # Expand the log
            hints['complex'] = False
            return log(abs).expand(deep, **hints) + S.ImaginaryUnit * arg
        else:
            return log(abs) + S.ImaginaryUnit * arg

    def _eval_is_real(self):
        return self.args[0].is_positive

    def _eval_is_bounded(self):
        arg = self.args[0]
        if arg.is_infinitesimal:
            return False
        return arg.is_bounded

    def _eval_is_positive(self):
        arg = self.args[0]
        if arg.is_positive:
            if arg.is_unbounded: return True
            if arg.is_infinitesimal: return False
            if arg.is_Number:
                return arg>1

    def _eval_is_zero(self):
        # XXX This is not quite useless. Try evaluating log(0.5).is_negative
        #     without it. There's probably a nicer way though.
        return (self.args[0] is S.One)

    def as_numer_denom(self):
        n, d = self.args[0].as_numer_denom()
        if d is S.One:
            return self.func(n), d
        return (self.func(n) - self.func(d)).as_numer_denom()

    def _eval_nseries(self, x, x0, n):
        from sympy import powsimp
        arg = self.args[0]
        k, l = Wild("k"), Wild("l")
        r = arg.match(k*x**l)
        if r is not None:
            k, l = r[k], r[l]
            if l != 0 and not l.has(x) and not k.has(x):
                r = log(k) + l*log(x)
                return r
        order = C.Order(x**n, x)
        arg = self.args[0]
        x = order.symbols[0]
        ln = C.log
        use_lt = not C.Order(1,x).contains(arg)
        if not use_lt:
            arg0 = arg.limit(x, 0)
            use_lt = (arg0 is S.Zero)
        if use_lt: # singularity, #example: self = log(sin(x))
            # arg = (arg / lt) * lt
            lt = arg.as_leading_term(x) # arg = sin(x); lt = x
            a = (arg/lt).expand() # a = sin(x)/x
            # the idea is to recursively call ln(a).series(), but one needs to
            # make sure that ln(sin(x)/x) doesn't get "simplified" to
            # -log(x)+ln(sin(x)) and an infinite recursion occurs, see also the
            # issue 252.
            obj = ln(lt) + ln(a)._eval_nseries(x, x0, n)
        else:
            # arg -> arg0 + (arg - arg0) -> arg0 * (1 + (arg/arg0 - 1))
            z = (arg/arg0 - 1)
            x = order.symbols[0]
            ln = C.log
            o = C.Order(z, x)
            if o is S.Zero:
                return ln(1+z)+ ln(arg0)
            if o.expr.is_number:
                e = ln(order.expr*x)/ln(x)
            else:
                e = ln(order.expr)/ln(o.expr)
            n = e.limit(x,0) + 1
            if n.is_unbounded:
                # requested accuracy gives infinite series,
                # order is probably nonpolynomial e.g. O(exp(-1/x), x).
                return ln(1+z)+ ln(arg0)
            try:
                n = int(n)
            except TypeError:
                #well, the n is something more complicated (like 1+log(2))
                n = int(n.evalf()) + 1
            assert n>=0,`n`
            l = []
            g = None
            for i in xrange(n+2):
                g = ln.taylor_term(i, z, g)
                g = g.nseries(x, x0, n)
                l.append(g)
            obj = C.Add(*l) + ln(arg0)
        obj2 = expand_log(powsimp(obj, deep=True, combine='exp'))
        if obj2 != obj:
            r = obj2.nseries(x, x0, n)
        else:
            r = obj
        if r == self:
            return self
        return r + order


    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)
        if arg is S.One:
            return (self.args[0] - 1).as_leading_term(x)
        return self.func(arg)

    #this is a lot faster:
    @classmethod
    def _eval_apply_evalf(cls, arg):
        arg = arg.evalf()
        if arg.is_number:
            import math
            from sympy import Real
            return Real(math.log(arg))

    def _sage_(self):
        import sage.all as sage
        return sage.log(self.args[0]._sage_())

# MrvLog is used by limit.py
class MrvLog(log):

    def _eval_subs(self, old, new):
        old = sympify(old)
        if old==self.func:
            arg = self.args[0]
            new = sympify(new)
            return new(arg._eval_subs(old, new))
        return self


class LambertW(Function):
    """Lambert W function, defined as the inverse function of
    x*exp(x). This function represents the principal branch
    of this inverse function, which like the natural logarithm
    is multivalued.

    For more information, see:
    http://en.wikipedia.org/wiki/Lambert_W_function
    """
    nargs = 1

    @classmethod
    @deprecated
    def canonize(cls, x):
        return cls.eval(x)

    @classmethod
    def eval(cls, x):
        if x == S.Zero: return S.Zero
        if x == S.Exp1: return S.One
        if x == -1/S.Exp1: return -S.One
        if x == -log(2)/2: return -log(2)
        if x == S.Infinity: return S.Infinity

    def fdiff(self, argindex=1):
        if argindex == 1:
            x = self.args[0]
            return LambertW(x)/(x*(1+LambertW(x)))
        else:
            raise ArgumentIndexError(self, argindex)

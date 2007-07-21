
from basic import Basic, S, cache_it, cache_it_immutable

from function import DefinedFunction, Apply, Lambda

class Exp(DefinedFunction):
    """ Exp() -> exp
    """
    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex==1:
            return self
        raise TypeError("argindex=%s is out of range [1,1] for %s" % (argindex,self))

    def inverse(self, argindex=1):
        return Log()

    def _eval_apply(self, arg):
        arg = Basic.sympify(arg)
        if isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.Zero):
                return S.One
            if isinstance(arg, Basic.One):
                return S.Exp1
            if isinstance(arg, Basic.Infinity):
                return S.Infinity
            if isinstance(arg, Basic.NegativeInfinity):
                return S.Zero
            if isinstance(arg, Basic.NaN):
                return S.NaN
        elif isinstance(arg, ApplyLog):
            return arg.args[0]
        elif isinstance(arg, (Basic.Add, Basic.Mul)):
            if isinstance(arg, Basic.Add):
                args = arg[:]
            else:
                args = [arg]
            l = []
            al = []
            for f in args:
                coeff, terms = f.as_coeff_terms()
                if isinstance(coeff, Basic.Infinity):
                    l.append(coeff**Basic.Mul(*terms))
                elif len(terms)==1 and isinstance(terms[0], Basic.ApplyLog):
                    l.append(terms[0].args[0]**coeff)
                else:
                    lt = None
                    cl = [coeff]
                    for t in terms:
                        if isinstance(t, Basic.ApplyLog):
                            if lt is None:
                                lt = t
                            else:
                                lt = None
                                break
                        elif t.is_comparable:
                            cl.append(t)
                        else:
                            lt = None
                            break
                    if lt is not None:
                        l.append(lt.args[0] ** Basic.Mul(*cl))
                    else:
                        al.append(f)
            if l:
                return Basic.Mul(*(l+[self(Basic.Add(*al))]))
                
    def _eval_apply_evalf(self, arg):
        arg = arg.evalf()
        if isinstance(arg, Basic.Number):
            return arg.exp()

    @cache_it_immutable
    def taylor_term(self, n, x, *previous_terms):
        if n<0: return S.Zero
        if n==0: return S.One
        x = Basic.sympify(x)
        if previous_terms:
            p = previous_terms[-1]
            if p is not None:
                return p * x / n
        return x**n/S.Factorial(n)

class ApplyExp(Apply):

    def precedence(self):
        b, e = self.as_base_exp()
        if e.is_negative: return 50 # same as default Mul
        return 70

    def tostr(self, level=0):
        p = self.precedence
        b, e = self.as_base_exp()
        if e.is_negative:
            r = '1/%s(%s)' % (self.func, -self.args[0])
        else:
            r = '%s(%s)' % (self.func, self.args[0])
        if p <= level:
            return '(%s)' % (r)
        return r

    def as_base_exp(self):
        #return Basic.Exp1(), self.args[0]
        coeff, terms = self.args[0].as_coeff_terms()
        return self.func(Basic.Mul(*terms)), coeff

    def as_coeff_terms(self, x=None):
        arg = self.args[0]
        if x is not None:
            c,f = arg.as_coeff_factors(x)
            return self.func(c), [self.func(a) for a in f]
        if isinstance(arg, Basic.Add):
            return Basic.One(), [self.func(a) for a in arg]
        return S.One,[self]

    def _eval_subs(self, old, new):
        if self==old: return new
        arg = self.args[0]
        o = old
        if isinstance(old, Basic.Pow): # handle (exp(3*log(x))).subs(x**2, z) -> z**(3/2)
            old = S.Exp(old.exp * S.Log(old.base))
        if isinstance(old, ApplyExp):
            b,e = self.as_base_exp()
            bo,eo = old.as_base_exp()
            if b==bo:
                return new ** (e/eo) # exp(2/3*x*3).subs(exp(3*x),y) -> y**(2/3)
            if isinstance(arg, Basic.Add): # exp(2*x+a).subs(exp(3*x),y) -> y**(2/3) * exp(a)
                # exp(exp(x) + exp(x**2)).subs(exp(exp(x)), w) -> w * exp(exp(x**2))
                oarg = old.args[0]
                new_l = []
                old_al = []
                coeff2,terms2 = oarg.as_coeff_terms()
                for a in arg:
                    a = a.subs(old, new)
                    coeff1,terms1 = a.as_coeff_terms()
                    if terms1==terms2:
                        new_l.append(new**(coeff1/coeff2))
                    else:
                        old_al.append(a.subs(old, new))
                if new_l:
                    new_l.append(self.func(Basic.Add(*old_al)))
                    r = Basic.Mul(*new_l)
                    return r
        old = o
        return self.func(arg.subs(old, new))

    def _eval_is_real(self):
        return self.args[0].is_real
    def _eval_is_positive(self):
        return self.args[0].is_real
    def _eval_is_bounded(self):
        arg = self.args[0]
        if arg.is_unbounded:
            if arg.is_negative: return True
            if arg.is_positive: return False
        if arg.is_bounded:
            return True

    def _eval_power(b, e):
        return b.func(b.args[0] * e)

    def _eval_oseries(self, order):
        arg = self.args[0]
        x = order.symbols[0]
        if not Basic.Order(1,x).contains(arg): # singularity
            arg0 = arg.as_leading_term(x)
            d = (arg-arg0).limit(x, S.Zero)
            if not isinstance(d, Basic.Zero):
                return S.Exp(arg)
        else:
            arg0 = arg.limit(x, S.Zero)
        o = order * S.Exp(-arg0)
        return self._compute_oseries(arg-arg0, o, S.Exp.taylor_term, S.Exp) * S.Exp(arg0)

    def _eval_as_leading_term(self, x):
        arg = self.args[0]
        if isinstance(arg, Basic.Add):
            return Basic.Mul(*[S.Exp(f).as_leading_term(x) for f in arg])
        arg = self.args[0].as_leading_term(x)
        if Basic.Order(1,x).contains(arg):
            return S.One
        return S.Exp(arg)


class Log(DefinedFunction):
    """ Log() -> log
    """
    is_comparable = True
    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex==1:
            s = Basic.Symbol('x',dummy=True)
            return Lambda(s**(-1),s)
        raise TypeError("argindex=%s is out of range [1,1] for %s" % (argindex,self))

    def inverse(self, argindex=1):
        return Exp()

    def _eval_apply(self, arg, base=None):
        if base is not None:
            base = Basic.sympify(base)
            if not isinstance(base, Basic.Exp1):
                return self(arg)/self(base)
        arg = Basic.sympify(arg)
        if isinstance(arg, Basic.Exp1):
            return S.One
        elif isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.One):
                return S.Zero
            if isinstance(arg, Basic.Infinity):
                return S.Infinity
            if isinstance(arg, Basic.NaN):
                return S.NaN
            if arg.is_negative:
                return S.Pi * S.ImaginaryUnit + self(-arg)
        elif isinstance(arg, Basic.Pow) and isinstance(arg.exp, Basic.Number):
            return arg.exp * self(arg.base)
        elif isinstance(arg, ApplyExp) and arg.args[0].is_real:
            return arg.args[0]
        elif isinstance(arg, Basic.Mul) and arg.is_real:
            return Basic.Add(*[self(a) for a in arg])

    def as_base_exp(self):
        return Exp(),Basic.Integer(-1)

    def _eval_apply_evalf(self, arg):
        arg = arg.evalf()
        if isinstance(arg, Basic.Number):
            return arg.log()
    
    def _calc_apply_positive(self, x):
        if x.is_positive and x.is_unbounded: return True

    def _calc_apply_unbounded(self, x):
        return x.is_unbounded

    @cache_it_immutable
    def taylor_term(self, n, x, *previous_terms): # of log(1+x)
        if n<0: return Basic.Zero()
        x = Basic.sympify(x)
        if n==0: return x
        if previous_terms:
            p = previous_terms[-1]
            if p is not None:
                return (-n) * p * x / (n+1)
        return (1-2*(n%2)) * x**(n+1)/(n+1)

class ApplyLog(Apply):

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
            if isinstance(arg, Basic.Number):
                return arg>1

    def as_numer_denom(self):
        n, d = self.args[0].as_numer_denom()
        if isinstance(d, Basic.One):
            return self.func(n), d
        return (self.func(n) - self.func(d)).as_numer_denom()

    def _eval_oseries(self, order):
        arg = self.args[0]
        x = order.symbols[0]
        ln = Basic.Log()
        use_lt = not Basic.Order(1,x).contains(arg)
        if not use_lt:
            arg0 = arg.limit(x, 0)
            use_lt = isinstance(arg0, Basic.Zero)
        if use_lt: # singularity
            # arg = (arg / lt) * lt
            lt = arg.as_leading_term(x)
            a = (arg/lt).expand()
            return ln(lt) + ln(a).oseries(order)
        # arg -> arg0 + (arg - arg0) -> arg0 * (1 + (arg/arg0 - 1))
        z = (arg/arg0 - 1)
        return self._compute_oseries(z, order, ln.taylor_term, lambda z: ln(1+z)) + ln(arg0)

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)
        if isinstance(arg, Basic.One):
            return (self.args[0] - 1).as_leading_term(x)
        return self.func(arg)

# MrvLog is used by limit.py
class MrvLog(Log):
    pass

class ApplyMrvLog(ApplyLog):

    def subs(self, old, new):
        old = Basic.sympify(old)
        if old==self.func:
            arg = self.args[0]
            new = Basic.sympify(new)
            return new(arg.subs(old, new))
        return self
#

class Sqrt(DefinedFunction):

    nofargs = 1
    
    def fdiff(self, argindex=1):
        if argindex==1:
            s = Basic.Symbol('x',dummy=True)
            return Lambda(Basic.Half() * s**(-Basic.Half()),s)
        raise TypeError("argindex=%s is out of range [1,1] for %s" % (argindex,self))

    def inverse(self, argindex=1):
        s = Basic.Symbol('x',dummy=True)
        return Lambda(s**2, s)

    def _eval_apply(self, arg):
        if isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.NaN):
                return S.NaN
            if isinstance(arg, Basic.Infinity):
                return S.Infinity
            if isinstance(arg, Basic.NegativeInfinity):
                return S.ImaginaryUnit * S.Infinity
            if isinstance(arg, Basic.Rational):
                factors = arg.factors()
                sqrt_factors = {}
                eval_factors = {}
                n = Basic.One()
                for k,v in factors.items():
                    n *= Basic.Integer(k) ** (v//2)
                    if v % 2:
                        n *= Basic.Integer(k) ** Basic.Half()
                return n
            return arg ** Basic.Half()
        coeff, terms = arg.as_coeff_terms()
        if not isinstance(coeff, Basic.One):
            return self(coeff) * self(Basic.Mul(*terms))
        base, exp = arg.as_base_exp()
        if isinstance(exp, Basic.Number):
            return base ** (exp/2)
        
    def _eval_apply_power(self, arg, exp):
        if isinstance(exp, Basic.Number):
            return arg ** (exp/2)

    def _eval_apply_evalf(self, arg):
        arg = arg.evalf()
        if isinstance(arg, Basic.Number):
            return arg.sqrt()

    def _eval_apply_subs(self, x, old, new):
        base, exp = old.as_base_exp()
        if base==x:
            return new ** (exp/2)

class ApplySqrt(Apply):

    def as_base_exp(self):
        return self.args[0], Basic.Half()

    def _eval_subs(self, old, new):
        if self==old: return new
        arg = self.args[0]
        func = self.func
        return func(arg.subs(old, new))


class Abs(DefinedFunction):

    nofargs = 1
    
    def fdiff(self, argindex=1):
        if argindex==1:
            raise NotImplementedError("Abs.fdiff()")
        raise TypeError("argindex=%s is out of range [1,1] for %s" % (argindex,self))

    def _eval_apply(self, arg):
        if isinstance(arg, Basic.NaN):
            return S.NaN
        if arg.is_positive: return arg
        if arg.is_negative: return -arg
        coeff, terms = arg.as_coeff_terms()
        if not isinstance(coeff, Basic.One):
            return self(coeff) * self(Basic.Mul(*terms))
        return

def Pi_coeff(expr):
    pi = Basic.Pi()
    if not expr.has(pi):
        return None
    x = Basic.Symbol('x',dummy=True)
    c = expr.subs(pi, x).diff(x)
    if c * pi == expr:
        return c
    return
 
class Sin(DefinedFunction):
    
    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex==1:
            return Cos()
        raise TypeError("argindex=%s is out of range [1,1] for %s" % (argindex,self))

    def inverse(self, argindex=1):
        return ASin()

    def _eval_apply(self, arg):
        if isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.NaN):
                return S.NaN
            if isinstance(arg, Basic.Zero):
                return arg
        c = Pi_coeff(arg)
        if c is not None:
            if isinstance(c, Basic.Integer):
                return Basic.Zero()
            if isinstance(c, Basic.Rational):
                c2 = 2 * c
                if isinstance(c2, Basic.Integer):
                    if (c2//2).is_even:
                        return Basic.One()
                    return Basic.NegativeOne()
        coeff, terms = arg.as_coeff_terms()
        if coeff.is_negative:
            return -self(-arg)
        return

    @cache_it_immutable
    def taylor_term(self, n, x, *previous_terms):
        if n<0: return Basic.Zero()
        x = Basic.sympify(x)
        if len(previous_terms)>1:
            p = previous_terms[-1]
            return -p * x**2 / ((2*n+1)*2*n)
        return (1-2*(n%2))*x**(2*n+1)/Basic.Factorial()(2*n+1)


class ApplySin(Apply):

    def _eval_expand(self):
        arg = self.args[0].expand()
        cos = Basic.Cos()
        sin = Basic.Sin()
        x = None
        if isinstance(arg, Basic.Add):
            x = arg[0]
            y = Basic.Add(*arg[1:])
        else:
            coeff, terms = arg.as_coeff_terms()
            if not isinstance(coeff, Basic.One) and isinstance(coeff, Basic.Integer):
                x = Basic.Mul(*terms)
                y = (coeff-1)*x
        if x is not None:
            return (sin(x)*cos(y) + sin(y)*cos(x)).expand()
        return sin(arg)

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)
        if Basic.Order(1,x).contains(arg):
            return arg
        return self.func(arg)

    def _eval_is_bounded(self):
        arg = self.args[0]
        if arg.is_real: return True

class Cos(DefinedFunction):
    
    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex==1:
            return -Sin()
        raise TypeError("argindex=%s is out of range [1,1] for %s" % (argindex,self))

    def inverse(self, argindex=1):
        return ACos()

    def _eval_apply(self, arg):
        if isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.NaN):
                return S.NaN
            if isinstance(arg, Basic.Zero):
                return Basic.One()
        c = Pi_coeff(arg)
        if c is not None:
            if isinstance(c, Basic.Integer):
                if c.is_even:
                    return Basic.One()
                return Basic.NegativeOne()
            if isinstance(c, Basic.Rational):
                c2 = 2 * c
                if isinstance(c2, Basic.Integer):
                    return Basic.Zero()
        coeff, terms = arg.as_coeff_terms()
        if coeff.is_negative:
            return self(-arg)
        return

    @cache_it_immutable
    def taylor_term(self, n, x, *previous_terms):
        if n<0: return Basic.Zero()
        x = Basic.sympify(x)
        if len(previous_terms)>1:
            p = previous_terms[-1]
            return -p * x**2 / ((2*n)*(2*n-1))
        return (1-2*(n%2))*x**(2*n)/Basic.Factorial()(2*n)

class ApplyCos(Apply):

    def _eval_expand(self):
        arg = self.args[0].expand()
        cos = Basic.Cos()
        sin = Basic.Sin()
        x = None
        if isinstance(arg, Basic.Add):
            x = arg[0]
            y = Basic.Add(*arg[1:])
            return (cos(x)*cos(y)-sin(y)*sin(x)).expand()
        else:
            coeff, terms = arg.as_coeff_terms()
            if not isinstance(coeff, Basic.One) and isinstance(coeff, Basic.Integer):
                x = Basic.Mul(*terms)
                return Basic.Chebyshev(coeff)(cos(x)).expand()
        return cos(arg)

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)
        if Basic.Order(1,x).contains(arg):
            return Basic.One()
        return self.func(arg)

    def _eval_is_bounded(self):
        arg = self.args[0]
        if arg.is_real: return True

class Tan(DefinedFunction):
    
    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex==1:
            return 1+self**2
        raise TypeError("argindex=%s is out of range [1,1] for %s" % (argindex,self))

    def inverse(self, argindex=1):
        return ATan()

    def _eval_apply(self, arg):
        if isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.NaN):
                return S.NaN
            if isinstance(arg, Basic.Zero):
                return arg
        return

    def _eval_apply_leadterm(self, x, arg):
        raise
        c0, e0 = arg.leadterm(x)
        if isinstance(e0, Basic.Zero):
            # tan(5+x) -> tan(5)
            c0 = self(c0)
            if not isinstance(c0, Basic.Zero):
                return c0, e0
            # tan(Pi+x) -> tan(-x)
            raise NotImplementedError("compute leading term %s(%s) at %s=0" % (self, arg, x))
        if e0.is_positive:
            # tan(2*x) -> 2 * x
            return c0, e0
        # tan(1/x)
        raise ValueError("unable to compute leading term %s(%s) at %s=0" % (self, arg, x))

class Sign(DefinedFunction):

    nofargs = 1

    def _eval_apply(self, arg):
        if isinstance(arg, Basic.NaN):
            return S.NaN
        if isinstance(arg, Basic.Zero): return Basic.One()
        if arg.is_positive: return Basic.One()
        if arg.is_negative: return -Basic.One()
        if isinstance(arg, Basic.Mul):
            coeff, terms = arg.as_coeff_terms()
            if not isinstance(coeff, Basic.One):
                return self(coeff) * self(Basic.Mul(*terms))

class ApplySign(Apply):

    is_bounded = True

class Conjugate(DefinedFunction):

    nofargs = 1

    def _eval_apply(self, arg):
        obj = arg._eval_conjugate()
        if obj is not None:
            return obj

class ApplyConjugate(Apply):

    def _eval_conjugate(self):
        return self.args[0]

class Max(DefinedFunction):

    nofargs = 2
    def _eval_apply(self, x, y):
        if isinstance(x, Basic.Number) and isinstance(y, Basic.Number):
            return max(x, y)
        if x.is_positive:
            if y.is_negative:
                return x
            if y.is_positive:
                if x.is_unbounded:
                    if y.is_unbounded:
                        return
                    return x
        elif x.is_negative:
            if y.is_negative:
                if y.is_unbounded:
                    if x.is_unbounded:
                        return
                    return x

class Min(DefinedFunction):

    nofargs = 2
    def _eval_apply(self, x, y):
        if isinstance(x, Basic.Number) and isinstance(y, Basic.Number):
            return min(x, y)


Basic.singleton['exp'] = Exp
Basic.singleton['log'] = Log
Basic.singleton['ln'] = Log
Basic.singleton['sin'] = Sin
Basic.singleton['cos'] = Cos
Basic.singleton['tan'] = Tan
Basic.singleton['sqrt'] = Sqrt
Basic.singleton['abs_'] = Abs
Basic.singleton['max_'] = Max
Basic.singleton['min_'] = Min
Basic.singleton['sign'] = Sign
Basic.singleton['conjugate'] = Conjugate

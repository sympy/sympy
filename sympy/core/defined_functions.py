
from basic import Basic, S, cache_it, cache_it_immutable
from function import DefinedFunction, Apply, Lambda
from symbol import Symbol

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

    def _eval_complex_expand(self):
        re, im = self.args[0].as_real_imag()

        exp, cos, sin = Exp()(re), Cos()(im), Sin()(im)
        return exp * cos + S.ImaginaryUnit * exp * sin

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
        if self.args[0].is_real:
            return True
    def _eval_is_bounded(self):
        arg = self.args[0]
        if arg.is_unbounded:
            if arg.is_negative: return True
            if arg.is_positive: return False
        if arg.is_bounded:
            return True
    def _eval_is_zero(self):
        return isinstance(self.args[0], Basic.NegativeInfinity)

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

    def _eval_expand(self):
        arg = self.args[0].expand()
        if isinstance(arg, Basic.Add):
            expr = 1
            for x in arg:
                expr *= self.func(x).expand()
            return expr
        return self.func(arg)


class Log(DefinedFunction):
    """ Log() -> log
    """
    is_comparable = True

    def fdiff(self, argindex=1):
        if argindex == 1:
            s = Basic.Symbol('x', dummy=True)
            return Lambda(s**(-1), s)
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return Exp()

    def _eval_apply(self, arg, base=None):
        if base is not None:
            base = Basic.sympify(base)

            if not isinstance(base, Basic.Exp1):
                return self(arg)/self(base)

        arg = Basic.sympify(arg)

        if isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.Zero):
                return S.NegativeInfinity
            elif isinstance(arg, Basic.One):
                return S.Zero
            elif isinstance(arg, Basic.Infinity):
                return S.Infinity
            elif isinstance(arg, Basic.NegativeInfinity):
                return S.Infinity
            elif isinstance(arg, Basic.NaN):
                return S.NaN
            elif arg.is_negative:
                return S.Pi * S.ImaginaryUnit + self(-arg)
        elif isinstance(arg, Basic.Exp1):
            return S.One
        elif isinstance(arg, ApplyExp) and arg.args[0].is_real:
            return arg.args[0]
        elif isinstance(arg, Basic.Pow):
            if isinstance(arg.exp, Basic.Number) or \
               isinstance(arg.exp, Basic.NumberSymbol):
                return arg.exp * self(arg.base)
        elif isinstance(arg, Basic.Mul) and arg.is_real:
            return Basic.Add(*[self(a) for a in arg])
        elif not isinstance(arg, Basic.Add):
            coeff = arg.as_coefficient(S.ImaginaryUnit)

            if coeff is not None:
                if isinstance(coeff, Basic.Infinity):
                    return S.Infinity
                elif isinstance(coeff, Basic.NegativeInfinity):
                    return S.Infinity
                elif isinstance(coeff, Basic.Rational):
                    if coeff.is_nonnegative:
                        return S.Pi * S.ImaginaryUnit * S.Half + self(coeff)
                    else:
                        return -S.Pi * S.ImaginaryUnit * S.Half + self(-coeff)

    def as_base_exp(self):
        return Exp(), Basic.Integer(-1)

    def _eval_apply_evalf(self, arg):
        arg = arg.evalf()
        if isinstance(arg, Basic.Number):
            return arg.log()

    def _calc_apply_positive(self, x):
        if x.is_positive and x.is_unbounded: return True

    def _calc_apply_unbounded(self, x):
        import pdb; pdb.set_trace()
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

    def _eval_complex_expand(self):
        re, im = self.args[0].as_real_imag()

        return Log()(Sqrt()(re) + Sqrt()(im)) + \
               S.ImaginaryUnit * Arg()(self.args[0])

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
    def _eval_is_zero(self): # USELESS
        return isinstance(self.args[0], Basic.One)

    def as_numer_denom(self):
        n, d = self.args[0].as_numer_denom()
        if isinstance(d, Basic.One):
            return self.func(n), d
        return (self.func(n) - self.func(d)).as_numer_denom()

    # similar code must be added to other functions with have singularites
    # in their domains eg. cot(), tan() ...
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

    def _eval_expand(self): # WRONG
        arg = self.args[0]
        if isinstance(arg, Basic.Mul):
            expr = 0
            for x in arg:
                expr += self.func(x).expand()
            return expr
        elif isinstance(arg, Basic.Pow):
            return arg.exp * self.func(arg.base).expand()
        return self

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
        if arg.is_nonnegative:
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
    def _eval_is_zero(self):
        return isinstance(self.args, Basic.Zero)

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

    def _eval_is_zero(self):
        return isinstance(self.args[0], Basic.Zero)

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

    def _eval_is_zero(self):
        return isinstance(self.args[0], Basic.Zero)

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
Basic.singleton['sqrt'] = Sqrt
Basic.singleton['abs_'] = Abs
Basic.singleton['max_'] = Max
Basic.singleton['min_'] = Min
Basic.singleton['sign'] = Sign

###############################################################################
######################### FLOOR and CEILING FUNCTIONS #########################
###############################################################################

class Floor(DefinedFunction):
    """Floor is a univariate function which returns the largest integer
       value not greater than its argument. However this implementaion
       generalizes floor to complex numbers.

       More information can be found in "Concrete mathematics" by Graham,
       pp. 87 or visit http://mathworld.wolfram.com/FloorFunction.html.

       >>> from sympy import *

       >>> floor(17)
       17

       >>> floor(Rational(23, 10))
       2

       >>> floor(2*E)
       5

       >>> floor(-Real(0.567))
       -1

       >>> floor(-I/2)
       -I

    """

    nofargs = 1

    def _eval_apply(self, arg):
        arg = Basic.sympify(arg)

        if arg.is_integer:
            return arg
        elif isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.Infinity):
                return S.Infinity
            elif isinstance(arg, Basic.NegativeInfinity):
                return S.NegativeInfinity
            elif isinstance(arg, Basic.NaN):
                return S.NaN
            elif isinstance(arg, Basic.Integer):
                return arg
            elif isinstance(arg, Basic.Rational):
                return Basic.Integer(arg.p // arg.q)
            elif isinstance(arg, Basic.Real):
                return Basic.Integer(int(arg.floor()))
        elif isinstance(arg, Basic.NumberSymbol):
            return arg.approximation_interval(Basic.Integer)[0]
        elif isinstance(arg, Basic.ImaginaryUnit):
            return S.ImaginaryUnit
        elif isinstance(arg, Basic.Add):
            included, excluded = [], []

            for term in arg:
                coeff = term.as_coefficient(S.ImaginaryUnit)

                if coeff is not None and coeff.is_real:
                    excluded.append(self(coeff)*S.ImaginaryUnit)
                elif term.is_real:
                    if term.is_integer:
                        excluded.append(term)
                    else:
                        included.append(term)
                else:
                    return

            if excluded:
                return self(Basic.Add(*included)) + Basic.Add(*excluded)
        else:
            coeff, terms = arg.as_coeff_terms(S.ImaginaryUnit)

            if not terms and not arg.atoms(type=Basic.Symbol):
                if arg.is_negative:
                    return -Ceiling()(-arg)
                else:
                    return self(arg.evalf())
            elif terms == [ S.ImaginaryUnit ] and coeff.is_real:
                return self(coeff)*S.ImaginaryUnit

class ApplyFloor(Apply):

    def _eval_is_bounded(self):
        return self.args[0].is_bounded

    def _eval_is_real(self):
        return self.args[0].is_real

    def _eval_is_integer(self):
        return self.args[0].is_real

class Ceiling(DefinedFunction):
    """Ceiling is a univariate function which returns the smallest integer
       value not less than its argument. Ceiling function is generalized
       in this implementation to complex numbers.

       More information can be found in "Concrete mathematics" by Graham,
       pp. 87 or visit http://mathworld.wolfram.com/CeilingFunction.html.

       >>> from sympy import *

       >>> ceiling(17)
       17

       >>> ceiling(Rational(23, 10))
       3

       >>> ceiling(2*E)
       6

       >>> ceiling(-Real(0.567))
       0

       >>> ceiling(I/2)
       I

    """

    nofargs = 1

    def _eval_apply(self, arg):
        arg = Basic.sympify(arg)

        if arg.is_integer:
            return arg
        elif isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.Infinity):
                return S.Infinity
            elif isinstance(arg, Basic.NegativeInfinity):
                return S.NegativeInfinity
            elif isinstance(arg, Basic.NaN):
                return S.NaN
            elif isinstance(arg, Basic.Integer):
                return arg
            elif isinstance(arg, Basic.Rational):
                return Basic.Integer(arg.p // arg.q + 1)
            elif isinstance(arg, Basic.Real):
                return Basic.Integer(int(arg.ceiling()))
        elif isinstance(arg, Basic.NumberSymbol):
            return arg.approximation_interval(Basic.Integer)[1]
        elif isinstance(arg, Basic.ImaginaryUnit):
            return S.ImaginaryUnit
        elif isinstance(arg, Basic.Add):
            included, excluded = [], []

            for term in arg:
                coeff = term.as_coefficient(S.ImaginaryUnit)

                if coeff is not None and coeff.is_real:
                    excluded.append(self(coeff)*S.ImaginaryUnit)
                elif term.is_real:
                    if term.is_integer:
                        excluded.append(term)
                    else:
                        included.append(term)
                else:
                    return

            if excluded:
                return self(Basic.Add(*included)) + Basic.Add(*excluded)
        else:
            coeff, terms = arg.as_coeff_terms(S.ImaginaryUnit)

            if not terms and not arg.atoms(type=Basic.Symbol):
                if arg.is_negative:
                    return -Floor()(-arg)
                else:
                    return self(arg.evalf())
            elif terms == [ S.ImaginaryUnit ] and coeff.is_real:
                return self(coeff)*S.ImaginaryUnit

class ApplyCeiling(Apply):

    def _eval_is_bounded(self):
        return self.args[0].is_bounded

    def _eval_is_real(self):
        return self.args[0].is_real

    def _eval_is_integer(self):
        return self.args[0].is_real

Basic.singleton['floor'] = Floor
Basic.singleton['ceiling'] = Ceiling

###############################################################################
######################## RISING and FALLING FACTORIALS ########################
###############################################################################

class RisingFactorial(DefinedFunction):
    """Rising factorial (also called Pochhammer symbol) is a double valued
       function arising in concrete mathematics, hypergeometric functions
       and series expanansions. It is defined by

                   rf(x, k) = x * (x+1) * ... * (x + k-1)

       where 'x' can be arbitrary expression and 'k' is an integer. For
       more information check "Concrete mathematics" by Graham, pp. 66
       or visit http://mathworld.wolfram.com/RisingFactorial.html page.

       >>> from sympy import *
       >>> x = Symbol('x')

       >>> rf(x, 0)
       1

       >>> rf(1, 5)
       120

       >>> rf(x, 5)
       x*(1 + x)*(2 + x)*(3 + x)*(4 + x)

    """

    nofargs = 2

    def _eval_apply(self, x, k):
        x = Basic.sympify(x)
        k = Basic.sympify(k)

        if isinstance(x, Basic.NaN):
            return Basic.NaN()
        elif isinstance(k, Basic.Integer):
            if isinstance(k, Basic.NaN):
                return Basic.NaN()
            if isinstance(k, Basic.Zero):
                return Basic.One()
            else:
                if k.is_positive:
                    if isinstance(x, Basic.Infinity):
                        return Basic.Infinity()
                    elif isinstance(x, Basic.NegativeInfinity):
                        if k.is_odd:
                            return Basic.NegativeInfinity()
                        else:
                            return Basic.Infinity()
                    else:
                        return reduce(lambda r, i: r*(x+i), xrange(0, int(k)), 1)
                else:
                    if isinstance(x, Basic.Infinity):
                        return Basic.Infinity()
                    elif isinstance(x, Basic.NegativeInfinity):
                        return Basic.Infinity()
                    else:
                        return 1/reduce(lambda r, i: r*(x-i), xrange(1, abs(int(k))+1), 1)

class ApplyRisingFactorial(Apply):

    def tostr(self, level=0):
        r = 'rf(%s)' % ', '.join([a.tostr() for a in self.args])

        if self.precedence <= level:
            return '(%s)' % (r)
        else:
            return r

class FallingFactorial(DefinedFunction):
    """Falling factorial (related to rising factorial) is a double valued
       function arising in concrete mathematics, hypergeometric functions
       and series expanansions. It is defined by

                   ff(x, k) = x * (x-1) * ... * (x - k+1)

       where 'x' can be arbitrary expression and 'k' is an integer. For
       more information check "Concrete mathematics" by Graham, pp. 66
       or visit http://mathworld.wolfram.com/FallingFactorial.html page.

       >>> from sympy import *
       >>> x = Symbol('x')

       >>> ff(x, 0)
       1

       >>> ff(5, 5)
       120

       >>> ff(x, 5)
       x*(1 - x)*(2 - x)*(3 - x)*(4 - x)

    """

    nofargs = 2

    def _eval_apply(self, x, k):
        x = Basic.sympify(x)
        k = Basic.sympify(k)

        if isinstance(x, Basic.NaN):
            return Basic.NaN()
        elif isinstance(k, Basic.Integer):
            if isinstance(k, Basic.NaN):
                return Basic.NaN()
            if isinstance(k, Basic.Zero):
                return Basic.One()
            else:
                result = Basic.One()

                if k.is_positive:
                    if isinstance(x, Basic.Infinity):
                        return Basic.Infinity()
                    elif isinstance(x, Basic.NegativeInfinity):
                        if k.is_odd:
                            return Basic.NegativeInfinity()
                        else:
                            return Basic.Infinity()
                    else:
                        return reduce(lambda r, i: r*(x-i), xrange(0, int(k)), 1)
                else:
                    if isinstance(x, Basic.Infinity):
                        return Basic.Infinity()
                    elif isinstance(x, Basic.NegativeInfinity):
                        return Basic.Infinity()
                    else:
                        return 1/reduce(lambda r, i: r*(x+i), xrange(1, abs(int(k))+1), 1)

class ApplyFallingFactorial(Apply):

    def tostr(self, level=0):
        r = 'ff(%s)' % ', '.join([a.tostr() for a in self.args])

        if self.precedence <= level:
            return '(%s)' % (r)
        else:
            return r

Basic.singleton['pochhammer'] = RisingFactorial

Basic.singleton['rf'] = RisingFactorial
Basic.singleton['ff'] = FallingFactorial

###############################################################################
############## REAL and IMAGINARY PARTS, ARGUMENT, CONJUGATION ################
###############################################################################

class Re(DefinedFunction):
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

       >>> re(x*I)
       -im(x)

       >>> re(im(x) + x*I + 2)
       2

    """

    nofargs = 1

    is_real = True

    def _eval_apply(self, arg):
        arg = Basic.sympify(arg)

        if isinstance(arg, Basic.NaN):
            return S.NaN
        elif arg.is_real:
            return arg
        else:
            if not isinstance(arg, Basic.Add):
                arg = [arg]

            included, reverted, excluded = [], [], []

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
                a, b, c = map(lambda xs: Basic.Add(*xs),
                    [included, reverted, excluded])

                return self(a) - Im()(b) + c

class ApplyRe(Apply):

    def _eval_is_real(self):
        return True

    def _eval_complex_expand(self):
        return self.func(self.args[0].as_real_imag()[0])

class Im(DefinedFunction):
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

    nofargs = 1

    is_real = True

    def _eval_apply(self, arg):
        arg = Basic.sympify(arg)

        if isinstance(arg, Basic.NaN):
            return S.NaN
        elif arg.is_real:
            return S.Zero
        else:
            if not isinstance(arg, Basic.Add):
                arg = [arg]

            included, reverted, excluded = [], [], []

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
                a, b, c = map(lambda xs: Basic.Add(*xs),
                    [included, reverted, excluded])

                return self(a) + Re()(b) + c

class ApplyIm(Apply):

    def _eval_is_real(self):
        return True

    def _eval_complex_expand(self):
        return self.func(self.args[0].as_real_imag()[1])

class Arg(DefinedFunction):

    nofargs = 1

    is_real = True

    def _eval_apply(self, arg):
        return

class ApplyArg(Apply):

    def _eval_is_real(self):
        return True

Basic.singleton['re'] = Re
Basic.singleton['im'] = Im
Basic.singleton['arg'] = Arg

class Conjugate(DefinedFunction):

    nofargs = 1

    def _eval_apply(self, arg):
        obj = arg._eval_conjugate()

        if obj is not None:
            return obj

class ApplyConjugate(Apply):

    def _eval_conjugate(self):
        return self.args[0]

Basic.singleton['conjugate'] = Conjugate

###############################################################################
########################## TRIGONOMETRIC FUNCTIONS ############################
###############################################################################

class Sin(DefinedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return Cos()
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return ASin()

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
                return S.ImaginaryUnit * Sinh()(i_coeff)
            else:
                pi_coeff = arg.as_coefficient(S.Pi)

                if pi_coeff is not None:
                    if pi_coeff.is_integer:
                        return S.Zero
                    elif isinstance(pi_coeff, Basic.Rational):
                        cst_table = {
                            2 : Basic.One(),
                            3 : Basic.Half()*Sqrt()(3),
                            4 : Basic.Half()*Sqrt()(2),
                            6 : Basic.Half(),
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
            return Basic.Zero()
        else:
            x = Basic.sympify(x)

            if len(previous_terms) > 2:
                p = previous_terms[-2]
                return -p * x**2 / (n*(n-1))
            else:
                return (-1)**(n//2) * x**(n)/Basic.Factorial()(n)

class ApplySin(Apply):

    def _eval_complex_expand(self):
        re, im = self.args[0].as_real_imag()

        return Sin()(re)*Cosh()(im) + \
            S.ImaginaryUnit*Cos()(re)*Sinh()(im)

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
            if not isinstance(coeff, Basic.One) and isinstance(coeff, Basic.Integer) and terms:
                x = Basic.Mul(*terms)
                y = (coeff-1)*x
        if x is not None:
            return (sin(x)*cos(y) + sin(y)*cos(x)).expand()
        return sin(arg)

    #def _eval_complex_expand(self):

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_is_bounded(self):
        arg = self.args[0]

        if arg.is_real:
            return True

class Cos(DefinedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -Sin()
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return ACos()

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
                return Cosh()(i_coeff)
            else:
                pi_coeff = arg.as_coefficient(S.Pi)

                if pi_coeff is not None:
                    if isinstance(pi_coeff, Basic.Rational):
                        cst_table = {
                            1 : Basic.One(),
                            2 : Basic.Zero(),
                            3 : Basic.Half(),
                            4 : Basic.Half()*Sqrt()(2),
                            6 : Basic.Half()*Sqrt()(3),
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
            return Basic.Zero()
        else:
            x = Basic.sympify(x)

            if len(previous_terms) > 2:
                p = previous_terms[-2]
                return -p * x**2 / (n*(n-1))
            else:
                return (-1)**(n//2)*x**(n)/Basic.Factorial()(n)

class ApplyCos(Apply):

    def _eval_complex_expand(self):
        re, im = self.args[0].as_real_imag()

        return Cos()(re)*Cosh()(im) - \
            S.ImaginaryUnit*Sin()(re)*Sinh()(im)

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
            if not isinstance(coeff, Basic.One) and isinstance(coeff, Basic.Integer) and terms:
                x = Basic.Mul(*terms)
                return Basic.Chebyshev(coeff)(cos(x)).expand()
        return cos(arg)

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return Basic.One()
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
            return 1/Cos()**2
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return ATan()

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
                return S.ImaginaryUnit * Tanh()(i_coeff)
            else:
                pi_coeff = arg.as_coefficient(S.Pi)

                if pi_coeff is not None:
                    if pi_coeff.is_integer:
                        return S.Zero
                    elif isinstance(pi_coeff, Basic.Rational):
                        cst_table = {
                            3 : Sqrt()(3),
                            4 : Basic.One(),
                            6 : 1 / Sqrt()(3),
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
            return Basic.Zero()
        else:
            x = Basic.sympify(x)

            a, b = ((n-1)//2), 2**(n+1)

            B = Basic.Bernoulli()(n+1)
            F = Basic.Factorial()(n+1)

            return (-1)**a * b*(b-1) * B/F * x**n

class ApplyTan(Apply):

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return Basic.One()
        else:
            return self.func(arg)

    def _eval_is_bounded(self):
        arg = self.args[0]

        if arg.is_imaginary:
            return True

class Cot(DefinedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -1/Sin()**2
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return ACot()

    def _eval_apply(self, arg):
        arg = Basic.sympify(arg)

        if isinstance(arg, Basic.Number):
            if isinstance(arg, Basic.NaN):
                return S.NaN
            elif arg.is_negative:
                return -self(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return -S.ImaginaryUnit * Coth()(i_coeff)
            else:
                pi_coeff = arg.as_coefficient(S.Pi)

                if pi_coeff is not None:
                    if isinstance(pi_coeff, Basic.Rational):
                        cst_table = {
                            2 : Basic.Zero(),
                            3 : 1 / Sqrt()(3),
                            4 : Basic.One(),
                            6 : Sqrt()(3)
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
            return Basic.Zero()
        else:
            x = Basic.sympify(x)

            B = Basic.Bernoulli()(n+1)
            F = Basic.Factorial()(n+1)

            return (-1)**((n+1)//2) * 2**(n+1) * B/F * x**n

class ApplyCot(Apply):

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return Basic.One()
        else:
            return self.func(arg)

Basic.singleton['sin'] = Sin
Basic.singleton['cos'] = Cos
Basic.singleton['tan'] = Tan
Basic.singleton['cot'] = Cot

###############################################################################
########################### HYPERBOLIC FUNCTIONS ##############################
###############################################################################

class Sinh(DefinedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return Cosh()
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return ASinh()

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
                return S.ImaginaryUnit * Sin()(i_coeff)
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
                return x**(n) / Basic.Factorial()(n)

class ApplySinh(Apply):

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return Basic.One()
        else:
            return self.func(arg)

    def _eval_is_bounded(self):
        arg = self.args[0]

        if arg.is_imaginary:
            return True

class Cosh(DefinedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return Sinh()
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return ACosh()

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
                return Cos()(i_coeff)
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
            return Basic.Zero()
        else:
            x = Basic.sympify(x)

            if len(previous_terms) > 2:
                p = previous_terms[-2]
                return p * x**2 / (n*(n-1))
            else:
                return x**(n)/Basic.Factorial()(n)

class ApplyCosh(Apply):

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return Basic.One()
        else:
            return self.func(arg)

    def _eval_is_bounded(self):
        arg = self.args[0]

        if arg.is_imaginary:
            return True

class Tanh(DefinedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/Cosh()**2
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return ATanh()

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
                return S.ImaginaryUnit * Tan()(i_coeff)
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
            return Basic.Zero()
        else:
            x = Basic.sympify(x)

            a = 2**(n+1)

            B = Basic.Bernoulli()(n+1)
            F = Basic.Factorial()(n+1)

            return a*(a-1) * B/F * x**n

class ApplyTanh(Apply):

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return Basic.One()
        else:
            return self.func(arg)

    def _eval_is_bounded(self):
        arg = self.args[0]

        if arg.is_real:
            return True

class Coth(DefinedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/Sinh()**2
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return ACoth()

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
                return -S.ImaginaryUnit * Cot()(i_coeff)
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
            return Basic.Zero()
        else:
            x = Basic.sympify(x)

            B = Basic.Bernoulli()(n+1)
            F = Basic.Factorial()(n+1)

            return 2**(n+1) * B/F * x**n

class ApplyCoth(Apply):

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return Basic.One()
        else:
            return self.func(arg)

Basic.singleton['sinh'] = Sinh
Basic.singleton['cosh'] = Cosh
Basic.singleton['tanh'] = Tanh
Basic.singleton['coth'] = Coth

###############################################################################
########################### TRIGONOMETRIC INVERSES ############################
###############################################################################

class ASin(DefinedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            s = Basic.Symbol('x', dummy=True)
            return Lambda((1 - s**2)**(-Basic.Half()), s)
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
                    Sqrt()(2)/2  : 4,
                    -Sqrt()(2)/2 : -4,
                    1/Sqrt()(2)  : 4,
                    -1/Sqrt()(2) : -4,
                    Sqrt()(3)/2  : 3,
                    -Sqrt()(3)/2 : -3,
                }

                if arg in cst_table:
                    return S.Pi / cst_table[arg]
                elif arg.is_negative:
                    return -self(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return S.ImaginaryUnit * ASinh()(i_coeff)
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
            return Basic.Zero()
        else:
            x = Basic.sympify(x)

            if len(previous_terms) > 2:
                p = previous_terms[-2]
                return p * (n-2)**2/(k*(k-1)) * x**2
            else:
                k = (n - 1) // 2

                R = Basic.RisingFactorial()(Basic.Half(), k)
                F = Basic.Factorial()(k)

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
            return Lambda(-(1 - s**2)**(-Basic.Half()), s)
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
                    Sqrt()(2)/2  : S.Pi/4,
                    -Sqrt()(2)/2 : 3*S.Pi/4,
                    1/Sqrt()(2)  : S.Pi/4,
                    -1/Sqrt()(2) : 3*S.Pi/4,
                    Sqrt()(3)/2  : S.Pi/6,
                    -Sqrt()(3)/2 : 5*S.Pi/6,
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
            return Basic.Zero()
        else:
            x = Basic.sympify(x)

            if len(previous_terms) > 2:
                p = previous_terms[-2]
                return p * (n-2)**2/(k*(k-1)) * x**2
            else:
                k = (n - 1) // 2

                R = Basic.RisingFactorial()(Basic.Half(), k)
                F = Basic.Factorial()(k)

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
                    Sqrt()(3)/3  : 6,
                    -Sqrt()(3)/3 : -6,
                    1/Sqrt()(3)  : 6,
                    -1/Sqrt()(3) : -6,
                    Sqrt()(3)    : 3,
                    -Sqrt()(3)   : -3,
                }

                if arg in cst_table:
                    return S.Pi / cst_table[arg]
                elif arg.is_negative:
                    return -self(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return S.ImaginaryUnit * ATanh()(i_coeff)
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
            return Basic.Zero()
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

class ACot(DefinedFunction):

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
                    Sqrt()(3)/3  : 3,
                    -Sqrt()(3)/3 : -3,
                    1/Sqrt()(3)  : 3,
                    -1/Sqrt()(3) : -3,
                    Sqrt()(3)    : 6,
                    -Sqrt()(3)   : -6,
                }

                if arg in cst_table:
                    return S.Pi / cst_table[arg]
                elif arg.is_negative:
                    return -self(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return -S.ImaginaryUnit * ACoth()(i_coeff)
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
            return Basic.Zero()
        else:
            x = Basic.sympify(x)
            return (-1)**((n+1)//2) * x**n / n

class ApplyACot(Apply):

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)

        if Basic.Order(1,x).contains(arg):
            return arg
        else:
            return self.func(arg)

Basic.singleton['asin'] = ASin
Basic.singleton['acos'] = ACos
Basic.singleton['atan'] = ATan
Basic.singleton['acot'] = ACot

###############################################################################
############################# HYPERBOLIC INVERSES #############################
###############################################################################

class ASinh(DefinedFunction):

    nofargs = 1

    def fdiff(self, argindex=1):
        if argindex == 1:
            z = Basic.Symbol('z', dummy=True)
            return Lambda((1 + z**2)**(-Basic.Half()), z)
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
                return Log()(Sqrt()(2) + 2)
            elif isinstance(arg, Basic.NegativeOne):
                return Log()(Sqrt()(2) - 2)
            elif arg.is_negative:
                return -self(-arg)
        else:
            i_coeff = arg.as_coefficient(S.ImaginaryUnit)

            if i_coeff is not None:
                return S.ImaginaryUnit * ASin()(i_coeff)
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
            return Basic.Zero()
        else:
            x = Basic.sympify(x)

            if len(previous_terms) > 2:
                p = previous_terms[-2]
                return -p * (n-2)**2/(k*(k-1)) * x**2
            else:
                k = (n - 1) // 2

                R = Basic.RisingFactorial()(Basic.Half(), k)
                F = Basic.Factorial()(k)

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
                    Sqrt()(2)/2  : S.Pi/4,
                    -Sqrt()(2)/2 : 3*S.Pi/4,
                    1/Sqrt()(2)  : S.Pi/4,
                    -1/Sqrt()(2) : 3*S.Pi/4,
                    Sqrt()(3)/2  : S.Pi/6,
                    -Sqrt()(3)/2 : 5*S.Pi/6,
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
            return Basic.Zero()
        else:
            x = Basic.sympify(x)

            if len(previous_terms) > 2:
                p = previous_terms[-2]
                return p * (n-2)**2/(k*(k-1)) * x**2
            else:
                k = (n - 1) // 2

                R = Basic.RisingFactorial()(Basic.Half(), k)
                F = Basic.Factorial()(k)

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
                return S.ImaginaryUnit * ATan()(i_coeff)
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
            return Basic.Zero()
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
                return -S.ImaginaryUnit * ACot()(i_coeff)
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
            return Basic.Zero()
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

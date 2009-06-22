
from basic import Basic, S, C
from sympify import _sympify
from cache import cacheit

from sympy import mpmath

from symbol import Symbol, Wild, Temporary
# from numbers import Number, Rational, Integer     /cyclic/
# from add import Add   /cyclic/
# from mul import Mul   /cyclic/

from math import exp as _exp
from math import log as _log

def integer_nthroot(y, n):
    """
    Return a tuple containing x = floor(y**(1/n))
    and a boolean indicating whether the result is exact (that is,
    whether x**n == y).

    >>> integer_nthroot(16,2)
    (4, True)
    >>> integer_nthroot(26,2)
    (5, False)
    """
    if y < 0: raise ValueError("y must be nonnegative")
    if n < 1: raise ValueError("n must be positive")
    if y in (0, 1): return y, True
    if n == 1: return y, True
    if n == 2:
        x, rem = mpmath.libmpf.sqrtrem(y)
        return int(x), not rem
    if n > y: return 1, False
    # Get initial estimate for Newton's method. Care must be taken to
    # avoid overflow
    try:
        guess = int(y ** (1./n)+0.5)
    except OverflowError:
        expt = _log(y,2)/n
        if expt > 53:
            shift = int(expt-53)
            guess = int(2.0**(expt-shift)+1) << shift
        else:
            guess = int(2.0**expt)
    #print n
    if guess > 2**50:
        # Newton iteration
        xprev, x = -1, guess
        while 1:
            t = x**(n-1)
            #xprev, x = x, x - (t*x-y)//(n*t)
            xprev, x = x, ((n-1)*x + y//t)//n
            #print n, x-xprev, abs(x-xprev) < 2
            if abs(x - xprev) < 2:
                break
    else:
        x = guess
    # Compensate
    t = x**n
    while t < y:
        x += 1
        t = x**n
    while t > y:
        x -= 1
        t = x**n
    return x, t == y

class Pow(Basic):

    is_Pow = True

    __slots__ = ['is_commutative']

    @cacheit
    def __new__(cls, a, b, **assumptions):
        a = _sympify(a)
        b = _sympify(b)
        if assumptions.get('evaluate') is False:
            return Basic.__new__(cls, a, b, **assumptions)
        if b is S.Zero:
            return S.One
        if b is S.One:
            return a
        obj = a._eval_power(b)
        if obj is None:
            obj = Basic.__new__(cls, a, b, **assumptions)
            obj.is_commutative = (a.is_commutative and b.is_commutative)
        return obj

    @property
    def base(self):
        return self._args[0]

    @property
    def exp(self):
        return self._args[1]

    def _eval_power(self, other):
        if other == S.NegativeOne:
            return Pow(self.base, self.exp * other)
        if self.exp.is_integer and other.is_integer:
            return Pow(self.base, self.exp * other)
        if self.base.is_nonnegative and self.exp.is_real and other.is_real:
            return Pow(self.base, self.exp * other)
        if self.exp.is_even and self.base.is_real:
            return Pow(abs(self.base), self.exp * other)
        if self.exp.is_real and other.is_real and abs(self.exp) < S.One:
            return Pow(self.base, self.exp * other)
        return

    def _eval_is_comparable(self):
        c1 = self.base.is_comparable
        if c1 is None: return
        c2 = self.exp.is_comparable
        if c2 is None: return
        return c1 and c2

    def _eval_is_even(self):
        if self.exp.is_integer and self.exp.is_positive:
            if self.base.is_even:
                return True
            if self.base.is_integer:
                return False

    def _eval_is_positive(self):
        if self.base.is_positive:
            if self.exp.is_real:
                return True
        elif self.base.is_negative:
            if self.exp.is_even:
                return True
            if self.exp.is_odd:
                return False
        elif self.base.is_nonpositive:
            if self.exp.is_odd:
                return False

    def _eval_is_negative(self):
        if self.base.is_negative:
            if self.exp.is_odd:
                return True
            if self.exp.is_even:
                return False
        elif self.base.is_positive:
            if self.exp.is_real:
                return False
        elif self.base.is_nonnegative:
            if self.exp.is_real:
                return False
        elif self.base.is_nonpositive:
            if self.exp.is_even:
                return False
        elif self.base.is_real:
            if self.exp.is_even:
                return False

    def _eval_is_integer(self):
        c1 = self.base.is_integer
        c2 = self.exp.is_integer
        if c1 is None or c2 is None:
            return None
        if not c1:
            if self.exp.is_nonnegative:
                return False
        if c1 and c2:
            if self.exp.is_nonnegative or self.exp.is_positive:
                return True
            if self.exp.is_negative:
                return False

    def _eval_is_real(self):
        c1 = self.base.is_real
        if c1 is None: return
        c2 = self.exp.is_real
        if c2 is None: return
        if c1 and c2:
            if self.base.is_positive:
                return True
            else:   # negative or zero (or positive)
                if self.exp.is_integer:
                    return True
                elif self.base.is_negative:
                    if self.exp.is_Rational:
                        return False

    def _eval_is_odd(self):
        if not (self.base.is_integer and self.exp.is_nonnegative): return
        return self.base.is_odd

    def _eval_is_bounded(self):
        if self.exp.is_negative:
            if self.base.is_infinitesimal:
                return False
            if self.base.is_unbounded:
                return True
        c1 = self.base.is_bounded
        if c1 is None: return
        c2 = self.exp.is_bounded
        if c2 is None: return
        if c1 and c2:
            if self.exp.is_nonnegative:
                return True

    def _eval_subs(self, old, new):
        if self==old: return new
        if isinstance(old, self.__class__) and self.base==old.base:
            coeff1,terms1 = self.exp.as_coeff_terms()
            coeff2,terms2 = old.exp.as_coeff_terms()
            if terms1==terms2: return new ** (coeff1/coeff2) # (x**(2*y)).subs(x**(3*y),z) -> z**(2/3*y)
        if old.func is C.exp:
            coeff1,terms1 = old.args[0].as_coeff_terms()
            coeff2,terms2 = (self.exp * C.log(self.base)).as_coeff_terms()
            if terms1==terms2: return new ** (coeff1/coeff2) # (x**(2*y)).subs(exp(3*y*log(x)),z) -> z**(2/3*y)
        return self.base._eval_subs(old, new) ** self.exp._eval_subs(old, new)

    def as_powers_dict(self):
        return { self.base : self.exp }

    def as_base_exp(self):
        if self.base.is_Rational and self.base.p==1:
            return 1/self.base, -self.exp
        return self.base, self.exp

    def _eval_conjugate(self):
        from sympy.functions.elementary.complexes import conjugate as c
        return c(self.base)**self.exp

    def _eval_expand_basic(self, deep=True, **hints):
        sargs, terms = self.args[:], []
        for term in sargs:
            if hasattr(term, '_eval_expand_basic'):
                newterm = term._eval_expand_basic(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

    def _eval_expand_power_exp(self, deep=True, *args, **hints):
        """a**(n+m) -> a**n*a**m"""
        if deep:
            b = self.base.expand(deep=deep, **hints)
            e = self.exp.expand(deep=deep, **hints)
        else:
            b = self.base
            e = self.exp
        if e.is_Add:
            expr = 1
            for x in e.args:
                if deep:
                    x = x.expand(deep=deep, **hints)
                expr *= (self.base**x)
            return expr
        return b**e

    def _eval_expand_power_base(self, deep=True, **hints):
        """(a*b)**n -> a**n * b**n"""
        b = self.base
        if deep:
            e = self.exp.expand(deep=deep, **hints)
        else:
            e = self.exp
        if b.is_Mul:
            if deep:
                return Mul(*(Pow(t.expand(deep=deep, **hints), e)\
                for t in b.args))
            else:
                return Mul(*(Pow(t, e) for t in b.args))
        else:
            return b**e

    def _eval_expand_mul(self, deep=True, **hints):
        sargs, terms = self.args[:], []
        for term in sargs:
            if hasattr(term, '_eval_expand_mul'):
                newterm = term._eval_expand_mul(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

    def _eval_expand_multinomial(self, deep=True, **hints):
        """(a+b+..) ** n -> a**n + n*a**(n-1)*b + .., n is positive integer"""
        if deep:
            b = self.base.expand(deep=deep, **hints)
            e = self.exp.expand(deep=deep, **hints)
        else:
            b = self.base
            e = self.exp

        if b is None:
            base = self.base
        else:
            base = b

        if e is None:
            exp = self.exp
        else:
            exp = e

        if e is not None or b is not None:
            result = base**exp

            if result.is_Pow:
                base, exp = result.base, result.exp
            else:
                return result
        else:
            result = None

        if exp.is_Integer and exp.p > 0 and base.is_Add:
            n = int(exp)

            if base.is_commutative:
                order_terms, other_terms = [], []

                for order in base.args:
                    if order.is_Order:
                        order_terms.append(order)
                    else:
                        other_terms.append(order)

                if order_terms:
                    # (f(x) + O(x^n))^m -> f(x)^m + m*f(x)^{m-1} *O(x^n)
                    f = Add(*other_terms)
                    g = (f**(n-1)).expand()

                    return (f*g).expand() + n*g*Add(*order_terms)

                if base.is_number:
                    # Efficiently expand expressions of the form (a + b*I)**n
                    # where 'a' and 'b' are real numbers and 'n' is integer.
                    a, b = base.as_real_imag()

                    if a.is_Rational and b.is_Rational:
                        if not a.is_Integer:
                            if not b.is_Integer:
                                k = (a.q * b.q) ** n
                                a, b = a.p*b.q, a.q*b.p
                            else:
                                k = a.q ** n
                                a, b = a.p, a.q*b
                        elif not b.is_Integer:
                            k = b.q ** n
                            a, b = a*b.q, b.p
                        else:
                            k = 1

                        a, b, c, d = int(a), int(b), 1, 0

                        while n:
                            if n & 1:
                                c, d = a*c-b*d, b*c+a*d
                                n -= 1
                            a, b = a*a-b*b, 2*a*b
                            n //= 2

                        I = S.ImaginaryUnit

                        if k == 1:
                            return c + I*d
                        else:
                            return Integer(c)/k + I*d/k

                p = other_terms
                # (x+y)**3 -> x**3 + 3*x**2*y + 3*x*y**2 + y**3
                # in this particular example:
                # p = [x,y]; n = 3
                # so now it's easy to get the correct result -- we get the
                # coefficients first:
                from sympy import multinomial_coefficients
                expansion_dict = multinomial_coefficients(len(p), n)
                # in our example: {(3, 0): 1, (1, 2): 3, (0, 3): 1, (2, 1): 3}
                # and now construct the expression.

                # An elegant way would be to use Poly, but unfortunately it is
                # slower than the direct method below, so it is commented out:
                #b = {}
                #for k in expansion_dict:
                #    b[k] = Integer(expansion_dict[k])
                #return Poly(b, *p).as_basic()

                from sympy.polys.polynomial import multinomial_as_basic
                result = multinomial_as_basic(expansion_dict, *p)
                return result
            else:
                if n == 2:
                    return Add(*[f*g for f in base.args for g in base.args])
                else:
                    return Mul(base, Pow(base, n-1).expand()).expand()
        elif exp.is_Add and base.is_Number:
            #  a + b      a  b
            # n      --> n  n  , where n, a, b are Numbers

            coeff, tail = S.One, S.Zero

            for term in exp.args:
                if term.is_Number:
                    coeff *= base**term
                else:
                    tail += term

            return coeff * base**tail
        else:
            return result

    def _eval_expand_log(self, deep=True, **hints):
        sargs, terms = self.args[:], []
        for term in sargs:
            if hasattr(term, '_eval_expand_log'):
                newterm = term._eval_expand_log(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

    def _eval_expand_complex(self, deep=True, **hints):
        if self.exp.is_Integer:
            exp = self.exp
            re, im = self.base.as_real_imag()
            if exp >= 0:
                base = re + S.ImaginaryUnit*im
            else:
                mag = re**2 + im**2
                base = re/mag - S.ImaginaryUnit*(im/mag)
                exp = -exp
            return (base**exp).expand()
        elif self.exp.is_Rational:
            # NOTE: This is not totally correct since for x**(p/q) with
            #       x being imaginary there are actually q roots, but
            #       only a single one is returned from here.
            re, im = self.base.as_real_imag()

            r = (re**2 + im**2)**S.Half
            t = C.atan2(im, re)

            rp, tp = r**self.exp, t*self.exp

            return rp*C.cos(tp) + rp*C.sin(tp)*S.ImaginaryUnit
        else:
            if deep:
                hints['complex'] = False
                return C.re(self.expand(deep, **hints)) + \
                S.ImaginaryUnit*C.im(self. expand(deep, **hints))
            else:
                return C.re(self) + S.ImaginaryUnit*C.im(self)
            return C.re(self) + S.ImaginaryUnit*C.im(self)

    def _eval_expand_trig(self, deep=True, **hints):
        sargs, terms = self.args[:], []
        for term in sargs:
            if hasattr(term, '_eval_expand_trig'):
                newterm = term._eval_expand_trig(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

    def _eval_expand_func(self, deep=True, **hints):
        sargs, terms = self.args[:], []
        for term in sargs:
            if hasattr(term, '_eval_expand_func'):
                newterm = term._eval_expand_func(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

    def _eval_derivative(self, s):
        dbase = self.base.diff(s)
        dexp = self.exp.diff(s)
        return self * (dexp * C.log(self.base) + dbase * self.exp/self.base)

    def _eval_evalf(self, prec):
        base, exp = self.as_base_exp()
        base = base._evalf(prec)
        if not exp.is_Integer:
            exp = exp._evalf(prec)
        if exp < 0 and not base.is_real:
            base = base.conjugate() / (base * base.conjugate())._evalf(prec)
            exp = -exp
        return (base ** exp).expand()

    @cacheit
    def count_ops(self, symbolic=True):
        if symbolic:
            return Add(*[t.count_ops(symbolic) for t in self.args]) + Symbol('POW')
        return Add(*[t.count_ops(symbolic) for t in self.args]) + 1

    def _eval_is_polynomial(self, syms):
        if self.exp.has(*syms):
            return False

        if self.base.has(*syms):
            # it would be nice to have is_nni working
            return self.base._eval_is_polynomial(syms) and \
                   self.exp.is_nonnegative and \
                   self.exp.is_integer
        else:
            return True

    def as_numer_denom(self):
        base, exp = self.as_base_exp()
        c,t = exp.as_coeff_terms()
        n,d = base.as_numer_denom()
        negate = False
        if exp.is_integer != True:
            if d.is_negative == True:
                # Roots need to take care that negative denominators behave
                # differently than the rest of the complex plane.
                negate = True
            elif d.is_negative is None:
                # Can make no conclusions.
                return self, S(1)
        if c.is_negative == True:
            exp = -exp
            n,d = d,n
        num = n ** exp
        den = d ** exp
        if negate:
            num = -num
        return num, den

    def matches(pattern, expr, repl_dict={}, evaluate=False):
        if evaluate:
            pat = pattern
            for old,new in repl_dict.items():
                pat = pat.subs(old, new)
            if pat!=pattern:
                return pat.matches(expr, repl_dict)

        expr = _sympify(expr)
        b, e = expr.as_base_exp()

        # special case, pattern = 1 and expr.exp can match to 0
        if expr is S.One:
            d = repl_dict.copy()
            d = pattern.exp.matches(S.Zero, d, evaluate=False)
            if d is not None:
                return d

        d = repl_dict.copy()
        d = pattern.base.matches(b, d, evaluate=False)
        if d is None:
            return None

        d = pattern.exp.matches(e, d, evaluate=True)
        if d is None:
            return Basic.matches(pattern, expr, repl_dict, evaluate)
        return d

    def _eval_nseries(self, x, x0, n):
        from sympy import powsimp, collect

        def geto(e):
            "Returns the O(..) symbol, or None if there is none."
            if e.is_Order:
                return e
            if e.is_Add:
                for x in e.args:
                    if x.is_Order:
                        return x

        def getn(e):
            """
            Returns the order of the expression "e".

            The order is determined either from the O(...) term. If there
            is no O(...) term, it returns None.

            Example:
            >>> getn(1+x+O(x**2))
            2
            >>> getn(1+x)
            >>>
            """
            o = geto(e)
            if o is None:
                return None
            else:
                o = o.expr
                if o.is_Symbol:
                    return Integer(1)
                if o.is_Pow:
                    return o.args[1]
                n, d = o.as_numer_denom()
                if isinstance(d, log):
                    # i.e. o = x**2/log(x)
                    if n.is_Symbol:
                        return Integer(1)
                    if n.is_Pow:
                        return n.args[1]

            raise NotImplementedError()

        base, exp = self.args
        if exp.is_Integer:
            if exp > 0:
                # positive integer powers are easy to expand, e.g.:
                # sin(x)**4 = (x-x**3/3+...)**4 = ...
                return (base.nseries(x, x0, n) ** exp).expand()
            elif exp == -1:
                # this is also easy to expand using the formula:
                # 1/(1 + x) = 1 + x + x**2 + x**3 ...
                # so we need to rewrite base to the form "1+x"
                from sympy import log
                if base.has(log(x)):
                    # we need to handle the log(x) singularity:
                    assert x0 == 0
                    y = Symbol("y", dummy=True)
                    p = self.subs(log(x), -1/y)
                    if not p.has(x):
                        p = p.nseries(y, x0, n)
                        p = p.subs(y, -1/log(x))
                        return p

                base = base.nseries(x, x0, n)
                if base.has(log(x)):
                    # we need to handle the log(x) singularity:
                    assert x0 == 0
                    y = Symbol("y", dummy=True)
                    self0 = 1/base
                    p = self0.subs(log(x), -1/y)
                    if not p.has(x):
                        p = p.nseries(y, x0, n)
                        p = p.subs(y, -1/log(x))
                        return p
                prefactor = base.as_leading_term(x)
                # express "rest" as: rest = 1 + k*x**l + ... + O(x**n)
                rest = powsimp(((base-prefactor)/prefactor).expand(),\
                deep=True, combine='exp')
                if rest == 0:
                    # if prefactor == w**4 + x**2*w**4 + 2*x*w**4, we need to
                    # factor the w**4 out using collect:
                    return 1/collect(prefactor, x)
                if rest.is_Order:
                    return ((1+rest)/prefactor).expand()
                n2 = getn(rest)
                if n2 is not None:
                    n = n2

                term2 = collect(rest.as_leading_term(x), x)
                k, l = Wild("k"), Wild("l")
                r = term2.match(k*x**l)
                k, l = r[k], r[l]
                if l.is_Integer and l>0:
                    l = int(l)
                elif l.is_number and l>0:
                    l = float(l)
                else:
                    raise NotImplementedError()

                s = 1
                m = 1
                while l * m < n:
                    s += ((-rest)**m).expand()
                    m += 1
                r = (s/prefactor).expand()
                if n2 is None:
                    # Append O(...) because it is not included in "r"
                    from sympy import O
                    r += O(x**n)
                return powsimp(r, deep=True, combine='exp')
            else:
                # negative powers are rewritten to the cases above, for example:
                # sin(x)**(-4) = 1/( sin(x)**4) = ...
                # and expand the denominator:
                denominator = (base**(-exp)).nseries(x, x0, n)
                if 1/denominator == self:
                    return self
                # now we have a type 1/f(x), that we know how to expand
                return (1/denominator).nseries(x, x0, n)

        if exp.has(x):
            import sympy
            return sympy.exp(exp*sympy.log(base)).nseries(x, x0, n)

        if base == x:
            return powsimp(self, deep=True, combine='exp')

        order = C.Order(x**n, x)
        x = order.symbols[0]
        e = self.exp
        b = self.base
        ln = C.log
        exp = C.exp
        if e.has(x):
            return exp(e * ln(b)).nseries(x, x0, n)
        if b==x:
            return self
        b0 = b.limit(x,0)
        if b0 is S.Zero or b0.is_unbounded:
            lt = b.as_leading_term(x)
            o = order * lt**(1-e)
            bs = b.nseries(x, x0, n-e)
            if bs.is_Add:
                bs = bs.removeO()
            if bs.is_Add:
                # bs -> lt + rest -> lt * (1 + (bs/lt - 1))
                return (lt**e * ((bs/lt).expand()**e).nseries(x,
                        x0, n-e)).expand() + order

            return bs**e+order
        o2 = order * (b0**-e)
        # b -> b0 + (b-b0) -> b0 * (1 + (b/b0-1))
        z = (b/b0-1)
        #r = self._compute_oseries3(z, o2, self.taylor_term)
        x = o2.symbols[0]
        ln = C.log
        o = C.Order(z, x)
        if o is S.Zero:
            r = (1+z)
        else:
            if o.expr.is_number:
                e2 = ln(o2.expr*x)/ln(x)
            else:
                e2 = ln(o2.expr)/ln(o.expr)
            n = e2.limit(x,0) + 1
            if n.is_unbounded:
                # requested accuracy gives infinite series,
                # order is probably nonpolynomial e.g. O(exp(-1/x), x).
                r = (1+z)
            else:
                try:
                    n = int(n)
                except TypeError:
                    #well, the n is something more complicated (like 1+log(2))
                    n = int(n.evalf()) + 1
                assert n>=0,`n`
                l = []
                g = None
                for i in xrange(n+2):
                    g = self.taylor_term(i, z, g)
                    g = g.nseries(x, x0, n)
                    l.append(g)
                r = Add(*l)
        return r * b0**e + order

    def _eval_as_leading_term(self, x):
        if not self.exp.has(x):
            return self.base.as_leading_term(x) ** self.exp
        return C.exp(self.exp * C.log(self.base)).as_leading_term(x)

    @cacheit
    def taylor_term(self, n, x, *previous_terms): # of (1+x)**e
        if n<0: return S.Zero
        x = _sympify(x)
        return C.Binomial(self.exp, n) * x**n

    def _sage_(self):
        return self.args[0]._sage_() ** self.args[1]._sage_()


# /cyclic/
import basic as _
_.Pow =     Pow
del _

import mul as _
_.Pow =     Pow
del _

import numbers as _
_.Pow =     Pow
del _


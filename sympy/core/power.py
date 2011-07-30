from math import log as _log

from sympify import _sympify
from cache import cacheit
from core import C
from singleton import S
from expr import Expr

from sympy import mpmath
from sympy.utilities.iterables import sift

def integer_nthroot(y, n):
    """
    Return a tuple containing x = floor(y**(1/n))
    and a boolean indicating whether the result is exact (that is,
    whether x**n == y).

    >>> from sympy import integer_nthroot
    >>> integer_nthroot(16,2)
    (4, True)
    >>> integer_nthroot(26,2)
    (5, False)

    """
    y, n = int(y), int(n)
    if y < 0: raise ValueError("y must be nonnegative")
    if n < 1: raise ValueError("n must be positive")
    if y in (0, 1): return y, True
    if n == 1: return y, True
    if n == 2:
        x, rem = mpmath.libmp.sqrtrem(y)
        return int(x), not rem
    if n > y: return 1, False
    # Get initial estimate for Newton's method. Care must be taken to
    # avoid overflow
    try:
        guess = int(y**(1./n) + 0.5)
    except OverflowError:
        expt = _log(y, 2)/n
        if expt > 53:
            shift = int(expt - 53)
            guess = int(2.0**(expt-shift) + 1) << shift
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

class Pow(Expr):

    is_Pow = True

    __slots__ = ['is_commutative']

    @cacheit
    def __new__(cls, b, e, evaluate=True):
        b = _sympify(b)
        e = _sympify(e)
        if evaluate:
            if e is S.Zero:
                return S.One
            elif e is S.One:
                return b
            else:
                obj = b._eval_power(e)
                if obj is not None:
                    return obj
        obj = Expr.__new__(cls, b, e)
        obj.is_commutative = (b.is_commutative and e.is_commutative)
        return obj

    @property
    def base(self):
        return self._args[0]

    @property
    def exp(self):
        return self._args[1]

    @classmethod
    def class_key(cls):
        return 3, 2, cls.__name__

    def _eval_power(self, other):
        b, e = self.as_base_exp()
        if other.is_integer:
            return Pow(b, e * other)
        if b.is_nonnegative and (e.is_real or other.is_real):
            return Pow(b, e * other)
        if e.is_even and b.is_real: # hence b is pos and e is real
            return Pow(abs(b), e * other)
        if abs(e) < S.One and other.is_real:
            return Pow(b, e * other)

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
        real_b = self.base.is_real
        if real_b is None: return
        real_e = self.exp.is_real
        if real_e is None: return
        if real_b and real_e:
            if self.base.is_positive:
                return True
            else:   # negative or zero (or positive)
                if self.exp.is_integer:
                    return True
                elif self.base.is_negative:
                    if self.exp.is_Rational:
                        return False
        im_b = self.base.is_imaginary
        im_e = self.exp.is_imaginary
        if im_b:
            if self.exp.is_integer:
                if self.exp.is_even:
                    return True
                elif self.exp.is_odd:
                    return False
            elif (self.exp in [S.ImaginaryUnit, -S.ImaginaryUnit] and
                  self.base in [S.ImaginaryUnit, -S.ImaginaryUnit]):
                return True
        if real_b and im_e:
            c = self.exp.coeff(S.ImaginaryUnit)
            if c:
                ok = (c*C.log(self.base)/S.Pi).is_Integer
                if ok is not None:
                    return ok

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
            if self.exp.is_nonnegative or self.base.is_nonzero:
                return True

    def _eval_subs(self, old, new):
        if self == old:
            return new
        if old.func is self.func and self.base == old.base:
            coeff1, terms1 = self.exp.as_coeff_mul()
            coeff2, terms2 = old.exp.as_coeff_mul()
            if terms1 == terms2:
                pow = coeff1/coeff2
                if pow.is_Integer or self.base.is_commutative:
                    return Pow(new, pow) # (x**(2*y)).subs(x**(3*y),z) -> z**(2/3)
        if old.func is C.exp:
            coeff1, terms1 = old.args[0].as_coeff_mul()
            coeff2, terms2 = (self.exp*C.log(self.base)).as_coeff_mul()
            if terms1 == terms2:
                pow = coeff1/coeff2
                if pow.is_Integer or self.base.is_commutative:
                    return Pow(new, pow) # (x**(2*y)).subs(x**(3*y),z) -> z**(2/3)
        return Pow(self.base._eval_subs(old, new), self.exp._eval_subs(old, new))

    def as_base_exp(self):
        if self.base.is_Rational and self.base.p==1 and self.base is not S.Infinity:
            return 1/self.base, -self.exp
        return self.base, self.exp

    def _eval_conjugate(self):
        from sympy.functions.elementary.complexes import conjugate as c
        return c(self.base)**self.exp

    def _eval_expand_basic(self, deep=True, **hints):
        sargs, terms = self.args, []
        for term in sargs:
            if hasattr(term, '_eval_expand_basic'):
                newterm = term._eval_expand_basic(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.func(*terms)

    def _eval_expand_power_exp(self, deep=True, *args, **hints):
        """a**(n+m) -> a**n*a**m"""
        if deep:
            b = self.base.expand(deep=deep, **hints)
            e = self.exp.expand(deep=deep, **hints)
        else:
            b = self.base
            e = self.exp
        if e.is_Add and e.is_commutative:
            expr = []
            for x in e.args:
                if deep:
                    x = x.expand(deep=deep, **hints)
                expr.append(Pow(self.base, x))
            return Mul(*expr)
        return Pow(b, e)

    def _eval_expand_power_base(self, deep=True, **hints):
        """(a*b)**n -> a**n * b**n"""
        force = hints.get('force', False)
        b, ewas = self.args
        if deep:
            e = self.exp.expand(deep=deep, **hints)
        else:
            e = self.exp
        if b.is_Mul:
            bargs = b.args
            if force or e.is_integer:
                nonneg = bargs
                other = []
            elif ewas.is_Rational or len(bargs) == 2 and bargs[0] is S.NegativeOne:
                # the Rational exponent was already expanded automatically
                # if there is a negative Number * foo, foo must be unknown
                #    or else it, too, would have automatically expanded;
                #    sqrt(-Number*foo) -> sqrt(Number)*sqrt(-foo); then
                #    sqrt(-foo) -> unchanged if foo is not positive else
                #               -> I*sqrt(foo)
                #    So...if we have a 2 arg Mul and the first is a Number
                #    that number is -1 and there is nothing more than can
                #    be done without the force=True hint
                nonneg= []
            else:
                # this is just like what is happening automatically, except
                # that now we are doing it for an arbitrary exponent for which
                # no automatic expansion is done
                sifted = sift(b.args, lambda x: x.is_nonnegative)
                nonneg = sifted.get(True, [])
                other = sifted.get(None, [])
                neg = sifted.get(False, [])

                # make sure the Number gets pulled out
                if neg and neg[0].is_Number and neg[0] is not S.NegativeOne:
                    nonneg.append(-neg[0])
                    neg[0] = S.NegativeOne

                # leave behind a negative sign
                oddneg = len(neg) % 2
                if oddneg:
                    other.append(S.NegativeOne)

                # negate all negatives and append to nonneg
                nonneg += [-n for n in neg]

            if nonneg: # then there's a new expression to return
                other = [Pow(Mul(*other), e)]
                if deep:
                    return Mul(*([Pow(b.expand(deep=deep, **hints), e)\
                    for b in nonneg] + other))
                else:
                    return Mul(*([Pow(b, e) for b in nonneg] + other))
        return Pow(b, e)

    def _eval_expand_mul(self, deep=True, **hints):
        sargs, terms = self.args, []
        for term in sargs:
            if hasattr(term, '_eval_expand_mul'):
                newterm = term._eval_expand_mul(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.func(*terms)

    def _eval_expand_multinomial(self, deep=True, **hints):
        """(a+b+..) ** n -> a**n + n*a**(n-1)*b + .., n is nonzero integer"""
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
            result = Pow(base, exp)

            if result.is_Pow:
                base, exp = result.base, result.exp
            else:
                return result
        else:
            result = None

        if exp.is_Rational and exp.p > 0 and base.is_Add:
            if not exp.is_Integer:
                n = Integer(exp.p // exp.q)

                if not n:
                    return Pow(base, exp)
                else:
                    radical, result = Pow(base, exp - n), []

                    for term in Add.make_args(Pow(base, n)._eval_expand_multinomial(deep=False)):
                        result.append(term*radical)

                    return Add(*result)

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
                                k = Pow(a.q * b.q, n)
                                a, b = a.p*b.q, a.q*b.p
                            else:
                                k = Pow(a.q, n)
                                a, b = a.p, a.q*b
                        elif not b.is_Integer:
                            k = Pow(b.q, n)
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
                #return Poly(b, *p).as_expr()

                from sympy.polys.polyutils import basic_from_dict
                result = basic_from_dict(expansion_dict, *p)
                return result
            else:
                if n == 2:
                    return Add(*[f*g for f in base.args for g in base.args])
                else:
                    multi = (base**(n-1))._eval_expand_multinomial(deep=False)
                    if multi.is_Add:
                        return Add(*[f*g for f in base.args for g in multi.args])
                    else:
                        return Add(*[f*multi for f in base.args])
        elif exp.is_Rational and exp.p < 0 and base.is_Add and abs(exp.p) > exp.q:
            return 1 / Pow(base, -exp)._eval_expand_multinomial(deep=False)
        elif exp.is_Add and base.is_Number:
            #  a + b      a  b
            # n      --> n  n  , where n, a, b are Numbers

            coeff, tail = S.One, S.Zero

            for term in exp.args:
                if term.is_Number:
                    coeff *= Pow(base, term)
                else:
                    tail += term

            return coeff * Pow(base, tail)
        else:
            return result

    def _eval_expand_log(self, deep=True, **hints):
        sargs, terms = self.args, []
        for term in sargs:
            if hasattr(term, '_eval_expand_log'):
                newterm = term._eval_expand_log(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.func(*terms)

    def _eval_expand_complex(self, deep=True, **hints):
        re_part, im_part = self.as_real_imag(deep=deep, **hints)
        return re_part + im_part*S.ImaginaryUnit

    def as_real_imag(self, deep=True, **hints):
        from sympy.core.symbol import symbols
        from sympy.polys.polytools import poly
        from sympy.core.function import expand_multinomial
        if self.exp.is_Integer:
            exp = self.exp
            re, im = self.base.as_real_imag(deep=deep)
            a, b = symbols('a b', cls=Dummy)
            if exp >= 0:
                if re.is_Number and im.is_Number:
                    # We can be more efficient in this case
                    expr = expand_multinomial(self.base**exp)
                    return expr.as_real_imag()

                expr = poly((a + b)**exp) # a = re, b = im; expr = (a + b*I)**exp
            else:
                mag = re**2 + im**2
                re, im = re/mag, -im/mag
                if re.is_Number and im.is_Number:
                    # We can be more efficient in this case
                    expr = expand_multinomial((re + im*S.ImaginaryUnit)**-exp)
                    return expr.as_real_imag()

                expr = poly((a + b)**-exp)

            # Terms with even b powers will be real
            r = [i for i in expr.terms() if not i[0][1] % 2]
            re_part = Add(*[cc*a**aa*b**bb for (aa, bb), cc in r])
            # Terms with odd b powers will be imaginary
            r = [i for i in expr.terms() if i[0][1] % 4 == 1]
            im_part1 = Add(*[cc*a**aa*b**bb for (aa, bb), cc in r])
            r = [i for i in expr.terms() if i[0][1] % 4 == 3]
            im_part3 = Add(*[cc*a**aa*b**bb for (aa, bb), cc in r])

            return (re_part.subs({a: re, b: S.ImaginaryUnit*im}),
            im_part1.subs({a: re, b: im}) + im_part3.subs({a: re, b: -im}))

        elif self.exp.is_Rational:
            # NOTE: This is not totally correct since for x**(p/q) with
            #       x being imaginary there are actually q roots, but
            #       only a single one is returned from here.
            re, im = self.base.as_real_imag(deep=deep)

            r = Pow(Pow(re, 2) + Pow(im, 2), S.Half)
            t = C.atan2(im, re)

            rp, tp = Pow(r, self.exp), t*self.exp

            return (rp*C.cos(tp), rp*C.sin(tp))
        else:

            if deep:
                hints['complex'] = False
                return (C.re(self.expand(deep, complex=False)),
                C.im(self. expand(deep, **hints)))
            else:
                return (C.re(self), C.im(self))

    def _eval_expand_trig(self, deep=True, **hints):
        sargs, terms = self.args, []
        for term in sargs:
            if hasattr(term, '_eval_expand_trig'):
                newterm = term._eval_expand_trig(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.func(*terms)

    def _eval_expand_func(self, deep=True, **hints):
        sargs, terms = self.args, []
        for term in sargs:
            if hasattr(term, '_eval_expand_func'):
                newterm = term._eval_expand_func(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.func(*terms)

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
        return Pow(base, exp).expand()

    def _eval_is_polynomial(self, syms):
        if self.exp.has(*syms):
            return False

        if self.base.has(*syms):
            return self.base._eval_is_polynomial(syms) and \
                   self.exp.is_Integer and \
                   self.exp >= 0
        else:
            return True

    def _eval_is_rational_function(self, syms):
        if self.exp.has(*syms):
            return False

        if self.base.has(*syms):
            return self.base._eval_is_rational_function(syms) and \
                   self.exp.is_Integer
        else:
            return True

    def as_numer_denom(self):
        base, exp = self.as_base_exp()
        n, d = base.as_numer_denom()
        if d.is_negative and n.is_negative:
            n, d = -n, -d
        if exp.is_Integer:
            if exp.is_negative:
                n, d = d, n
                exp = -exp
            return Pow(n, exp), Pow(d, exp)
        elif exp.is_Rational or d.is_positive:
            if d.is_negative is None:
                # we won't split up the base
                if exp.is_negative:
                    return S.One, Pow(base, -exp)
                else:
                    return self, S.One
            if d.is_negative:
                n = -n
                d = -d
            c, t = exp.as_coeff_mul()
            if c.is_negative:
                n, d = d, n
                exp = -exp
            return Pow(n, exp), Pow(d, exp)
        else:
            c, t = exp.as_coeff_mul()
            if c.is_negative:
                return 1, base**-exp
        # unprocessed Float and NumberSymbol
        return self, S.One

    def matches(self, expr, repl_dict={}, evaluate=False):
        if evaluate:
            return self.subs(repl_dict).matches(expr, repl_dict)

        expr = _sympify(expr)
        b, e = expr.as_base_exp()

        # special case, pattern = 1 and expr.exp can match to 0
        if expr is S.One:
            d = repl_dict.copy()
            d = self.exp.matches(S.Zero, d)
            if d is not None:
                return d

        d = repl_dict.copy()
        d = self.base.matches(b, d)
        if d is None:
            return None

        d = self.exp.subs(d).matches(e, d)
        if d is None:
            return Expr.matches(self, expr, repl_dict, evaluate)
        return d

    def _eval_nseries(self, x, n, logx):
        # NOTE! This function is an important part of the gruntz algorithm
        #       for computing limits. It has to return a generalised power series
        #       with coefficients in C(log, log(x)). In more detail:
        # It has to return an expression
        #     c_0*x**e_0 + c_1*x**e_1 + ... (finitely many terms)
        # where e_i are numbers (not necessarily integers) and c_i are expression
        # involving only numbers, the log function, and log(x).
        from sympy import powsimp, collect, exp, log, O, ceiling
        b, e = self.args
        if e.is_Integer:
            if e > 0:
                # positive integer powers are easy to expand, e.g.:
                # sin(x)**4 = (x-x**3/3+...)**4 = ...
                return Pow(b._eval_nseries(x, n=n, logx=logx), e
                           )._eval_expand_multinomial(deep = False)
            elif e is S.NegativeOne:
                # this is also easy to expand using the formula:
                # 1/(1 + x) = 1 - x + x**2 - x**3 ...
                # so we need to rewrite base to the form "1+x"

                b = b._eval_nseries(x, n=n, logx=logx)
                prefactor = b.as_leading_term(x)
                # express "rest" as: rest = 1 + k*x**l + ... + O(x**n)
                rest = ((b - prefactor)/prefactor)._eval_expand_mul()
                if rest == 0:
                    # if prefactor == w**4 + x**2*w**4 + 2*x*w**4, we need to
                    # factor the w**4 out using collect:
                    return 1/collect(prefactor, x)
                if rest.is_Order:
                    return 1/prefactor + rest/prefactor
                n2 = rest.getn()
                if n2 is not None:
                    n = n2
                    # remove the O - powering this is slow
                    if logx is not None:
                        rest = rest.removeO()

                k, l = rest.leadterm(x)
                if l.is_Rational and l > 0:
                    pass
                elif l.is_number and l > 0:
                    l = l.evalf()
                else:
                    raise NotImplementedError()

                terms = [1/prefactor]
                for m in xrange(1, ceiling(n/l)):
                    new_term = terms[-1]*(-rest)
                    if new_term.is_Pow:
                        new_term = new_term._eval_expand_multinomial(deep = False)
                    else:
                        new_term = new_term._eval_expand_mul(deep = False)
                    terms.append(new_term)

                # Append O(...), we know the order.
                if n2 is None or logx is not None:
                    terms.append(O(x**n))
                return powsimp(Add(*terms), deep=True, combine='exp')
            else:
                # negative powers are rewritten to the cases above, for example:
                # sin(x)**(-4) = 1/( sin(x)**4) = ...
                # and expand the denominator:
                denominator = (b**(-e))._eval_nseries(x, n=n, logx=logx)
                if 1/denominator == self:
                    return self
                # now we have a type 1/f(x), that we know how to expand
                return (1/denominator)._eval_nseries(x, n=n, logx=logx)

        if e.has(Symbol):
            return exp(e*log(b))._eval_nseries(x, n=n, logx=logx)

        if b == x:
            return powsimp(self, deep=True, combine='exp')

        # work for b(x)**e where e is not an Integer and does not contain x
        # and hopefully has no other symbols

        def e2int(e):
            """return the integer value (if possible) of e and a
            flag indicating whether it is bounded or not."""
            n = e.limit(x, 0)
            unbounded = n.is_unbounded
            if not unbounded:
                # XXX was int or floor intended? int used to behave like floor
                # so int(-Rational(1, 2)) returned -1 rather than int's 0
                try:
                    n = int(n)
                except TypeError:
                    #well, the n is something more complicated (like 1+log(2))
                    try:
                        n = int(n.evalf()) + 1 # XXX why is 1 being added?
                    except TypeError:
                        pass # hope that base allows this to be resolved
                n = _sympify(n)
            return n, unbounded

        order = O(x**n, x)
        ei, unbounded = e2int(e)
        b0 = b.limit(x, 0)
        if unbounded and (b0 is S.One or b0.has(Symbol)):
            # XXX what order
            if b0 is S.One:
                resid = (b - 1)
                if resid.is_positive:
                    return S.Infinity
                elif resid.is_negative:
                    return S.Zero
                raise ValueError('cannot determine sign of %s' % resid)

            return b0**ei


        if (b0 is S.Zero or b0.is_unbounded):
            if unbounded is not False:
                return b0**e # XXX what order

            if not ei.is_number: # if not, how will we proceed?
                raise ValueError('expecting numerical exponent but got %s' % ei)

            nuse = n - ei
            lt = b.compute_leading_term(x, logx=logx) # arg = sin(x); lt = x
            #  XXX o is not used -- was this to be used as o and o2 below to compute a new e?
            o = order*lt**(1 - e)
            bs = b._eval_nseries(x, n=nuse, logx=logx)
            if bs.is_Add:
                bs = bs.removeO()
            if bs.is_Add:
                # bs -> lt + rest -> lt*(1 + (bs/lt - 1))
                return ((Pow(lt, e)*
                         Pow((bs/lt).expand(), e).
                         nseries(x, n=nuse, logx=logx)).expand() +
                         order)

            return bs**e + order

        # either b0 is bounded but neither 1 nor 0 or e is unbounded
        # b -> b0 + (b-b0) -> b0 * (1 + (b/b0-1))
        o2 = order*(b0**-e)
        z = (b/b0 - 1)
        o = O(z, x)
        #r = self._compute_oseries3(z, o2, self.taylor_term)
        if o is S.Zero or o2 is S.Zero:
            unbounded = True
        else:
            if o.expr.is_number:
                e2 = log(o2.expr*x)/log(x)
            else:
                e2 = log(o2.expr)/log(o.expr)
            n, unbounded = e2int(e2)
        if unbounded:
            # requested accuracy gives infinite series,
            # order is probably nonpolynomial e.g. O(exp(-1/x), x).
            r = 1 + z
        else:
            l = []
            g = None
            for i in xrange(n + 2):
                g = self.taylor_term(i, z, g)
                g = g.nseries(x, n=n, logx=logx)
                l.append(g)
            r = Add(*l)
        return r*b0**e + order

    def _eval_as_leading_term(self, x):
        if not self.exp.has(x):
            return Pow(self.base.as_leading_term(x), self.exp)
        return C.exp(self.exp * C.log(self.base)).as_leading_term(x)

    @cacheit
    def taylor_term(self, n, x, *previous_terms): # of (1+x)**e
        if n<0: return S.Zero
        x = _sympify(x)
        return C.binomial(self.exp, n) * Pow(x, n)

    def _sage_(self):
        return self.args[0]._sage_()**self.args[1]._sage_()

from add import Add
from numbers import Integer
from mul import Mul
from symbol import Symbol, Dummy


from basic import Basic, S, cache_it, cache_it_immutable
from operations import AssocOp
from methods import RelMeths, ArithMeths

class Mul(AssocOp, RelMeths, ArithMeths):

    @classmethod
    def flatten(cls, seq):
        # apply associativity, separate commutative part of seq
        c_part = []
        nc_part = []
        coeff = Basic.One()
        c_seq = []
        nc_seq = seq
        c_powers = {}
        lambda_args = None
        order_symbols = None
        while c_seq or nc_seq:
            if c_seq:
                # first process commutative objects
                o = c_seq.pop(0)
                if isinstance(o, Basic.Function):
                    if o.nofargs is not None:
                        o, lambda_args = o.with_dummy_arguments(lambda_args)
                if isinstance(o, Basic.Order):
                    o, order_symbols = o.as_expr_symbols(order_symbols)
                if o.__class__ is cls:
                    # associativity
                    c_seq = list(o._args) + c_seq
                    continue
                if isinstance(o, Basic.Number):
                    coeff *= o
                    continue
                if isinstance(o, Basic.ApplyExp):
                    # exp(x) / exp(y) -> exp(x-y)
                    b = Basic.Exp1()
                    e = o.args[0]
                else:
                    b, e = o.as_base_exp()
                if isinstance(b, Basic.Add) and isinstance(e, Basic.Number):
                    c, t = b.as_coeff_terms()
                    if not isinstance(c, Basic.One):
                        coeff *= c ** e
                        assert len(t)==1,`t`
                        b = t[0]

                if b in c_powers:
                    c_powers[b] += e
                else:
                    c_powers[b] = e
            else:
                o = nc_seq.pop(0)
                if isinstance(o, Basic.WildFunction):
                    pass
                elif isinstance(o, Basic.Function):
                    if o.nofargs is not None:
                        o, lambda_args = o.with_dummy_arguments(lambda_args)
                elif isinstance(o, Basic.Order):
                    o, order_symbols = o.as_expr_symbols(order_symbols)

                if o.is_commutative:
                    # separate commutative symbols
                    c_seq.append(o)
                    continue
                if o.__class__ is cls:
                    # associativity
                    nc_seq = list(o._args) + nc_seq
                    continue
                if not nc_part:
                    nc_part.append(o)
                    continue
                # try to combine last terms: a**b * a ** c -> a ** (b+c)
                o1 = nc_part.pop()
                b1,e1 = o1.as_base_exp()
                b2,e2 = o.as_base_exp()
                if b1==b2:
                    nc_seq.insert(0, b1 ** (e1 + e2))
                else:
                    nc_part.append(o1)
                    nc_part.append(o)

        for b, e in c_powers.items():
            if isinstance(e, Basic.Zero):
                continue

            if isinstance(e, Basic.One):
                if isinstance(b, Basic.Number):
                    coeff *= b
                else:
                    c_part.append(b)
            elif isinstance(e, Basic.Integer) and isinstance(b, Basic.Number):
                coeff *= b ** e
            else:
                c_part.append(Basic.Pow(b, e))

        if isinstance(coeff, (Basic.Infinity, Basic.NegativeInfinity)):
            new_c_part = []
            for t in c_part:
                if t.is_positive:
                    continue
                if t.is_negative:
                    coeff = -coeff
                    continue
                new_c_part.append(t)
            c_part = new_c_part
            new_nc_part = []
            for t in nc_part:
                if t.is_positive:
                    continue
                if t.is_negative:
                    coeff = -coeff
                    continue
                new_nc_part.append(t)
            nc_part = new_nc_part
            c_part.insert(0, coeff)
        elif isinstance(coeff, (Basic.Zero, Basic.NaN)):
            c_part, nc_part = [coeff], []
        elif not isinstance(coeff, Basic.One):
            c_part.insert(0, coeff)

        c_part.sort(Basic.compare)
        if len(c_part)==2 and isinstance(c_part[0], Basic.Number) and isinstance(c_part[1], Basic.Add):
            # 2*(1+a) -> 2 + 2 * a
            coeff = c_part[0]
            c_part = [Basic.Add(*[coeff*f for f in c_part[1]])]
        return c_part, nc_part, lambda_args, order_symbols

    def _eval_power(b, e):
        if isinstance(e, Basic.Number):
            if b.is_commutative:
                if isinstance(e, Basic.Integer):
                    # (a*b)**2 -> a**2 * b**2
                    return Mul(*[s**e for s in b])

                if e.is_rational and not b.is_nonnegative:
                    if not isinstance(b, (Basic.Pow, Basic.Number)):
                        return

                coeff, rest = b.as_coeff_terms()
                if not isinstance(coeff, Basic.One):
                    # (2*a)**3 -> 2**3 * a**3
                    return coeff**e * Mul(*[s**e for s in rest])
            elif isinstance(e, Basic.Integer):
                coeff, rest = b.as_coeff_terms()
                l = [s**e for s in rest]
                if e.is_negative:
                    l.reverse()
                return coeff**e * Mul(*l)

        c,t = b.as_coeff_terms()
        if e.is_even and isinstance(c, Basic.Number) and c < 0:
            return (-c * Basic.Mul(*t)) ** e

        if e.atoms(Basic.Wild):
            return Mul(*[t**e for t in b])

    @property
    def precedence(self):
        coeff, rest = self.as_coeff_terms()
        if coeff.is_negative: return Basic.Add_precedence
        return Basic.Mul_precedence

    def tostr(self, level = 0):
        precedence = self.precedence
        coeff, terms = self.as_coeff_terms()
        if coeff.is_negative:
            coeff = -coeff
            if not isinstance(coeff, Basic.One):
                terms.insert(0, coeff)
            r = '-' + '*'.join([t.tostr(precedence) for t in terms])
        else:
            r = '*'.join([t.tostr(precedence) for t in self])
        r = r.replace('*1/', '/')
        if precedence <= level:
            return '(%s)' % r
        return r

        numer,denom = self.as_numer_denom()
        if isinstance(denom, Basic.One):
            delim = '*'
            coeff, rest = self.as_coeff_terms()
            r = delim.join([s.tostr(precedence) for s in rest])
            if isinstance(coeff, Basic.One):
                pass
            elif isinstance(-coeff, Basic.One):
                r = '-' + r
            elif coeff.is_negative:
                r = '-' + (-coeff).tostr() + delim + r
            else:
                r = coeff.tostr() + delim + r
        else:
            if len(denom[:])>1:
                r = '(' + numer.tostr() + ') / (' + denom.tostr() + ')'
            else:
                r = '(' + numer.tostr() + ') / ' + denom.tostr()
        if precedence<=level:
            return '(%s)' % r
        return r

    @cache_it
    def as_coeff_terms(self, x=None):
        if x is not None:
            l1 = []
            l2 = []
            for f in self:
                if f.has(x):
                    l2.append(f)
                else:
                    l1.append(f)
            return Mul(*l1), l2
        coeff = self[0]
        if isinstance(coeff, Basic.Number):
            return coeff, list(self[1:])
        return Basic.One(), list(self[:])

    def _eval_expand(self):
        """
        (a + b + ..) * c -> a * c + b * c + ..
        """
        seq = [Basic.One()]
        for t in self:
            t = t.expand()
            if isinstance(t, Basic.Add):
                seq = [f1*f2 for f1 in seq for f2 in t]
            else:
                seq = [f*t for f in seq]
        return Basic.Add(*seq, **self._assumptions)

    def _eval_derivative(self, s):
        terms = list(self)
        factors = []
        for i in xrange(len(terms)):
            t = terms[i].diff(s)
            if isinstance(t, Basic.Zero):
                continue
            factors.append(Mul(*(terms[:i]+[t]+terms[i+1:])))
        return Basic.Add(*factors)

    def _matches_simple(pattern, expr, repl_dict):
        # handle (w*3).matches('x*5') -> {w: x*5/3}
        coeff, terms = pattern.as_coeff_terms()
        if len(terms)==1:
            return terms[0].matches(expr / coeff, repl_dict)
        return

    def matches(pattern, expr, repl_dict={}, evaluate=False):
        expr = Basic.sympify(expr)
        if pattern.is_commutative and expr.is_commutative:
            return AssocOp._matches_commutative(pattern, expr, repl_dict, evaluate)
        # todo for commutative parts, until then use the default matches method for non-commutative products
        return Basic.matches(pattern, expr, repl_dict, evaluate)

    @staticmethod
    def _combine_inverse(lhs, rhs):
        if lhs == rhs:
            return Basic.One()
        return lhs / rhs

    def as_numer_denom(self):
        numers,denoms = [],[]
        for t in self:
            if isinstance(t, Basic.Number):
                numers.append(t)
                continue
            n,d = t.as_numer_denom()
            numers.append(n)
            denoms.append(d)
        return Mul(*numers), Mul(*denoms)

    @cache_it_immutable
    def count_ops(self, symbolic=True):
        if symbolic:
            return Basic.Add(*[t.count_ops(symbolic) for t in self[:]]) + Basic.Symbol('MUL') * (len(self[:])-1)
        return Basic.Add(*[t.count_ops(symbolic) for t in self[:]]) + (len(self[:])-1)

    def _eval_integral(self, s):
        coeffs = []
        terms = []
        for t in self:
            if not t.has(s): coeffs.append(t)
            else: terms.append(t)
        if coeffs:
            return Mul(*coeffs) * Mul(*terms).integral(s)
        u = self[0].integral(s)
        v = Mul(*(self[1:]))
        uv = u * v
        return uv - (u*v.diff(s)).integral(s)

    def _eval_defined_integral(self, s, a, b):
        coeffs = []
        terms = []
        for t in self:
            if not t.has(s): coeffs.append(t)
            else: terms.append(t)
        if coeffs:
            return Mul(*coeffs) * Mul(*terms).integral(s==[a,b])
        # (u'v) -> (uv) - (uv')
        u = self[0].integral(s)
        v = Mul(*(self[1:]))
        uv = u * v
        return (uv.subs(s,b) - uv.subs(s,a)) - (u*v.diff(s)).integral(s==[a,b])

    def _eval_is_polynomial(self, syms):
        for term in self:
            if not term._eval_is_polynomial(syms):
                return False
        return True

    _eval_is_real = lambda self: self._eval_template_is_attr('is_real')
    _eval_is_bounded = lambda self: self._eval_template_is_attr('is_bounded')
    _eval_is_commutative = lambda self: self._eval_template_is_attr('is_commutative')
    _eval_is_integer = lambda self: self._eval_template_is_attr('is_integer')
    _eval_is_comparable = lambda self: self._eval_template_is_attr('is_comparable')

    def _eval_is_irrational(self):
        for t in self:
            a = t.is_irrational
            if a: return True
            if a is None: return
        return False

    def _eval_is_positive(self):
        terms = [t for t in self if not t.is_positive]
        if not terms:
            return True
        c = terms[0]
        if len(terms)==1:
            if c.is_nonpositive:
                return False
            return
        r = Mul(*terms[1:])
        if c.is_negative and r.is_negative:
            return True
        if r.is_negative and c.is_negative:
            return True
        # check for nonpositivity, <=0
        if c.is_negative and r.is_nonnegative:
            return False
        if r.is_negative and c.is_nonnegative:
            return False
        if c.is_nonnegative and r.is_nonpositive:
            return False
        if r.is_nonnegative and c.is_nonpositive:
            return False


    def _eval_is_negative(self):
        terms = [t for t in self if not t.is_positive]
        if not terms:
            return None
        c = terms[0]
        if len(terms)==1:
            return c.is_negative
        r = Mul(*terms[1:])
        # check for nonnegativity, >=0
        if c.is_negative and r.is_nonpositive:
            return False
        if r.is_negative and c.is_nonpositive:
            return False
        if c.is_nonpositive and r.is_nonpositive:
            return False
        if c.is_nonnegative and r.is_nonnegative:
            return False

    def _eval_is_odd(self):
        if self.is_integer:
            r = True
            for t in self:
                if t.is_even:
                    return False
                if t.is_odd is None:
                    r = None
            return r

    def _eval_subs(self, old, new):
        if self==old:
            return new
        coeff1,terms1 = self.as_coeff_terms()
        coeff2,terms2 = old.as_coeff_terms()
        if terms1==terms2: # (2*a).subs(3*a,y) -> 2/3*y
            return new * coeff1/coeff2
        l1,l2 = len(terms1),len(terms2)
        if l2<l1: # (a*b*c*d).subs(b*c,x) -> a*x*d
            for i in xrange(l1-l2+1):
                if terms2==terms1[i:i+l2]:
                    m1 = Mul(*terms1[:i]).subs(old,new)
                    m2 = Mul(*terms1[i+l2:]).subs(old,new)
                    return Mul(*([coeff1/coeff2,m1,new,m2]))
        return self.__class__(*[s.subs(old, new) for s in self])

    def _eval_oseries(self, order):
        x = order.symbols[0]
        l = []
        r = []
        lt = []
        for t in self:
            if not t.has(x):
                l.append(t)
                continue
            r.append(t)
            lt.append(t.as_leading_term(x))
        if not r:
            if order.contains(1,x): return S.Zero
            return Mul(*l)
        if len(r)==1:
            return Mul(*(l + [r[0].oseries(order)]))
        for i in xrange(len(r)):
            m = Mul(*(lt[:i]+lt[i+1:]))
            o = order/m
            l.append(r[i].oseries(o))
        return Mul(*l)

    def _eval_as_leading_term(self, x):
        return Mul(*[t.as_leading_term(x) for t in self])

    def _eval_conjugate(self):
        return Mul(*[t.conjugate() for t in self])

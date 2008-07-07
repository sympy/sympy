
from basic import Basic, S, C, sympify
from operations import AssocOp
from methods import RelMeths, ArithMeths
from cache import cacheit

from symbol import Symbol, Wild
# from function import FunctionClass, WildFunction /cyclic/
# from numbers import Number, Integer, Real /cyclic/
# from add   import Add /cyclic/
# from power import Pow /cyclic/

class Mul(AssocOp, RelMeths, ArithMeths):

    __slots__ = []

    is_Mul = True

    @classmethod
    def flatten(cls, seq):
        # apply associativity, separate commutative part of seq
        c_part = []
        nc_part = []

        c_seq = []
        nc_seq = seq

        coeff = S.One       # standalone term
                            # e.g. 3 * ...

        c_powers = {}       # base -> exp                      z
                            # e.g. (x+y) -> z  for  ... * (x+y)  * ...

        exp_dict = {}       # num-base -> exp           y
                            # e.g.  3 -> y  for  ... * 3  * ...

        inv_exp_dict = {}   # exp -> Mul(num-bases)     x    x
                            # e.g.  x -> 6  for  ... * 2  * 3  * ...

        lambda_args = None
        order_symbols = None


        while c_seq or nc_seq:

            # COMMUTATIVE
            if c_seq:
                # first process commutative objects
                o = c_seq.pop(0)
                if isinstance(o, FunctionClass):
                    if o.nargs is not None:
                        o, lambda_args = o.with_dummy_arguments(lambda_args)

                # O(x)
                if o.is_Order:
                    o, order_symbols = o.as_expr_symbols(order_symbols)

                # Mul([...])
                if o.is_Mul:
                    # associativity
                    c_seq = list(o.args[:]) + c_seq
                    continue

                # 3
                if o.is_Number:
                    coeff *= o
                    continue

                #  y
                # x
                if o.is_Pow:
                    base, exponent = o.as_base_exp()

                    #  y
                    # 3
                    if base.is_Number:

                        # let's collect factors with numeric base
                        if base in exp_dict:
                            exp_dict[base] += exponent
                        else:
                            exp_dict[base] = exponent
                        continue

                # exp(x)
                if o.func is C.exp:
                    # exp(x) / exp(y) -> exp(x-y)
                    b = S.Exp1
                    e = o.args[0]

                # everything else
                else:
                    b, e = o.as_base_exp()

                # now we have
                # o = b**e

                #         n          n          n
                # (-3 + y)   ->  (-1)  * (3 - y)
                if b.is_Add and e.is_Number:
                    #found factor (x+y)**number; split off initial coefficient
                    c, t = b.as_coeff_terms()
                    #last time I checked, Add.as_coeff_terms returns One or NegativeOne
                    #but this might change
                    if c.is_negative and not e.is_integer:
                        # extracting root from negative number: ignore sign
                        if c is not S.NegativeOne:
                            # make c positive (probably never occurs)
                            coeff *= (-c) ** e
                            assert len(t)==1,`t`
                            b = -t[0]
                        #else: ignoring sign from NegativeOne: nothing to do!
                    elif c is not S.One:
                        coeff *= c ** e
                        assert len(t)==1,`t`
                        b = t[0]
                    #else: c is One, so pass

                # let's collect factors with the same base, so e.g.
                #  y    z     y+z
                # x  * x  -> x
                if b in c_powers:
                    c_powers[b] += e
                else:
                    c_powers[b] = e


            # NON-COMMUTATIVE
            else:
                o = nc_seq.pop(0)
                if isinstance(o, WildFunction):
                    pass
                elif isinstance(o, FunctionClass):
                    if o.nargs is not None:
                        o, lambda_args = o.with_dummy_arguments(lambda_args)
                elif o.is_Order:
                    o, order_symbols = o.as_expr_symbols(order_symbols)

                # -> commutative
                if o.is_commutative:
                    # separate commutative symbols
                    c_seq.append(o)
                    continue

                # Mul([...])
                if o.__class__ is cls:
                    # associativity
                    nc_seq = list(o.args) + nc_seq
                    continue
                if not nc_part:
                    nc_part.append(o)
                    continue

                #                             b    c       b+c
                # try to combine last terms: a  * a   ->  a
                o1 = nc_part.pop()
                b1,e1 = o1.as_base_exp()
                b2,e2 = o.as_base_exp()
                if b1==b2:
                    nc_seq.insert(0, b1 ** (e1 + e2))
                else:
                    nc_part.append(o1)
                    nc_part.append(o)


        # ................................
        # now we have:
        # - coeff:
        # - c_powers:   b -> e
        # - exp_dict:   3 -> e

        # XXX
        for b, e in c_powers.items():
            if e is S.Zero:
                continue

            if e is S.One:
                if b.is_Number:
                    coeff *= b
                else:
                    c_part.append(b)
            elif e.is_Integer and b.is_Number:
                coeff *= b ** e
            else:
                c_part.append(Pow(b, e))

        #  x    x     x
        # 2  * 3  -> 6
        for b,e in exp_dict.items():
            if e in inv_exp_dict:
                inv_exp_dict[e] *= b
            else:
                inv_exp_dict[e] = b

        for e,b in inv_exp_dict.items():
            if e is S.Zero:
                continue

            if e is S.One:
                if b.is_Number:
                    coeff *= b
                else:
                    c_part.append(b)
            elif e.is_Integer and b.is_Number:
                coeff *= b ** e
            else:
                obj = b**e
                if obj.is_Number:
                    coeff *= obj
                else:
                    c_part.append(obj)


        # deal with
        # (oo|nan|zero) * ...
        if (coeff is S.Infinity) or (coeff is S.NegativeInfinity):
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
        elif (coeff is S.Zero) or (coeff is S.NaN):
            c_part, nc_part = [coeff], []
        elif coeff.is_Real:
            if coeff == Real(0):
                c_part, nc_part = [coeff], []
            elif coeff != Real(1):
                c_part.insert(0, coeff)
        elif coeff is not S.One:
            c_part.insert(0, coeff)

        # order commutative part canonically
        c_part.sort(Basic.compare)

        # we are done
        if len(c_part)==2 and c_part[0].is_Number and c_part[1].is_Add:
            # 2*(1+a) -> 2 + 2 * a
            coeff = c_part[0]
            c_part = [Add(*[coeff*f for f in c_part[1].args])]

        return c_part, nc_part, lambda_args, order_symbols


    def _eval_power(b, e):
        if e.is_Number:
            if b.is_commutative:
                if e.is_Integer:
                    # (a*b)**2 -> a**2 * b**2
                    return Mul(*[s**e for s in b.args])

                if e.is_rational:
                    coeff, rest = b.as_coeff_terms()
                    if coeff == -1:
                        return None
                    elif coeff < 0:
                        return (-coeff)**e * Mul(*((S.NegativeOne,) +rest))**e
                    else:
                        return coeff**e * Mul(*[s**e for s in rest])


                coeff, rest = b.as_coeff_terms()
                if coeff is not S.One:
                    # (2*a)**3 -> 2**3 * a**3
                    return coeff**e * Mul(*[s**e for s in rest])
            elif e.is_Integer:
                coeff, rest = b.as_coeff_terms()
                l = [s**e for s in rest]
                if e.is_negative:
                    l.reverse()
                return coeff**e * Mul(*l)

        c,t = b.as_coeff_terms()
        if e.is_even and c.is_Number and c < 0:
            return (-c * Mul(*t)) ** e

        #if e.atoms(Wild):
        #    return Mul(*[t**e for t in b])

    def _eval_evalf(self):
        return self.__class__(*[s.evalf() for s in self.args]).expand()

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
            if coeff is not S.One:
                terms = (coeff,) + terms
            if isinstance(terms, Basic):
                terms = terms.args
            r = '-' + '*'.join([t.tostr(precedence) for t in terms])
        else:
            r = '*'.join([t.tostr(precedence) for t in self.args])
        r = r.replace('*1/', '/')
        if precedence <= level:
            return '(%s)' % r
        return r

        numer,denom = self.as_numer_denom()
        if denom is S.One:
            delim = '*'
            coeff, rest = self.as_coeff_terms()
            r = delim.join([s.tostr(precedence) for s in rest.args])
            if coeff is S.One:
                pass
            elif -coeff is S.One:
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

    @cacheit
    def as_two_terms(self):
        if len(self.args) == 1:
            return S.One, self
        return self.args[0], Mul(*self.args[1:])

    @cacheit
    def as_coeff_terms(self, x=None):
        if x is not None:
            l1 = []
            l2 = []
            for f in self.args:
                if f.has(x):
                    l2.append(f)
                else:
                    l1.append(f)
            return Mul(*l1), tuple(l2)
        coeff = self.args[0]
        if coeff.is_Number:
            return coeff, self.args[1:]
        return S.One, self.args

    @staticmethod
    def _expandsums(sums):
        L = len(sums)
        if len(sums) == 1:
            return sums[0]
        terms = []
        left = Mul._expandsums(sums[:L//2])
        right = Mul._expandsums(sums[L//2:])
        if isinstance(right, Basic):
            right = right.args
        if isinstance(left, Basic):
            left = left.args

        if len(left) == 1 and len(right) == 1:
            # no expansion needed, bail out now to avoid infinite recursion
            return [Mul(left[0], right[0])]

        terms = []
        for a in left:
            for b in right:
                terms.append(Mul(a,b).expand())
            added = Add(*terms)
            if added.is_Add:
                terms = list(added.args)
            else:
                terms = [added]
        return terms

    def _eval_expand_basic(self):
        plain, sums, rewrite = [], [], False

        for factor in self.args:
            terms = factor._eval_expand_basic()

            if terms is not None:
                factor = terms

            if factor.is_Add:
                sums.append(factor)
                rewrite = True
            else:
                if factor.is_commutative:
                    plain.append(factor)
                else:
                    sums.append([factor])

                if terms is not None:
                    rewrite = True

        if not rewrite:
            return None
        else:
            if sums:
                terms = Mul._expandsums(sums)

                if isinstance(terms, Basic):
                    terms = terms.args

                plain = Mul(*plain)

                return Add(*(Mul(plain, term) for term in terms), **self._assumptions)
            else:
                return Mul(*plain, **self._assumptions)

    def _eval_derivative(self, s):
        terms = list(self.args)
        factors = []
        for i in xrange(len(terms)):
            t = terms[i].diff(s)
            if t is S.Zero:
                continue
            factors.append(Mul(*(terms[:i]+[t]+terms[i+1:])))
        return Add(*factors)

    def _matches_simple(pattern, expr, repl_dict):
        # handle (w*3).matches('x*5') -> {w: x*5/3}
        coeff, terms = pattern.as_coeff_terms()
        if len(terms)==1:
            return terms[0].matches(expr / coeff, repl_dict)
        return

    def matches(pattern, expr, repl_dict={}, evaluate=False):
        expr = sympify(expr)
        if pattern.is_commutative and expr.is_commutative:
            return AssocOp._matches_commutative(pattern, expr, repl_dict, evaluate)
        # todo for commutative parts, until then use the default matches method for non-commutative products
        return Basic.matches(pattern, expr, repl_dict, evaluate)

    @staticmethod
    def _combine_inverse(lhs, rhs):
        if lhs == rhs:
            return S.One
        return lhs / rhs

    def as_powers_dict(self):
        return dict([ term.as_base_exp() for term in self ])

    def as_numer_denom(self):
        numers, denoms = [],[]
        for t in self.args:
            n,d = t.as_numer_denom()
            numers.append(n)
            denoms.append(d)
        return Mul(*numers), Mul(*denoms)

    @cacheit
    def count_ops(self, symbolic=True):
        if symbolic:
            return Add(*[t.count_ops(symbolic) for t in self[:]]) + Symbol('MUL') * (len(self[:])-1)
        return Add(*[t.count_ops(symbolic) for t in self.args[:]]) + (len(self.args)-1)

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
        for term in self.args:
            if not term._eval_is_polynomial(syms):
                return False
        return True

    _eval_is_real = lambda self: self._eval_template_is_attr('is_real')
    _eval_is_bounded = lambda self: self._eval_template_is_attr('is_bounded')
    _eval_is_commutative = lambda self: self._eval_template_is_attr('is_commutative')
    _eval_is_integer = lambda self: self._eval_template_is_attr('is_integer')
    _eval_is_comparable = lambda self: self._eval_template_is_attr('is_comparable')

    def _eval_is_irrational(self):
        for t in self.args:
            a = t.is_irrational
            if a: return True
            if a is None: return
        return False

    def _eval_is_positive(self):
        terms = [t for t in self.args if not t.is_positive]
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
        terms = [t for t in self.args if not t.is_positive]
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
            for t in self.args:
                if t.is_even:
                    return False
                if t.is_odd is None:
                    r = None
            return r

    def _eval_subs(self, old, new):
        if self==old:
            return new
        if isinstance(old, FunctionClass):
            return self.__class__(*[s.subs(old, new) for s in self.args ])
        coeff1,terms1 = self.as_coeff_terms()
        coeff2,terms2 = old.as_coeff_terms()
        if terms1==terms2: # (2*a).subs(3*a,y) -> 2/3*y
            return new * coeff1/coeff2
        l1,l2 = len(terms1),len(terms2)
        if l2 == 0:
            # if old is just a number, go through the self.args one by one
            return Mul(*[x.subs(old, new) for x in self.args])
        elif l2<l1:
            # old is some something more complex, like:
            # (a*b*c*d).subs(b*c,x) -> a*x*d
            # then we need to search where in self.args the "old" is, and then
            # correctly substitute both terms and coefficients.
            for i in xrange(l1-l2+1):
                if terms2==terms1[i:i+l2]:
                    m1 = Mul(*terms1[:i]).subs(old,new)
                    m2 = Mul(*terms1[i+l2:]).subs(old,new)
                    return Mul(*([coeff1/coeff2,m1,new,m2]))
        return self.__class__(*[s.subs(old, new) for s in self.args])

    def _eval_oseries(self, order):
        x = order.symbols[0]
        l = []
        r = []
        lt = []
        #separate terms containing "x" (r) and the rest (l)
        for t in self.args:
            if not t.has(x):
                l.append(t)
                continue
            r.append(t)
        #if r is empty or just one term, it's easy:
        if not r:
            if order.contains(1,x): return S.Zero
            return Mul(*l)
        if len(r)==1:
            return Mul(*(l + [r[0].oseries(order)]))
        #otherwise, we need to calculate how many orders we need to calculate
        #in each term. Currently this is done using as_leading_term, but this
        #is fragile and slow, because this involves limits. Let's find some
        #more clever approach.
        lt = [t.as_leading_term(x) for t in r]
        for i in xrange(len(r)):
            m = Mul(*(lt[:i]+lt[i+1:]))
            #calculate how many orders we want
            o = order/m
            #expand each term and multiply things together
            l.append(r[i].oseries(o))
        #shouldn't we rather expand everything? This seems to me to leave
        #things as (x+x**2+...)*(x-x**2+...) etc.:
        return Mul(*l)

    def nseries(self, x, x0, n):
        terms = [t.nseries(x, x0, n) for t in self.args]
        return Mul(*terms).expand()


    def _eval_as_leading_term(self, x):
        return Mul(*[t.as_leading_term(x) for t in self.args])

    def _eval_conjugate(self):
        return Mul(*[t.conjugate() for t in self.args])

    def _sage_(self):
        s = 1
        for x in self.args:
            s *= x._sage_()
        return s


# /cyclic/
import basic as _
_.Mul       = Mul
del _

import methods as _
_.Mul       = Mul
del _

import add as _
_.Mul       = Mul
del _

import power as _
_.Mul       = Mul
del _

import numbers as _
_.Mul       = Mul
del _

import operations as _
_.Mul       = Mul
del _

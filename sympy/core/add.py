
from basic import Basic, S, C
from operations import AssocOp
from methods import RelMeths, ArithMeths
from cache import cacheit

from symbol import Symbol, Wild, Temporary
# from numbers import Number    /cyclic/
# from mul import Mul    /cyclic/
# from function import FunctionClass    /cyclic/

class Add(AssocOp, RelMeths, ArithMeths):

    precedence = Basic.Add_precedence

    __slots__ = []

    is_Add = True

    @classmethod
    def flatten(cls, seq):
        """
        Takes the sequence "seq" of nested Adds and returns a flatten list.

        Returns: (commutative_part, noncommutative_part, lambda_args,
            order_symbols)

        Applies associativity, all terms are commutable with respect to
        addition.
        """
        terms = {}      # term -> coeff
                        # e.g. x**2 -> 5   for ... + 5*x**2 + ...

        coeff = S.Zero  # standalone term
                        # e.g. 3 + ...
        lambda_args = None
        order_factors = []
        while seq:
            o = seq.pop(0)
            #if o.is_Function:
            #    if o.nargs is not None:
            #        o, lambda_args = o.with_dummy_arguments(lambda_args)

            # O(x)
            if o.is_Order:
                for o1 in order_factors:
                    if o1.contains(o):
                        o = None
                        break
                if o is None:
                    continue
                order_factors = [o]+[o1 for o1 in order_factors if not o.contains(o1)]
                continue

            # 3
            if o.is_Number:
                coeff += o
                continue

            # Add([...])
            if o.is_Add:
                seq = list(o.args) + seq
                continue

            # Mul([...])
            if o.is_Mul:
                c = o.args[0]

                # 3*...
                if c.is_Number:
                    if c is S.One:
                        s = o
                    else:
                        s = Mul(*o.args[1:])
                else:
                    c = S.One
                    s = o

            # everything else
            else:
                c = S.One
                s = o


            # now we have:
            # o = c*s, where
            #
            # c is a Number
            # s is an expression with number factor extracted

            # let's collect terms with the same s, so e.g.
            # 2*x**2 + 3*x**2  ->  5*x**2
            if terms.has_key(s):
                terms[s] += c
            else:
                terms[s] = c


        # now let's construct new args:
        # [2*x**2, x**3, 7*x**4, pi, ...]
        newseq = []
        noncommutative = False
        for s,c in terms.items():
            # 0*s
            if c is S.Zero:
                continue
            # 1*s
            elif c is S.One:
                newseq.append(s)
            # c*s
            else:
                newseq.append(Mul(c,s))
            noncommutative = noncommutative or not s.is_commutative

        # deal with nan, oo, -oo, etc...
        if coeff is S.NaN:
            newseq = [coeff]
        elif (coeff is S.Infinity) or (coeff is S.NegativeInfinity):
            newseq = [coeff] + [f for f in newseq if not f.is_real]
        elif coeff is not S.Zero:
            newseq.insert(0, coeff)

        # process O(x)
        if order_factors:
            newseq2 = []
            for t in newseq:
                for o in order_factors:
                    if o.contains(t):
                        t = None
                        break
                if t is not None:
                    newseq2.append(t)
            newseq = newseq2 + order_factors

        # order args canonically
        # Currently we sort things using hashes, as it is quite fast. A better
        # solution is not to sort things at all - but this needs some more
        # fixing.
        newseq.sort(key=hash)

        # we are done
        if noncommutative:
            return [],newseq,lambda_args,None
        else:
            return newseq,[],lambda_args,None

    @staticmethod
    def compare_terms(a, b):
        """
        Is a>b in the sense of ordering in printing?

        yes ..... return 1
        no ...... return -1
        equal ... return 0

        Strategy:

        It uses Basic.compare as a fallback, but improves it in many cases,
        like x**3, x**4, O(x**3) etc. In those simple cases, it just parses the
        expression and returns the "sane" ordering such as:

        1 < x < x**2 < x**3 < O(x**4) etc.

        """
        from sympy.series.order import Order
        if isinstance(a, Order) and not isinstance(b, Order):
            return 1
        if not isinstance(a, Order) and isinstance(b, Order):
            return -1
        p1, p2, p3 = Wild("p1"), Wild("p2"), Wild("p3")
        r_a = a.match(p1 * p2**p3)
        r_b = b.match(p1 * p2**p3)
        if r_a is not None and r_b is not None:
            c = Basic.compare(r_a[p3], r_b[p3])
            if c!=0:
                return c
        return Basic.compare(a,b)

    def tostr(self, level=0):
        coeff, rest = self.as_coeff_factors()

        # Now we need to sort the factors in Add, which are in "rest". Any
        # ordering is fine, but some ordering looks better and some looks bad.

        # This particular solution is slow, but it ensures a sane ordering. It
        # can of course be improved:
        rest = list(rest)
        rest.sort(self.compare_terms)

        l = []
        precedence = self.precedence
        if coeff is not S.Zero:
            l.append(coeff.tostr(precedence))
        for factor in rest:
            f = factor.tostr()
            if f.startswith('-'):
                if l:
                    l.extend(['-',f[1:]])
                else:
                    l.append(f)
            else:
                l.extend(['+',f])
        if l[0]=='+': l.pop(0)
        r = ' '.join(l)
        if precedence<=level:
            return '(%s)' % r
        return r

    @cacheit
    def as_coeff_factors(self, x=None):
        if x is not None:
            l1 = []
            l2 = []
            for f in self.args:
                if f.has(x):
                    l2.append(f)
                else:
                    l1.append(f)
            return Add(*l1), tuple(l2)
        coeff = self.args[0]
        if coeff.is_Number:
            return coeff, self.args[1:]
        return S.Zero, self.args

    def _eval_derivative(self, s):
        return Add(*[f.diff(s) for f in self.args])

    def nseries(self, x, x0, n):
        terms = [t.nseries(x, x0, n) for t in self.args]
        return Add(*terms)

    def _matches_simple(pattern, expr, repl_dict):
        # handle (w+3).matches('x+5') -> {w: x+2}
        coeff, factors = pattern.as_coeff_factors()
        if len(factors)==1:
            return factors[0].matches(expr - coeff, repl_dict)
        return

    matches = AssocOp._matches_commutative

    @staticmethod
    def _combine_inverse(lhs, rhs):
        return lhs - rhs

    @cacheit
    def as_two_terms(self):
        if len(self.args) == 1:
            return S.Zero, self
        return self.args[0], Add(*self.args[1:])

    def as_numer_denom(self):
        numers, denoms = [],[]
        for n,d in [f.as_numer_denom() for f in self.args]:
            numers.append(n)
            denoms.append(d)
        r = xrange(len(numers))
        return Add(*[Mul(*(denoms[:i]+[numers[i]]+denoms[i+1:])) for i in r]),Mul(*denoms)

    def _calc_splitter(self, d):
        if d.has_key(self):
            return d[self]
        coeff, factors = self.as_coeff_factors()
        r = self.__class__(*[t._calc_splitter(d) for t in factors])
        if r.is_Add and 0:
            for e,t in d.items():
                w = Wild('w')
                d1 = r.match(e+w)
                if d1 is not None:
                    r1 = t + d1[w]
                    break
        if d.has_key(r):
            return coeff + d[r]
        s = d[r] = Temporary()
        return s + coeff

    def count_ops(self, symbolic=True):
        if symbolic:
            return Add(*[t.count_ops(symbolic) for t in self[:]]) + Symbol('ADD') * (len(self[:])-1)
        return Add(*[t.count_ops(symbolic) for t in self.args[:]]) + (len(self.args[:])-1)

    def _eval_integral(self, s):
        return Add(*[f.integral(s) for f in self])

    def _eval_defined_integral(self, s,a,b):
        return Add(*[f.integral(s==[a,b]) for f in self])

    def _eval_is_polynomial(self, syms):
        for term in self.args:
            if not term._eval_is_polynomial(syms):
                return False
        return True

    # assumption methods
    _eval_is_real = lambda self: self._eval_template_is_attr('is_real')
    _eval_is_bounded = lambda self: self._eval_template_is_attr('is_bounded')
    _eval_is_commutative = lambda self: self._eval_template_is_attr('is_commutative')
    _eval_is_integer = lambda self: self._eval_template_is_attr('is_integer')
    _eval_is_comparable = lambda self: self._eval_template_is_attr('is_comparable')

    def _eval_is_odd(self):
        l = [f for f in self.args if not (f.is_even==True)]
        if not l:
            return False
        if l[0].is_odd:
            return Add(*l[1:]).is_even

    def _eval_is_irrational(self):
        for t in self:
            a = t.is_irrational
            if a: return True
            if a is None: return
        return False

    def _eval_is_positive(self):
        c = self.args[0]
        r = Add(*self.args[1:])
        if c.is_positive and r.is_positive:
            return True
        if c.is_unbounded:
            if r.is_unbounded:
                # either c or r is negative
                return
            else:
                return c.is_positive
        elif r.is_unbounded:
            return r.is_positive
        if c.is_nonnegative and r.is_positive:
            return True
        if r.is_nonnegative and c.is_positive:
            return True
        if c.is_nonpositive and r.is_nonpositive:
            return False

    def _eval_is_negative(self):
        c = self.args[0]
        r = Add(*self.args[1:])
        if c.is_negative and r.is_negative:
            return True
        if c.is_unbounded:
            if r.is_unbounded:
                # either c or r is positive
                return
            else:
                return c.is_negative
        elif r.is_unbounded:
            return r.is_negative
        if c.is_nonpositive and r.is_negative:
            return True
        if r.is_nonpositive and c.is_negative:
            return True
        if c.is_nonnegative and r.is_nonnegative:
            return False

    def as_coeff_terms(self, x=None):
        # -2 + 2 * a -> -1, 2-2*a
        if self.args[0].is_Number and self.args[0].is_negative:
            return -S.One,(-self,)
        return S.One,(self,)

    def _eval_subs(self, old, new):
        if self==old: return new
        if isinstance(old, FunctionClass):
            return self.__class__(*[s.subs(old, new) for s in self.args ])
        coeff1,factors1 = self.as_coeff_factors()
        coeff2,factors2 = old.as_coeff_factors()
        if factors1==factors2: # (2+a).subs(3+a,y) -> 2-3+y
            return new + coeff1 - coeff2
        if old.is_Add:
            l1,l2 = len(factors1),len(factors2)
            if l2<l1: # (a+b+c+d).subs(b+c,x) -> a+x+d
                for i in xrange(l1-l2+1):
                    if factors2==factors1[i:i+l2]:
                        factors1 = list(factors1)
                        factors2 = list(factors2)
                        return Add(*([coeff1-coeff2]+factors1[:i]+[new]+factors1[i+l2:]))
        return self.__class__(*[s.subs(old, new) for s in self.args])

    def _eval_oseries(self, order):
        return Add(*[f.oseries(order) for f in self.args])

    @cacheit
    def extract_leading_order(self, *symbols):
        lst = []
        seq = [(f, C.Order(f, *symbols)) for f in self.args]
        for ef,of in seq:
            for e,o in lst:
                if o.contains(of):
                    of = None
                    break
            if of is None:
                continue
            new_lst = [(ef,of)]
            for e,o in lst:
                if of.contains(o):
                    continue
                new_lst.append((e,o))
            lst = new_lst
        return tuple(lst)

    def _eval_as_leading_term(self, x):
        coeff, factors = self.as_coeff_factors(x)
        has_unbounded = bool([f for f in self.args if f.is_unbounded])
        if has_unbounded:
            if isinstance(factors, Basic):
                factors = factors.args
            factors = [f for f in factors if not f.is_bounded]
        if coeff is not S.Zero:
            o = C.Order(x)
        else:
            o = C.Order(factors[0]*x,x)
        s = self.oseries(o)
        while s is S.Zero:
            o *= x
            s = self.oseries(o)
        if s.is_Add:
            lst = s.extract_leading_order(x)
            return Add(*[e for (e,f) in lst])
        return s.as_leading_term(x)

    def _eval_conjugate(self):
        return Add(*[t.conjugate() for t in self.args])

    def __neg__(self):
        return Add(*[-t for t in self.args])

    def _sage_(self):
        s = 0
        for x in self.args:
            s += x._sage_()
        return s


# /cyclic/
import basic as _
_.Add       = Add
del _

import methods as _
_.Add       = Add
del _

import mul as _
_.Add       = Add
del _

import power as _
_.Add       = Add
del _

import operations as _
_.Add       = Add
del _

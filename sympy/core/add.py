from basic import Basic, S, C
from operations import AssocOp
from cache import cacheit
from symbol import Symbol

class Add(AssocOp):

    __slots__ = []

    is_Add = True

    @classmethod
    def flatten(cls, seq):
        """
        Takes the sequence "seq" of nested Adds and returns a flatten list.

        Returns: (commutative_part, noncommutative_part, order_symbols)

        Applies associativity, all terms are commutable with respect to
        addition.
        """
        terms = {}      # term -> coeff
                        # e.g. x**2 -> 5   for ... + 5*x**2 + ...

        coeff = S.Zero  # standalone term
                        # e.g. 3 + ...
        order_factors = []

        for o in seq:

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
            elif o.is_Number:
                coeff += o
                continue

            # Add([...])
            elif o.is_Add:
                # NB: here we assume Add is always commutative
                seq.extend(o.args)  # TODO zerocopy?
                continue

            # Mul([...])
            elif o.is_Mul:
                c = o.args[0]

                # 3*...
                if c.is_Number:
                    if c is S.One:
                        s = o
                    else:
                        s = o.as_two_terms()[1]

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
            if s in terms:
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
                if s.is_Mul:
                    # Mul, already keeps it's arguments in perfect order.
                    # so we can simply put c in slot0 and go the fast way.
                    cs = s._new_rawargs(*((c,) + s.args))
                    newseq.append(cs)

                else:
                    # alternatively we have to call all Mul's machinery (slow)
                    newseq.append(Mul(c,s))

            noncommutative = noncommutative or not s.is_commutative

        # nan
        if coeff is S.NaN:
            # we know for sure the result will be nan
            return [S.NaN], [], None

        # oo, -oo
        elif (coeff is S.Infinity) or (coeff is S.NegativeInfinity):
            newseq = [f for f in newseq if not f.is_real]


        # process O(x)
        if order_factors:
            newseq2 = []
            for t in newseq:
                for o in order_factors:
                    # x + O(x) -> O(x)
                    if o.contains(t):
                        t = None
                        break
                # x + O(x**2) -> x + O(x**2)
                if t is not None:
                    newseq2.append(t)
            newseq = newseq2 + order_factors
            # 1 + O(1) -> O(1)
            for o in order_factors:
                if o.contains(coeff):
                    coeff = S.Zero
                    break


        # order args canonically
        # Currently we sort things using hashes, as it is quite fast. A better
        # solution is not to sort things at all - but this needs some more
        # fixing.
        newseq.sort(key=hash)

        # current code expects coeff to be always in slot-0
        if coeff is not S.Zero:
            newseq.insert(0, coeff)

        # we are done
        if noncommutative:
            return [], newseq, None
        else:
            return newseq, [], None


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

    def _eval_nseries(self, x, x0, n):
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
        """
        Returns lhs - rhs, but treats arguments like symbols, so things like
        oo - oo return 0, instead of a nan.
        """
        from sympy import oo, I, expand_mul
        if lhs == oo and rhs == oo or lhs == oo*I and rhs == oo*I:
            return S.Zero
        return expand_mul(lhs - rhs)

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

    def count_ops(self, symbolic=True):
        if symbolic:
            return Add(*[t.count_ops(symbolic) for t in self.args]) + \
                Symbol('ADD') * (len(self.args) - 1)
        return Add(*[t.count_ops(symbolic) for t in self.args]) + \
            (len(self.args) - 1)

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
        for t in self.args:
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
        if self == old: return new
        if isinstance(old, FunctionClass):
            return self.__class__(*[s._eval_subs(old, new) for s in self.args ])
        coeff_self, factors_self = self.as_coeff_factors()
        coeff_old, factors_old = old.as_coeff_factors()
        if factors_self == factors_old: # (2+a).subs(3+a,y) -> 2-3+y
            return Add(new, coeff_self, -coeff_old)
        if old.is_Add:
            if len(factors_old) < len(factors_self): # (a+b+c+d).subs(b+c,x) -> a+x+d
                self_set = set(factors_self)
                old_set = set(factors_old)
                if old_set < self_set:
                    ret_set = self_set - old_set
                    return Add(new, coeff_self, -coeff_old, *[s._eval_subs(old, new) for s in ret_set])
        return self.__class__(*[s._eval_subs(old, new) for s in self.args])

    @cacheit
    def extract_leading_order(self, *symbols):
        """
        Returns the leading term and it's order.

        Examples:

        >>> (x+1+1/x**5).extract_leading_order(x)
        ((1/x**5, O(1/x**5)),)
        >>> (1+x).extract_leading_order(x)
        ((1, O(1, x)),)
        >>> (x+x**2).extract_leading_order(x)
        ((x, O(x)),)
        """
        lst = []
        seq = [(f, C.Order(f, *symbols)) for f in self.args]
        for ef,of in seq:
            for e,o in lst:
                if o.contains(of) and o != of:
                    of = None
                    break
            if of is None:
                continue
            new_lst = [(ef,of)]
            for e,o in lst:
                if of.contains(o) and o != of:
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
        n = 1
        s = self.nseries(x, 0, n)
        while s.is_Order:
            n +=1
            s = self.nseries(x, 0, n)
        if s.is_Add:
            s = s.removeO()
        if s.is_Add:
            lst = s.extract_leading_order(x)
            return Add(*[e for (e,f) in lst])
        return s.as_leading_term(x)

    def _eval_power(self, other):
        #         n          n          n
        # (-3 + y)   ->  (-1)  * (3 - y)
        # similar to Mul.flatten()
        c, t = self.as_coeff_terms()
        if c.is_negative and not other.is_integer:
            if c is not S.NegativeOne:
                coeff = (-c) ** other
                assert len(t) == 1, 't'
                b = -t[0]
                return coeff*b**other
        elif c is not S.One:
            coeff = c ** other
            assert len(t) == 1, 't'
            b = t[0]
            return coeff*b**other
        return

    def _eval_conjugate(self):
        return Add(*[t.conjugate() for t in self.args])

    def _eval_expand_basic(self, deep=True, **hints):
        sargs, terms = self.args[:], []
        for term in sargs:
            if hasattr(term, '_eval_expand_basic'):
                newterm = term._eval_expand_basic(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

    def _eval_expand_power_exp(self, deep=True, **hints):
        sargs, terms = self.args[:], []
        for term in sargs:
            if hasattr(term, '_eval_expand_power_exp'):
                newterm = term._eval_expand_power_exp(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

    def _eval_expand_power_base(self, deep=True, **hints):
        sargs, terms = self.args[:], []
        for term in sargs:
            if hasattr(term, '_eval_expand_power_base'):
                newterm = term._eval_expand_power_base(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

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
        sargs, terms = self.args[:], []
        for term in sargs:
            if hasattr(term, '_eval_expand_multinomial'):
                newterm = term._eval_expand_multinomial(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

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
        sargs, terms = self.args[:], []
        for term in sargs:
            if hasattr(term, '_eval_expand_complex'):
                newterm = term._eval_expand_complex(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

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

    def __neg__(self):
        return Add(*[-t for t in self.args])

    def _sage_(self):
        s = 0
        for x in self.args:
            s += x._sage_()
        return s

from mul import Mul
from function import FunctionClass

from core import C
from basic import Basic
from singleton import S
from operations import AssocOp
from cache import cacheit
from expr import Expr

class Add(AssocOp):

    __slots__ = []

    is_Add = True

    #identity = S.Zero
    # cyclic import, so defined in numbers.py

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
                    # Mul, already keeps its arguments in perfect order.
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
    def as_coeff_add(self, *deps):
        if deps:
            l1 = []
            l2 = []
            for f in self.args:
                if f.has(*deps):
                    l2.append(f)
                else:
                    l1.append(f)
            return self._new_rawargs(*l1), tuple(l2)
        coeff, notrat = self.args[0].as_coeff_add()
        if not coeff is S.Zero:
            return coeff, notrat + self.args[1:]
        return S.Zero, self.args

    @cacheit
    def as_coeff_mul(self, *deps):
        # -2 + 2 * a -> -1, 2-2*a
        if self.args[0].is_Rational and self.args[0].is_negative:
            return S.NegativeOne, (-self,)
        return Expr.as_coeff_mul(self, *deps)

    def _eval_derivative(self, s):
        return Add(*[f.diff(s) for f in self.args])

    def _eval_nseries(self, x, n):
        terms = [t.nseries(x, n=n) for t in self.args]
        return Add(*terms)

    def _matches_simple(self, expr, repl_dict):
        # handle (w+3).matches('x+5') -> {w: x+2}
        coeff, terms = self.as_coeff_add()
        if len(terms)==1:
            return terms[0].matches(expr - coeff, repl_dict)
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
        """Return head and tail of self.

        This is the most efficient way to get the head and tail of an
        expression.

        - if you want only the head, use self.args[0];
        - if you want to process the arguments of the tail then use
          self.as_coef_add() which gives the head and a tuple containing
          the arguments of the tail when treated as an Add.
        - if you want the coefficient when self is treated as a Mul
          then use self.as_coeff_mul()[0]

        >>> from sympy.abc import x, y
        >>> (3*x*y).as_two_terms()
        (3, x*y)
        """
        if len(self.args) == 1:
            return S.Zero, self
        return self.args[0], self._new_rawargs(*self.args[1:])

    def as_numer_denom(self):
        numers, denoms = [],[]
        for n,d in [f.as_numer_denom() for f in self.args]:
            numers.append(n)
            denoms.append(d)
        r = xrange(len(numers))
        return Add(*[Mul(*(denoms[:i]+[numers[i]]+denoms[i+1:]))
                     for i in r]), Mul(*denoms)

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
            return self._new_rawargs(*l[1:]).is_even

    def _eval_is_irrational(self):
        for t in self.args:
            a = t.is_irrational
            if a: return True
            if a is None: return
        return False

    def _eval_is_positive(self):
        c, r = self.as_two_terms()
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
        c, r = self.as_two_terms()
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

    def _eval_subs(self, old, new):
        if self == old:
            return new
        if isinstance(old, FunctionClass):
            return self.__class__(*[s._eval_subs(old, new) for s in self.args ])
        coeff_self, terms_self = self.as_coeff_add()
        coeff_old, terms_old = old.as_coeff_add()
        if terms_self == terms_old: # (2+a).subs(3+a,y) -> 2-3+y
            return Add(new, coeff_self, -coeff_old)
        if old.is_Add:
            if len(terms_old) < len(terms_self): # (a+b+c+d).subs(b+c,x) -> a+x+d
                self_set = set(terms_self)
                old_set = set(terms_old)
                if old_set < self_set:
                    ret_set = self_set - old_set
                    return Add(new, coeff_self, -coeff_old, *[s._eval_subs(old, new) for s in ret_set])
        return self.__class__(*[s._eval_subs(old, new) for s in self.args])

    @cacheit
    def extract_leading_order(self, *symbols):
        """
        Returns the leading term and it's order.

        Examples:

        >>> from sympy.abc import x
        >>> (x+1+1/x**5).extract_leading_order(x)
        ((x**(-5), O(x**(-5))),)
        >>> (1+x).extract_leading_order(x)
        ((1, O(1)),)
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

    def as_real_imag(self, deep=True):
        sargs, terms = self.args, []
        re_part, im_part = [], []
        for term in sargs:
            re, im = term.as_real_imag(deep=deep)
            re_part.append(re)
            im_part.append(im)
        return (self.new(*re_part), self.new(*im_part))

    def _eval_as_leading_term(self, x):
        coeff, terms = self.as_coeff_add(x)
        has_unbounded = bool([f for f in self.args if f.is_unbounded])
        if has_unbounded:
            if isinstance(terms, Basic):
                terms = terms.args
            terms = [f for f in terms if not f.is_bounded]
        if coeff is not S.Zero:
            o = C.Order(x)
        else:
            o = C.Order(terms[0]*x,x)
        n = 1
        s = self.nseries(x, n=n)
        while s.is_Order:
            n +=1
            s = self.nseries(x, n=n)
        if s.is_Add:
            s = s.removeO()
        if s.is_Add:
            lst = s.extract_leading_order(x)
            return Add(*[e for (e,f) in lst])
        return s.as_leading_term(x)

    def _eval_power(self, other, terms=False):
        #         n          n          n
        # (-3 + y)   ->  (-1)  * (3 - y)
        #
        # If terms=True then return the arguments that should be
        # multiplied together rather than multiplying them.
        #
        # At present, as_coeff_terms return +/-1 but the
        # following should work even if that changes.
        if Basic.keep_sign:
            return None

        rv = None
        c, t = self.as_coeff_mul()
        if c.is_negative and not other.is_integer:
            if c is not S.NegativeOne and self.is_positive:
                coeff = C.Pow(-c, other)
                assert len(t) == 1, 't'
                b = -t[0]
                rv = (coeff, C.Pow(b, other))
        elif c is not S.One:
            coeff = C.Pow(c, other)
            assert len(t) == 1, 't'
            b = t[0]
            rv = (coeff, C.Pow(b, other))
        if not rv or terms:
            return rv
        else:
            return C.Mul(*rv)

    def _eval_conjugate(self):
        return Add(*[t.conjugate() for t in self.args])

    def _eval_expand_basic(self, deep=True, **hints):
        sargs, terms = self.args, []
        for term in sargs:
            if hasattr(term, '_eval_expand_basic'):
                newterm = term._eval_expand_basic(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

    def _eval_expand_power_exp(self, deep=True, **hints):
        sargs, terms = self.args, []
        for term in sargs:
            if hasattr(term, '_eval_expand_power_exp'):
                newterm = term._eval_expand_power_exp(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

    def _eval_expand_power_base(self, deep=True, **hints):
        sargs, terms = self.args, []
        for term in sargs:
            if hasattr(term, '_eval_expand_power_base'):
                newterm = term._eval_expand_power_base(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

    def _eval_expand_mul(self, deep=True, **hints):
        sargs, terms = self.args, []
        for term in sargs:
            if hasattr(term, '_eval_expand_mul'):
                newterm = term._eval_expand_mul(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

    def _eval_expand_multinomial(self, deep=True, **hints):
        sargs, terms = self.args, []
        for term in sargs:
            if hasattr(term, '_eval_expand_multinomial'):
                newterm = term._eval_expand_multinomial(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

    def _eval_expand_log(self, deep=True, **hints):
        sargs, terms = self.args, []
        for term in sargs:
            if hasattr(term, '_eval_expand_log'):
                newterm = term._eval_expand_log(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

    def _eval_expand_complex(self, deep=True, **hints):
        sargs, terms = self.args, []
        for term in sargs:
            if hasattr(term, '_eval_expand_complex'):
                newterm = term._eval_expand_complex(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

    def _eval_expand_trig(self, deep=True, **hints):
        sargs, terms = self.args, []
        for term in sargs:
            if hasattr(term, '_eval_expand_trig'):
                newterm = term._eval_expand_trig(deep=deep, **hints)
            else:
                newterm = term
            terms.append(newterm)
        return self.new(*terms)

    def _eval_expand_func(self, deep=True, **hints):
        sargs, terms = self.args, []
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


    def primitive(self):
        """
        Divide ``self`` by the GCD of coefficients of ``self``.

        Example
        =======

        >>> from sympy.abc import x, y

        >>> (2*x + 4*y).primitive()
        (2, x + 2*y)

        >>> (2*x/3 + 4*y/9).primitive()
        (2/9, 2*y + 3*x)

        >>> (2*x/3 + 4.1*y).primitive()
        (1, 2*x/3 + 4.1*y)

        """
        terms = []
        cont = S.Zero

        for term in self.args:
            coeff = term.as_coeff_mul()[0]

            if coeff.is_Rational:
                cont = cont.gcd(coeff)

                if cont is not S.One:
                    terms.append(term)
                    continue

            return S.One, self

        for i, term in enumerate(terms):
            # XXX: this is extremely slow
            terms[i] = term/cont

        return cont, self._new_rawargs(*terms)

from function import FunctionClass
from mul import Mul
from symbol import Symbol

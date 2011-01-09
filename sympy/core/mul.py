from basic import Basic, C
from singleton import S
from operations import AssocOp
from cache import cacheit
from logic import fuzzy_not
from compatibility import any

# internal marker to indicate:
#   "there are still non-commutative objects -- don't forget to process them"
class NC_Marker:
    is_Order    = False
    is_Mul      = False
    is_Number   = False
    is_Poly     = False

    is_commutative = False


class Mul(AssocOp):

    __slots__ = []

    is_Mul = True

    #identity = S.One
    # cyclic import, so defined in numbers.py

    @classmethod
    def flatten(cls, seq):

        # apply associativity, separate commutative part of seq
        c_part = []         # out: commutative factors
        nc_part = []        # out: non-commutative factors

        nc_seq = []

        coeff = S.One       # standalone term
                            # e.g. 3 * ...

        c_powers = []       # (base,exp)      n
                            # e.g. (x,n) for x

        num_exp = []        # (num-base, exp)           y
                            # e.g.  (3, y)  for  ... * 3  * ...

        order_symbols = None



        # --- PART 1 ---
        #
        # "collect powers and coeff":
        #
        # o coeff
        # o c_powers
        # o num_exp
        #
        # NOTE: this is optimized for all-objects-are-commutative case

        for o in seq:
            # O(x)
            if o.is_Order:
                o, order_symbols = o.as_expr_variables(order_symbols)

            # Mul([...])
            if o.is_Mul:
                if o.is_commutative:
                    seq.extend(o.args)    # XXX zerocopy?

                else:
                    # NCMul can have commutative parts as well
                    for q in o.args:
                        if q.is_commutative:
                            seq.append(q)
                        else:
                            nc_seq.append(q)

                    # append non-commutative marker, so we don't forget to
                    # process scheduled non-commutative objects
                    seq.append(NC_Marker)

                continue

            # 3
            elif o.is_Number:
                coeff *= o
                continue


            elif o.is_commutative:
                #      e
                # o = b
                b, e = o.as_base_exp()

                #  y
                # 3
                if o.is_Pow and b.is_Number:
                    # get all the factors with numeric base so they can be
                    # combined below, but don't combine negatives unless
                    # the exponent is an integer
                    if b.is_positive or e.is_integer:
                        num_exp.append((b, e))
                        continue

                #         n          n          n
                # (-3 + y)   ->  (-1)  * (3 - y)
                #
                # Give powers a chance to become a Mul if that's the
                # behavior obtained from Add._eval_power()
                if not Basic.keep_sign and b.is_Add and e.is_Number:
                    cb = b._eval_power(e, terms=1)
                    if cb:
                        c, b = cb
                        coeff *= c

                c_powers.append((b,e))

            # NON-COMMUTATIVE
            # TODO: Make non-commutative exponents not combine automatically
            else:
                if o is not NC_Marker:
                    nc_seq.append(o)

                # process nc_seq (if any)
                while nc_seq:
                    o = nc_seq.pop(0)
                    if not nc_part:
                        nc_part.append(o)
                        continue

                    #                             b    c       b+c
                    # try to combine last terms: a  * a   ->  a
                    o1 = nc_part.pop()
                    b1,e1 = o1.as_base_exp()
                    b2,e2 = o.as_base_exp()
                    new_exp = e1 + e2
                    # Only allow powers to combine if the new exponent is
                    # not an Add. This allow things like a**2*b**3 == a**5
                    # if a.is_commutative == False, but prohibits
                    # a**x*a**y and x**a*x**b from combining (x,y commute).
                    if b1==b2 and (not new_exp.is_Add):
                        o12 = b1 ** new_exp

                        # now o12 could be a commutative object
                        if o12.is_commutative:
                            seq.append(o12)
                            continue
                        else:
                            nc_seq.insert(0, o12)

                    else:
                        nc_part.append(o1)
                        nc_part.append(o)

        # We do want a combined exponent if it would not be an Add, such as
        #  y    2y     3y
        # x  * x   -> x
        # We determine this if two exponents have the same term in as_coeff_mul
        #
        # Unfortunately, this isn't smart enough to consider combining into
        # exponents that might already be adds, so things like:
        #  z - y    y
        # x      * x  will be left alone.  This is because checking every possible
        # combination can slow things down.

        # gather exponents of common bases...
        # in c_powers
        new_c_powers = []
        common_b = {} # b:e

        for b, e in c_powers:
            co = e.as_coeff_mul()
            common_b.setdefault(b, {}).setdefault(co[1], []).append(co[0])
        for b, d in common_b.items():
            for di, li in d.items():
                d[di] = Add(*li)

        for b, e in common_b.items():
            for t, c in e.items():
                new_c_powers.append((b,c*Mul(*t)))
        c_powers = new_c_powers

        # and in num_exp
        new_num_exp = []
        common_b = {} # b:e

        for b, e in num_exp:
            co = e.as_coeff_mul()
            common_b.setdefault(b, {}).setdefault(co[1], []).append(co[0])
        for b, d in common_b.items():
            for di, li in d.items():
                d[di] = Add(*li)

        for b, e in common_b.items():
            for t, c in e.items():
                new_num_exp.append((b,c*Mul(*t)))
        num_exp = new_num_exp

        # --- PART 2 ---
        #
        # o process collected powers  (x**0 -> 1; x**1 -> x; otherwise Pow)
        # o combine collected powers  (2**x * 3**x -> 6**x)
        #   with numeric base

        # ................................
        # now we have:
        # - coeff:
        # - c_powers:    (b, e)
        # - num_exp:     (2, e)

        #  0             1
        # x  -> 1       x  -> x
        for b, e in c_powers:
            if e is S.Zero:
                continue

            if e is S.One:
                if b.is_Number:
                    coeff *= b
                else:
                    c_part.append(b)
            elif e.is_Integer and b.is_Number:
                coeff *= Pow(b, e)
            else:
                c_part.append(Pow(b, e))

        #  x    x     x
        # 2  * 3  -> 6
        inv_exp_dict = {}   # exp:Mul(num-bases)     x    x
                            # e.g.  x:6  for  ... * 2  * 3  * ...
        for b, e in num_exp:
            inv_exp_dict.setdefault(e, []).append(b)
        for e, b in inv_exp_dict.items():
            inv_exp_dict[e] = Mul(*b)

        reeval = False

        for e,b in inv_exp_dict.items():
            if e is S.Zero:
                continue

            if e is S.One:
                if b.is_Number:
                    coeff *= b
                else:
                    c_part.append(b)
            elif e.is_Integer and b.is_Number:
                coeff *= Pow(b, e)
            else:
                obj = b**e
                if obj.is_Mul:
                    # We may have split out a number that needs to go in coeff
                    # e.g., sqrt(6)*sqrt(2) == 2*sqrt(3).  See issue 415.
                    reeval = True
                if obj.is_Number:
                    coeff *= obj
                else:
                    c_part.append(obj)


        # oo, -oo
        if (coeff is S.Infinity) or (coeff is S.NegativeInfinity):
            new_c_part = []
            coeff_sign = 1
            for t in c_part:
                if t.is_positive:
                    continue
                if t.is_negative:
                    coeff_sign *= -1
                    continue
                new_c_part.append(t)
            c_part = new_c_part
            new_nc_part = []
            for t in nc_part:
                if t.is_positive:
                    continue
                if t.is_negative:
                    coeff_sign *= -1
                    continue
                new_nc_part.append(t)
            nc_part = new_nc_part
            coeff *= coeff_sign

        # 0, nan
        elif (coeff is S.Zero) or (coeff is S.NaN):
            # we know for sure the result will be the same as coeff (0 or nan)
            return [coeff], [], order_symbols

        elif coeff.is_Real:
            if coeff == Real(0):
                c_part, nc_part = [coeff], []
            elif coeff == Real(1):
                # change it to One, so it doesn't get inserted to slot0
                coeff = S.One


        # order commutative part canonically
        c_part.sort(Basic.compare)

        # current code expects coeff to be always in slot-0
        if coeff is not S.One:
            c_part.insert(0, coeff)


        # we are done
        if len(c_part)==2 and c_part[0].is_Number and c_part[1].is_Add:
            # 2*(1+a) -> 2 + 2 * a
            coeff = c_part[0]
            c_part = [Add(*[coeff*f for f in c_part[1].args])]

        if reeval:
            c_part, _, _ = Mul.flatten(c_part)
        return c_part, nc_part, order_symbols


    def _eval_power(b, e):
        if e.is_Number:
            if b.is_commutative:
                if e.is_Integer:
                    # (a*b)**2 -> a**2 * b**2
                    return Mul(*[s**e for s in b.args])

                if e.is_rational:
                    coeff, rest = b.as_coeff_mul()
                    unk=[]
                    nonneg=[]
                    neg=[]
                    for bi in rest:
                        if not bi.is_negative is None: #then we know the sign
                            if bi.is_negative:
                                neg.append(bi)
                            else:
                                nonneg.append(bi)
                        else:
                            unk.append(bi)
                    if len(unk) == len(rest) or len(neg) == len(rest) == 1:
                        # if all terms were unknown there is nothing to pull
                        # out except maybe the coeff OR if there was only a
                        # single negative term then it shouldn't be pulled out
                        # either.
                        if coeff < 0:
                            coeff *= -1
                        if coeff is S.One:
                            return None
                        return Mul(Pow(coeff, e), Pow(b/coeff, e))

                    # otherwise return the new expression expanding out the
                    # known terms; those that are not known can be expanded
                    # out with separate() but this will introduce a lot of
                    # "garbage" that is needed to keep one on the same branch
                    # as the unexpanded expression. The negatives are brought
                    # out with a negative sign added and a negative left behind
                    # in the unexpanded terms.
                    if neg:
                        neg = [-w for w in neg]
                        if len(neg) % 2 and not coeff.is_negative:
                            unk.append(S.NegativeOne)
                        if coeff.is_negative:
                            coeff = -coeff
                            unk.append(S.NegativeOne)
                    return Mul(*[Pow(s, e) for s in nonneg + neg + [coeff]])* \
                       Pow(Mul(*unk), e)


                coeff, rest = b.as_coeff_mul()
                if coeff is not S.One:
                    # (2*a)**3 -> 2**3 * a**3
                    return Mul(Pow(coeff, e), Mul(*[s**e for s in rest]))
            elif e.is_Integer:
                coeff, rest = b.as_coeff_mul()
                if coeff == S.One:
                    return # the test below for even exponent needs coeff != 1
                else:
                    return Mul(Pow(coeff, e), Pow(Mul(*rest), e))

        c, t = b.as_coeff_mul()
        if e.is_even and c.is_Number and c < 0:
            return Pow((Mul(-c, Mul(*t))), e)

        #if e.has(Wild):
        #    return Mul(*[t**e for t in b])

    def _eval_evalf(self, prec):
        return AssocOp._eval_evalf(self, prec).expand()

    @cacheit
    def as_two_terms(self):
        """Return head and tail of self.

        This is the most efficient way to get the head and tail of an
        expression.

        - if you want only the head, use self.args[0];
        - if you want to process the arguments of the tail then use
          self.as_coef_mul() which gives the head and a tuple containing
          the arguments of the tail when treated as a Mul.
        - if you want the coefficient when self is treated as an Add
          then use self.as_coeff_add()[0]

        >>> from sympy.abc import x, y
        >>> (3*x*y).as_two_terms()
        (3, x*y)
        """
        args = self.args

        if len(args) == 1:
            return S.One, self
        elif len(args) == 2:
            return args

        else:
            return args[0], self._new_rawargs(*args[1:])

    @cacheit
    def as_coeff_mul(self, *deps):
        if deps:
            l1 = []
            l2 = []
            for f in self.args:
                if f.has(*deps):
                    l2.append(f)
                else:
                    l1.append(f)
            return self._new_rawargs(*l1), tuple(l2)
        coeff, notrat = self.args[0].as_coeff_mul()
        if not coeff is S.One:
            return coeff, notrat + self.args[1:]
        return S.One, self.args

    @staticmethod
    def _expandsums(sums):
        """
        Helper function for _eval_expand_mul.

        sums must be a list of instances of Basic.
        """

        L = len(sums)
        if L == 1:
            return sums[0].args
        terms = []
        left = Mul._expandsums(sums[:L//2])
        right = Mul._expandsums(sums[L//2:])

        terms = [Mul(a, b) for a in left for b in right]
        added = Add(*terms)
        return Add.make_args(added) #it may have collapsed down to one term

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
        plain, sums, rewrite = [], [], False
        for factor in self.args:
            if deep:
                term = factor.expand(deep=deep, **hints)
                if term != factor:
                    factor = term
                    rewrite = True

            if factor.is_Add:
                sums.append(factor)
                rewrite = True
            else:
                if factor.is_commutative:
                    plain.append(factor)
                else:
                    Wrapper = Basic
                    sums.append(Wrapper(factor))

        if not rewrite:
            return self
        else:
            if sums:
                terms = Mul._expandsums(sums)
                plain = Mul(*plain)
                return Add(*[Mul(plain, term) for term in terms],
                           **self.assumptions0)
            else:
                return Mul(*plain, **self.assumptions0)

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

    def _eval_derivative(self, s):
        terms = list(self.args)
        factors = []
        for i in xrange(len(terms)):
            t = terms[i].diff(s)
            if t is S.Zero:
                continue
            factors.append(Mul(*(terms[:i]+[t]+terms[i+1:])))
        return Add(*factors)

    def _matches_simple(self, expr, repl_dict):
        # handle (w*3).matches('x*5') -> {w: x*5/3}
        coeff, terms = self.as_coeff_mul()
        if len(terms) == 1:
            return terms[0].matches(expr/coeff, repl_dict)
        return

    def matches(self, expr, repl_dict={}, evaluate=False):
        expr = sympify(expr)
        if self.is_commutative and expr.is_commutative:
            return AssocOp._matches_commutative(self, expr, repl_dict, evaluate)
        # todo for commutative parts, until then use the default matches method for non-commutative products
        return self._matches(expr, repl_dict, evaluate)

    def _matches(self, expr, repl_dict={}, evaluate=False):
        if evaluate:
            return self.subs(repl_dict).matches(expr, repl_dict)

        # weed out negative one prefixes
        sign = 1
        a, b = self.as_two_terms()
        if a is S.NegativeOne:
            if b.is_Mul:
                sign = -sign
            else:
                # the remainder, b, is not a Mul anymore
                return b.matches(-expr, repl_dict, evaluate)
        expr = sympify(expr)
        if expr.is_Mul and expr.args[0] is S.NegativeOne:
            expr = -expr; sign = -sign

        if not expr.is_Mul:
            # expr can only match if it matches b and a matches +/- 1
            if len(self.args) == 2:
                # quickly test for equality
                if b == expr:
                    return a.matches(Rational(sign), repl_dict, evaluate)
                # do more expensive match
                dd = b.matches(expr, repl_dict, evaluate)
                if dd == None:
                    return None
                dd = a.matches(Rational(sign), dd, evaluate)
                return dd
            return None

        d = repl_dict.copy()

        # weed out identical terms
        pp = list(self.args)
        ee = list(expr.args)
        for p in self.args:
            if p in expr.args:
                ee.remove(p)
                pp.remove(p)

        # only one symbol left in pattern -> match the remaining expression
        if len(pp) == 1 and isinstance(pp[0], C.Wild):
            if len(ee) == 1:
                d[pp[0]] = sign * ee[0]
            else:
                d[pp[0]] = sign * (type(expr)(*ee))
            return d

        if len(ee) != len(pp):
            return None

        i = 0
        for p, e in zip(pp, ee):
            if i == 0 and sign != 1:
                try:
                    e = sign * e
                except TypeError:
                    return None
            d = p.matches(e, d, evaluate=not i)
            i += 1
            if d is None:
                return None
        return d


    @staticmethod
    def _combine_inverse(lhs, rhs):
        """
        Returns lhs/rhs, but treats arguments like symbols, so things like
        oo/oo return 1, instead of a nan.
        """
        if lhs == rhs:
            return S.One
        if lhs.is_Mul and rhs.is_Mul:
            a = list(lhs.args)
            b = [1]
            for x in rhs.args:
                if x in a:
                    a.remove(x)
                else:
                    b.append(x)
            return Mul(*a)/Mul(*b)
        return lhs/rhs

    def as_powers_dict(self):
        d = {}
        for term in self.args:
            b, e = term.as_base_exp()
            if b not in d:
                d[b] = e
            else:
                d[b] += e
        return d

    def as_numer_denom(self):
        numers, denoms = [],[]
        for t in self.args:
            n,d = t.as_numer_denom()
            numers.append(n)
            denoms.append(d)
        return Mul(*numers), Mul(*denoms)

    def as_coeff_Mul(self):
        """Efficiently extract the coefficient of a product. """
        coeff, args = self.args[0], self.args[1:]

        if coeff.is_Number:
            if len(args) == 1:
                return coeff, args[0]
            else:
                return coeff, self._new_rawargs(*args)
        else:
            return S.One, self

    def _eval_is_polynomial(self, syms):
        for term in self.args:
            if not term._eval_is_polynomial(syms):
                return False
        return True

    _eval_is_bounded = lambda self: self._eval_template_is_attr('is_bounded')
    _eval_is_commutative = lambda self: self._eval_template_is_attr('is_commutative')
    _eval_is_integer = lambda self: self._eval_template_is_attr('is_integer')
    _eval_is_comparable = lambda self: self._eval_template_is_attr('is_comparable')


    # I*I -> R,  I*I*I -> -I

    def _eval_is_real(self):
        im_count = 0
        is_neither = False
        for t in self.args:
            if t.is_imaginary:
                im_count += 1
                continue
            t_real = t.is_real
            if t_real:
                continue
            elif t_real is False:
                if is_neither:
                    return None
                else:
                    is_neither = True
            else:
                return None
        if is_neither:
            return False

        return (im_count % 2 == 0)

    def _eval_is_imaginary(self):
        im_count = 0
        is_neither = False
        for t in self.args:
            if t.is_imaginary:
                im_count += 1
                continue
            t_real = t.is_real
            if t_real:
                continue
            elif t_real is False:
                if is_neither:
                    return None
                else:
                    is_neither = True
            else:
                return None
        if is_neither:
            return False

        return (im_count % 2 == 1)


    def _eval_is_irrational(self):
        for t in self.args:
            a = t.is_irrational
            if a:
                return True
            if a is None:
                return
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
        r = self._new_rawargs(*terms[1:])
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
            # all terms are either positive -- 2*Symbol('n', positive=T)
            #               or     unknown  -- 2*Symbol('x')
            if self.is_positive:
                return False
            else:
                return None
        c = terms[0]
        if len(terms)==1:
            return c.is_negative
        r = self._new_rawargs(*terms[1:])
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
        is_integer = self.is_integer

        if is_integer:
            r = True
            for t in self.args:
                if t.is_even:
                    return False
                if t.is_odd is None:
                    r = None
            return r

        # !integer -> !odd
        elif is_integer == False:
            return False


    def _eval_is_even(self):
        is_integer = self.is_integer

        if is_integer:
            return fuzzy_not(self._eval_is_odd())

        elif is_integer == False:
            return False

    def _eval_subs(self, old, new):

        from sympy import sign
        from sympy.simplify.simplify import powdenest

        if self == old:
            return new

        def fallback():
            """Return this value when partial subs has failed."""

            return self.__class__(*[s._eval_subs(old, new) for s in
                                  self.args])

        def breakup(eq):
            """break up powers assuming (not checking) that eq is a Mul:
                   b**(Rational*e) -> b**e, Rational
                commutatives come back as a dictionary {b**e: Rational}
                noncommutatives come back as a list [(b**e, Rational)]
            """

            (c, nc) = (dict(), list())
            for (i, a) in enumerate(Mul.make_args(eq) or [eq]): # remove or [eq] after 2114 accepted
                a = powdenest(a)
                (b, e) = a.as_base_exp()
                if not e is S.One:
                    (co, _) = e.as_coeff_mul()
                    b = Pow(b, e/co)
                    e = co
                if a.is_commutative:
                    if b in c: # handle I and -1 like things where b, e for I is -1, 1/2
                        c[b] += e
                    else:
                        c[b] = e
                else:
                    nc.append([b, e])
            return (c, nc)

        def rejoin(b, co):
            """
            Put rational back with exponent; in general this is not ok, but
            since we took it from the exponent for analysis, it's ok to put
            it back.
            """

            (b, e) = b.as_base_exp()
            return Pow(b, e*co)

        def ndiv(a, b):
            """if b divides a in an extractive way (like 1/4 divides 1/2
            but not vice versa, and 2/5 does not divide 1/3) then return
            the integer number of times it divides, else return 0.
            """

            if not b.q % a.q or not a.q % b.q:
                return int(a/b)
            return 0

        if not old.is_Mul:
            return fallback()

        # handle the leading coefficient and use it to decide if anything
        # should even be started; we always know where to find the Rational
        # so it's a quick test

        coeff = S.One
        co_self = self.args[0]
        co_old = old.args[0]
        if co_old.is_Rational and co_self.is_Rational:
            co_xmul = co_self.extract_multiplicatively(co_old)
        elif co_old.is_Rational:
            co_xmul = None
        else:
            co_xmul = True

        if not co_xmul:
            return fallback()

        (c, nc) = breakup(self)
        (old_c, old_nc) = breakup(old)

        # update the coefficients if we had an extraction

        if getattr(co_xmul, 'is_Rational', False):
            c.pop(co_self)
            c[co_xmul] = S.One
            old_c.pop(co_old)

        # do quick tests to see if we can't succeed

        ok = True
        if (
            # more non-commutative terms
            len(old_nc) > len(nc)):
            ok = False
        elif (
            # more commutative terms
            len(old_c) > len(c)):
            ok = False
        elif (
            # unmatched non-commutative bases
            set(_[0] for _ in  old_nc).difference(set(_[0] for _ in nc))):
            ok = False
        elif (
            # unmatched commutative terms
            set(old_c).difference(set(c))):
            ok = False
        elif (
            # differences in sign
            any(sign(c[b]) != sign(old_c[b]) for b in old_c)):
            ok = False
        if not ok:
            return fallback()

        if not old_c:
            cdid = None
        else:
            rat = []
            for (b, old_e) in old_c.items():
                c_e = c[b]
                rat.append(ndiv(c_e, old_e))
                if not rat[-1]:
                    return fallback()
            cdid = min(rat)

        if not old_nc:
            ncdid = None
            for i in range(len(nc)):
                nc[i] = rejoin(*nc[i])
        else:
            ncdid = 0  # number of nc replacements we did
            take = len(old_nc)  # how much to look at each time
            limit = cdid or S.Infinity # max number that we can take
            failed = []  # failed terms will need subs if other terms pass
            i = 0
            while limit and i + take <= len(nc):
                hit = False

                # the bases must be equivalent in succession, and
                # the powers must be extractively compatible on the
                # first and last factor but equal inbetween.

                rat = []
                for j in range(take):
                    if nc[i + j][0] != old_nc[j][0]:
                        break
                    elif j == 0:
                        rat.append(ndiv(nc[i + j][1], old_nc[j][1]))
                    elif j == take - 1:
                        rat.append(ndiv(nc[i + j][1], old_nc[j][1]))
                    elif nc[i + j][1] != old_nc[j][1]:
                        break
                    else:
                        rat.append(1)
                    j += 1
                else:
                    ndo = min(rat)
                    if ndo:
                        if take == 1:
                            if cdid:
                                ndo = min(cdid, ndo)
                            nc[i] = Pow(new, ndo)*rejoin(nc[i][0],
                                    nc[i][1] - ndo*old_nc[0][1])
                        else:
                            ndo = 1

                            # the left residual

                            l = rejoin(nc[i][0], nc[i][1] - ndo*
                                    old_nc[0][1])

                            # eliminate all middle terms

                            mid = new

                            # the right residual (which may be the same as the middle if take == 2)

                            ir = i + take - 1
                            r = (nc[ir][0], nc[ir][1] - ndo*
                                 old_nc[-1][1])
                            if r[1]:
                                if i + take < len(nc):
                                    nc[i:i + take] = [l*mid, r]
                                else:
                                    r = rejoin(*r)
                                    nc[i:i + take] = [l*mid*r]
                            else:

                                # there was nothing left on the right

                                nc[i:i + take] = [l*mid]

                        limit -= ndo
                        ncdid += ndo
                        hit = True
                if not hit:

                    # do the subs on this failing factor

                    failed.append(i)
                i += 1
            else:

                if not ncdid:
                    return fallback()

                # although we didn't fail, certain nc terms may have
                # failed so we rebuild them after attempting a partial
                # subs on them

                failed.extend(range(i, len(nc)))
                for i in failed:
                    nc[i] = rejoin(*nc[i]).subs(old, new)

        # rebuild the expression

        if cdid is None:
            do = ncdid
        elif ncdid is None:
            do = cdid
        else:
            do = min(ncdid, cdid)

        margs = []
        for b in c:
            if b in old_c:

                # calculate the new exponent

                e = c[b] - old_c[b]*do
                margs.append(rejoin(b, e))
            else:
                margs.append(rejoin(b.subs(old, new), c[b]))
        if cdid and not ncdid:

            # in case we are replacing commutative with non-commutative,
            # we want the new term to come at the front just like the
            # rest of this routine

            margs = [Pow(new, cdid)] + margs
        return Mul(*margs)*Mul(*nc)

    def _eval_nseries(self, x, n):
        from sympy import powsimp
        terms = [t.nseries(x, n=n) for t in self.args]
        return powsimp(Mul(*terms).expand(), combine='exp', deep=True)


    def _eval_as_leading_term(self, x):
        return Mul(*[t.as_leading_term(x) for t in self.args])

    def _eval_conjugate(self):
        return Mul(*[t.conjugate() for t in self.args])

    def _sage_(self):
        s = 1
        for x in self.args:
            s *= x._sage_()
        return s

from power import Pow
from numbers import Real, Integer, Rational
from function import FunctionClass
from sympify import sympify
from add import Add

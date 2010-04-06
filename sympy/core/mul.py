from basic import Basic, S
from operations import AssocOp
from cache import cacheit
from logic import fuzzy_not
from symbol import Symbol, Wild
from sympy.utilities.iterables import make_list

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
                o, order_symbols = o.as_expr_symbols(order_symbols)

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
                    # combined below
                    num_exp.append((b,e))
                    continue


                #         n          n          n
                # (-3 + y)   ->  (-1)  * (3 - y)
                if not Basic.keep_sign and b.is_Add and e.is_Number:
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
                    if b1==b2:
                        o12 = b1 ** (e1 + e2)

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
        # We determine this if two exponents have the same term in as_coeff_terms
        #
        # Unfortunately, this isn't smart enough to consider combining into
        # exponents that might already be adds, so things like:
        #  z - y    y
        # x      * x  will be left alone.  This is because checking every possible
        # combination can slow things down.
        new_c_powers = []
        common_b = {} # b:e

        # First gather exponents of common bases
        for b, e in c_powers:
            co = e.as_coeff_terms()
            if b in common_b:
                if  co[1] in common_b[b]:
                    common_b[b][co[1]] += co[0]
                else:
                    common_b[b][co[1]] = co[0]
            else:
                common_b[b] = {co[1]:co[0]}

        for b,e, in common_b.items():
            for t, c in e.items():
                new_c_powers.append((b,c*Mul(*t)))
        c_powers = new_c_powers

        # And the same for numeric bases
        new_num_exp = []
        common_b = {} # b:e
        for b, e in num_exp:
            co = e.as_coeff_terms()
            if b in common_b:
                if  co[1] in common_b[b]:
                    common_b[b][co[1]] += co[0]
                else:
                    common_b[b][co[1]] = co[0]
            else:
                common_b[b] = {co[1]:co[0]}

        for b,e, in common_b.items():
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
                coeff *= b ** e
            else:
                c_part.append(Pow(b, e))

        #  x    x     x
        # 2  * 3  -> 6
        inv_exp_dict = {}   # exp:Mul(num-bases)     x    x
                            # e.g.  x:6  for  ... * 2  * 3  * ...
        for b,e in num_exp:
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


        # oo, -oo
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

        return c_part, nc_part, order_symbols


    def _eval_power(b, e):
        if e.is_Number:
            if b.is_commutative:
                if e.is_Integer:
                    # (a*b)**2 -> a**2 * b**2
                    return Mul(*[s**e for s in b.args])

                if e.is_rational:
                    coeff, rest = b.as_coeff_terms()
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
                            coeff = -coeff
                        if coeff == S.One:
                            return None
                        b = b / coeff
                        return coeff ** e * b ** e

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
                    return Mul(*[s**e for s in nonneg + neg + [coeff]])* \
                       Mul(*(unk)) ** e


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

        c, t = b.as_coeff_terms()
        if e.is_even and c.is_Number and c < 0:
            return (-c * Mul(*t)) ** e

        #if e.atoms(Wild):
        #    return Mul(*[t**e for t in b])

    def _eval_evalf(self, prec):
        return AssocOp._eval_evalf(self, prec).expand()

    @cacheit
    def as_two_terms(self):
        args = self.args

        if len(args) == 1:
            return S.One, self
        elif len(args) == 2:
            return args

        else:
            return args[0], self._new_rawargs(*args[1:])

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
        return make_list(added, Add) #it may have collapsed down to one term

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
        coeff, terms = self.as_coeff_terms()
        if len(terms)==1:
            return terms[0].matches(expr / coeff, repl_dict)
        return

    def matches(self, expr, repl_dict={}, evaluate=False):
        expr = sympify(expr)
        if self.is_commutative and expr.is_commutative:
            return AssocOp._matches_commutative(self, expr, repl_dict, evaluate)
        # todo for commutative parts, until then use the default matches method for non-commutative products
        return self._matches(expr, repl_dict, evaluate)

    def _matches(self, expr, repl_dict={}, evaluate=False):
        # weed out negative one prefixes
        sign = 1
        if self.args[0] == -1:
            self = -self; sign = -sign
        if expr.is_Mul and expr.args[0] == -1:
            expr = -expr; sign = -sign

        if evaluate:
            return self.subs(repl_dict).matches(expr, repl_dict)
        expr = sympify(expr)
        if not isinstance(expr, self.__class__):
            # if we can omit the first factor, we can match it to sign * one
            if Mul(*self.args[1:]) == expr:
                return self.args[0].matches(Rational(sign), repl_dict, evaluate)
            # two-factor product: if the 2nd factor matches, the first part must be sign * one
            if len(self.args[:]) == 2:
                dd = self.args[1].matches(expr, repl_dict, evaluate)
                if dd == None:
                    return None
                dd = self.args[0].matches(Rational(sign), dd, evaluate)
                return dd
            return None

        if len(self.args[:])==0:
            if self == expr:
                return repl_dict
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
        if len(pp) == 1 and isinstance(pp[0], Wild):
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
            return Add(*[t.count_ops(symbolic) for t in self.args]) + \
                Symbol('MUL') * (len(self.args) - 1)
        return Add(*[t.count_ops(symbolic) for t in self.args]) + \
            (len(self.args) - 1)

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
            # all terms are either positive -- 2*Symbol('n', positive=T)
            #               or     unknown  -- 2*Symbol('x')
            if self.is_positive:
                return False
            else:
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
        # base cases
        # simplest
        if self == old:
            return new
        # pass it off to its own class
        if isinstance(old, FunctionClass):
            return self.__class__(*[s._eval_subs(old, new) for s in self.args ])

        # break up self and old into terms
        coeff_self,terms_self = self.as_coeff_terms()
        coeff_old,terms_old = old.as_coeff_terms()

        # NEW - implementation of strict substitution
        # if the coefficients are not the same, do not substitute.
        # the only exception is if old has a coefficient of 1, then always to the sub.
        if coeff_self != coeff_old and coeff_old != 1:
            return self.__class__(*[s._eval_subs(old, new) for s in self.args])

        # break up powers, i.e., x**2 -> x*x
        def breakup(terms):
            temp = []
            for t in terms:
                if isinstance(t,Pow) and isinstance(t.exp, Integer):
                    if t.exp.is_positive:
                        temp.extend([t.base]*int(t.exp))
                    elif t.exp.is_negative:
                        temp.extend([1/t.base]*int(abs(t.exp)))
                else:
                    temp.append(t)
            return temp
        terms_old = breakup(terms_old)
        terms_self = breakup(terms_self)

        # break up old and self terms into commutative and noncommutative lists
        comm_old = []; noncomm_old = []
        comm_self = []; noncomm_self = []
        for o in terms_old:
            if o.is_commutative:
                comm_old.append(o)
            else:
                noncomm_old.append(o)
        for s in terms_self:
            if s.is_commutative:
                comm_self.append(s)
            else:
                noncomm_self.append(s)
        comm_old_len, noncomm_old_len = len(comm_old), len(noncomm_old)
        comm_self_len, noncomm_self_len = len(comm_self), len(noncomm_self)

        # if the noncommutative part of the 'to-be-replaced' expression is
        # smaller than the noncommutative part of the whole expression, scan
        # to see if the whole thing is there
        if noncomm_old_len <= noncomm_self_len and noncomm_old_len > 0:
            for i in range(noncomm_self_len):
                if noncomm_self[i] == noncomm_old[0]:
                    for j in range(noncomm_old_len):
                        # make sure each noncommutative term matches in order
                        if (i+j) < noncomm_self_len and \
                           noncomm_self[i+j] == noncomm_old[j]:
                            # we only care once we've reached the end of old's
                            # noncommutative part.
                            if j == noncomm_old_len-1:
                                # get rid of noncommutative terms and
                                # substitute new expression into total
                                # expression
                                noncomms_final = noncomm_self[:i] + \
                                                 noncomm_self[i+j+1:]
                                noncomms_final.insert(i,new)

                                myFlag = True
                                comms_final = comm_self[:]
                                # check commutative terms
                                for ele in comm_old:
                                    # flag to make sure all the commutative
                                    # terms in old are in self
                                    if ele not in comm_self:
                                        myFlag = False
                                    # collect commutative terms
                                    else:
                                        comms_final.remove(ele)

                                # continue only if all commutative terms in
                                # old are present
                                if myFlag == True:
                                    expr = comms_final+noncomms_final
                                    return Mul(coeff_self/coeff_old,
                                               Mul(*expr)._eval_subs(old,new))
                                               #*[e._eval_subs(old,new) for e in expr])

            return self.__class__(*[s._eval_subs(old, new) for s in self.args])

        # but what if the noncommutative lists subexpression and the whole
        # expression are both empty
        elif noncomm_old_len == noncomm_self_len == 0:
            # just check commutative parts then.
            if comm_old_len > 0 and comm_old_len<=comm_self_len:
                if comm_self == comm_old:
                    return Mul(coeff_self/coeff_old*new)
                myFlag = True
                comms_final = comm_self[:]
                # check commutative terms
                for ele in comm_old:
                    # flag to make sure all the commutative terms in old are
                    # in self
                    if ele not in comm_self:
                        myFlag = False
                    # collect commutative terms
                    else:
                        # needed if old has an element to an integer power
                        if ele in comms_final:
                            comms_final.remove(ele)
                        else:
                            myFlag = False

                # continue only if all commutative terms in old are present
                if myFlag == True:
                    return Mul(coeff_self/coeff_old, new,
                               Mul(*comms_final)._eval_subs(old,new))#*[c._eval_subs(old,new) for c in comms_final])
                else:
                    return self.__class__(*[s._eval_subs(old, new) for
                                            s in self.args])

        # else the subexpression isn't in the totaly expression
        return self.__class__(*[s._eval_subs(old, new) for s in self.args])

    def _eval_nseries(self, x, x0, n):
        from sympy import powsimp
        terms = [t.nseries(x, x0, n) for t in self.args]
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

    def as_Mul(self):
        """Returns `self` as it was `Mul` instance. """
        return list(self.args)

from power import Pow
from numbers import Real, Integer, Rational
from function import FunctionClass
from sympify import sympify
from add import Add


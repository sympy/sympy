
from basic import Basic, S, sympify
from operations import AssocOp
from cache import cacheit

from logic import fuzzy_not

from symbol import Symbol, Wild
# from function import FunctionClass, WildFunction /cyclic/
# from numbers import Number, Integer, Real /cyclic/
# from add   import Add /cyclic/
# from power import Pow /cyclic/

import sympy.mpmath as mpmath

# internal marker to indicate:
#   "there are still non-commutative objects -- don't forget to processe them"
class NC_Marker:
    is_Order    = False
    is_Mul      = False
    is_Number   = False

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

        cnum_powers = {}    # base:num-exp                  2
                            # e.g. (x+y):2  for  ... * (x+y)  * ...

        c_powers = []       # (base,exp)      z
                            # e.g. (x,z) for x

        num_exp = []     # (num-base, exp)           y
                            # e.g.  (3, y)  for  ... * 3  * ...

        order_symbols = None



        # --- PART 1 ---
        #
        # "collect powers and coeff":
        #
        # o coeff
        # o c_powers
        # o exp_dict
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

                # let's collect factors with the same base if the exponent is
                # a number, so e.g.
                #  2    3     5
                # x  * x  -> x
                if e.is_Number:
                    if b in cnum_powers:
                        cnum_powers[b] += e
                    else:
                        cnum_powers[b] = e
                else:
                    c_powers.append((b,e))


            # NON-COMMUTATIVE
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
        new_c_powers = []
        common_b = {} # b:e
        for b, e in c_powers:
            if b in common_b:
                if (e+common_b[b]).is_Add:
                    new_c_powers.append((b,e))
                else:
                    common_b[b] += e # Only add exponents if the sum is not an Add
            else:
                common_b[b] = e

        for b,e, in common_b.items():
            new_c_powers.append((b,e))
        c_powers = new_c_powers

        # And the same for numeric bases
        new_num_exp = []
        common_b = {} # b:e
        for b, e in num_exp:
            if b in common_b:
                if (e+common_b[b]).is_Add:
                    new_num_exp.append((b,e))
                else:
                    common_b[b] += e # Only add exponents if the sum is not an Add
            else:
                common_b[b] = e

        for b,e, in common_b.items():
            new_num_exp.append((b,e))
        num_exp = new_num_exp
        #

        # --- PART 2 ---
        #
        # o process collected powers  (x**0 -> 1; x**1 -> x; otherwise Pow)
        # o combine collected powers  (2**x * 3**x -> 6**x)
        #   with numeric base

        # ................................
        # now we have:
        # - coeff:
        # - cnum_powers:  b:3
        # - c_powers:    (b, e)
        # - num_exp:     (2, e)

        #  0             1
        # x  -> 1       x  -> x
        for b, e in cnum_powers.items() + c_powers:
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
                            # e.g.  x -> 6  for  ... * 2  * 3  * ...
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

                return Add(*(Mul(plain, term) for term in terms), **self.assumptions0)
            else:
                return Mul(*plain, **self.assumptions0)

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
        """
        Returns lhs/rhs, but treats arguments like symbols, so things like
        oo/oo return 1, instead of a nan.
        """
        if lhs == rhs:
            return S.One
        if lhs.is_Mul and rhs.is_Mul:
            a = list(lhs.args[:])
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
        re_not   = False

        for t in self.args:
            if t.is_imaginary:
                im_count += 1
                continue

            t_real = t.is_real
            if t_real:
                continue

            elif fuzzy_not(t_real):
                re_not = True

            else:
                return None

        if re_not:
            return False

        return (im_count % 2 == 0)


    def _eval_is_imaginary(self):
        im_count = 0

        for t in self.args:
            if t.is_imaginary:
                im_count += 1

            elif t.is_real:
                continue

            # real=F|U
            else:
                return None

        return (im_count % 2 == 1)



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
        if self == old:
            return new
        if isinstance(old, FunctionClass):
            return self.__class__(*[s._eval_subs(old, new) for s in self.args ])
        coeff_self,terms_self = self.as_coeff_terms()
        coeff_old,terms_old = old.as_coeff_terms()
        if terms_self == terms_old: # (2*a).subs(3*a,y) -> 2/3*y
            return new * coeff_self/coeff_old
        l1, l2 = len(terms_self), len(terms_old)
        if l2 == 0:
            # if old is just a number, go through the self.args one by one
            return Mul(*[x._eval_subs(old, new) for x in self.args])
        elif l2 < l1:
            # old is some something more complex, like:
            # (a*b*c*d).subs(b*c,x) -> a*x*d
            # then we need to search where in self.args the "old" is, and then
            # correctly substitute both terms and coefficients.
            self_set = set(terms_self)
            old_set = set(terms_old)
            if old_set < self_set:
                ret_set = self_set - old_set
                return Mul(new, coeff_self/coeff_old, *[s._eval_subs(old, new) for s in ret_set])
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


# /cyclic/
import basic as _
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

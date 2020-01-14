"""Provides the basic operations between the operators."""

from .operator import OpExpr, Op
from collections import defaultdict
from sympy.core.add import Add, _addsort
from sympy.core.evalf import PrecisionExhausted
from sympy.core.logic import fuzzy_not
from sympy.core.mul import Mul, NC_Marker, _mulsort, _keep_coeff, _unevaluated_Mul
from sympy.core.numbers import Rational
from sympy.core.parameters import global_parameters
from sympy.core.power import Pow
from sympy.core.singleton import S
from sympy.core.sympify import _sympify

class OpAdd(OpExpr, Add):
    """Added operator.

    Examples
    ========

    >>> from sympy import sin, cos
    >>> from sympy.abc import x

    >>> 1+cos+sin
    1 + cos + sin
    >>> (1+cos+sin)(x)
    sin(x) + cos(x) + 1
    """
    identity = S.Zero

    @classmethod
    def flatten(cls, seq):
        from sympy.calculus.util import AccumBounds
        from sympy.matrices.expressions import MatrixExpr
        from sympy.tensor.tensor import TensExpr
        rv = None
        if len(seq) == 2:
            a, b = seq
            if b.is_Rational:
                a, b = b, a
            if a.is_Rational:
                if b.is_Mul:
                    rv = [a, b], [], None
            if rv:
                if all(s.is_commutative for s in rv[0]):
                    return rv
                return [], rv[0], None

        terms = {}      # term -> coeff
                        # e.g. x**2 -> 5   for ... + 5*x**2 + ...

        coeff = S.Zero  # coefficient (Number or zoo) to always be in slot 0
                        # e.g. 3 + ...
        order_factors = []

        extra = []

        for o in seq:

            # O(x)
            if o.is_Order:
                for o1 in order_factors:
                    if o1.contains(o):
                        o = None
                        break
                if o is None:
                    continue
                order_factors = [o] + [
                    o1 for o1 in order_factors if not o.contains(o1)]
                continue

            # 3 or NaN
            elif o.is_Number:
                if (o is S.NaN or coeff is S.ComplexInfinity and
                        o.is_finite is False) and not extra:
                    # we know for sure the result will be nan
                    return [S.NaN], [], None
                if coeff.is_Number:
                    coeff += o
                    if coeff is S.NaN and not extra:
                        # we know for sure the result will be nan
                        return [S.NaN], [], None
                continue

            elif isinstance(o, AccumBounds):
                coeff = o.__add__(coeff)
                continue

            elif isinstance(o, MatrixExpr):
                # can't add 0 to Matrix so make sure coeff is not 0
                extra.append(o)
                continue

            elif isinstance(o, TensExpr):
                coeff = o.__add__(coeff) if coeff else o
                continue

            elif o is S.ComplexInfinity:
                if coeff.is_finite is False and not extra:
                    # we know for sure the result will be nan
                    return [S.NaN], [], None
                coeff = S.ComplexInfinity
                continue

            # Add([...])
            elif o.is_Add:
                # NB: here we assume Add is always commutative
                seq.extend(o.args)  # TODO zerocopy?
                continue

            # Mul([...])
            elif o.is_Mul:
                c, s = o.as_coeff_Mul()

            # check for unevaluated Pow, e.g. 2**3 or 2**(-1/2)
            elif o.is_Pow:
                b, e = o.as_base_exp()
                if b.is_Number and (e.is_Integer or
                                   (e.is_Rational and e.is_negative)):
                    seq.append(b**e)
                    continue
                c, s = S.One, o

            else:
                # everything else
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
                if terms[s] is S.NaN and not extra:
                    # we know for sure the result will be nan
                    return [S.NaN], [], None
            else:
                terms[s] = c

        # now let's construct new args:
        # [2*x**2, x**3, 7*x**4, pi, ...]
        newseq = []
        noncommutative = False
        for s, c in terms.items():
            # 0*s
            if c.is_zero:
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
                elif s.is_Add:
                    # we just re-create the unevaluated Mul
                    if not isinstance(s, OpExpr):
                        newseq.append(Mul(c, s, evaluate=False))
                    else:
                        newseq.append(OpMul(c, s, evaluate=False))
                else:
                    # alternatively we have to call all Mul's machinery (slow)
                    if not isinstance(s, OpExpr):
                        newseq.append(Mul(c, s))
                    else:
                        newseq.append(OpMul(c, s))

            noncommutative = noncommutative or not s.is_commutative

        # oo, -oo
        if coeff is S.Infinity:
            newseq = [f for f in newseq if not (f.is_extended_nonnegative or f.is_real)]

        elif coeff is S.NegativeInfinity:
            newseq = [f for f in newseq if not (f.is_extended_nonpositive or f.is_real)]

        if coeff is S.ComplexInfinity:
            # zoo might be
            #   infinite_real + finite_im
            #   finite_real + infinite_im
            #   infinite_real + infinite_im
            # addition of a finite real or imaginary number won't be able to
            # change the zoo nature; adding an infinite qualtity would result
            # in a NaN condition if it had sign opposite of the infinite
            # portion of zoo, e.g., infinite_real - infinite_real.
            newseq = [c for c in newseq if not (c.is_finite and
                                                c.is_extended_real is not None)]

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
        _addsort(newseq)

        # current code expects coeff to be first
        if coeff is not S.Zero:
            newseq.insert(0, coeff)

        if extra:
            newseq += extra
            noncommutative = True

        # we are done
        if noncommutative:
            return [], newseq, None
        else:
            return newseq, [], None

    def __new__(cls, *args, **kwargs):
        new_args = cls._process_operator(*args)
        obj = super(cls, OpAdd).__new__(cls, *new_args, **kwargs)
        if not isinstance(obj, OpExpr):
            return Op(obj)
        return obj

    def _eval_operation(self, *args, **kwargs):
        args = [a(*args, evaluate=True, **kwargs) if isinstance(a, OpExpr) else a for a in self.args]
        return Add(*args)

class OpMul(OpExpr, Mul):
    """Multiplied operator.

    Examples
    ========

    >>> from sympy import sin, cos
    >>> from sympy.abc import x

    >>> 2*cos
    2*cos
    >>> 1*sin
    sin
    >>> cos*sin
    cos*sin

    >>> (cos*sin*2)(x)
    2*sin(x)*cos(x)
    """
    identity = S.One

    def __new__(cls, *args, **kwargs):
        new_args = cls._process_operator(*args)
        obj = super(cls, OpMul).__new__(cls, *new_args, **kwargs)
        if not isinstance(obj, OpExpr):
            return Op(obj)
        return obj

    @classmethod
    def flatten(cls, seq):
        from sympy.calculus.util import AccumBounds
        from sympy.matrices.expressions import MatrixExpr
        rv = None
        if len(seq) == 2:
            a, b = seq
            if b.is_Rational:
                a, b = b, a
                seq = [a, b]
            assert not a is S.One
            if not a.is_zero and a.is_Rational:
                r, b = b.as_coeff_Mul()
                if b.is_Add:
                    if r is not S.One:  # 2-arg hack
                        # leave the Mul as a Mul
                        rv = [cls(a*r, b, evaluate=False)], [], None
                    elif global_parameters.distribute and b.is_commutative:
                        r, b = b.as_coeff_Add()
                        bargs = [_keep_coeff(a, bi) for bi in Add.make_args(b)]
                        _addsort(bargs)
                        ar = a*r
                        if ar:
                            bargs.insert(0, ar)
                        bargs = [Add._from_args(bargs)]
                        rv = bargs, [], None
            if rv:
                return rv

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

        neg1e = S.Zero      # exponent on -1 extracted from Number-based Pow and I

        pnum_rat = {}       # (num-base, Rat-exp)          1/2
                            # e.g.  (3, 1/2)  for  ... * 3     * ...

        order_symbols = None

        # --- PART 1 ---
        #
        # "collect powers and coeff":
        #
        # o coeff
        # o c_powers
        # o num_exp
        # o neg1e
        # o pnum_rat
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
                if o is S.NaN or coeff is S.ComplexInfinity and o.is_zero:
                    # we know for sure the result will be nan
                    return [S.NaN], [], None
                elif coeff.is_Number or isinstance(coeff, AccumBounds):  # it could be zoo
                    coeff *= o
                    if coeff is S.NaN:
                        # we know for sure the result will be nan
                        return [S.NaN], [], None
                continue

            elif isinstance(o, AccumBounds):
                coeff = o.__mul__(coeff)
                continue

            elif o is S.ComplexInfinity:
                if not coeff:
                    # 0 * zoo = NaN
                    return [S.NaN], [], None
                if coeff is S.ComplexInfinity:
                    # zoo * zoo = zoo
                    return [S.ComplexInfinity], [], None
                coeff = S.ComplexInfinity
                continue

            elif o is S.ImaginaryUnit:
                neg1e += S.Half
                continue

            elif o.is_commutative:
                #      e
                # o = b
                b, e = o.as_base_exp()

                #  y
                # 3
                if o.is_Pow:
                    if b.is_Number:

                        # get all the factors with numeric base so they can be
                        # combined below, but don't combine negatives unless
                        # the exponent is an integer
                        if e.is_Rational:
                            if e.is_Integer:
                                coeff *= Pow(b, e)  # it is an unevaluated power
                                continue
                            elif e.is_negative:    # also a sign of an unevaluated power
                                seq.append(Pow(b, e))
                                continue
                            elif b.is_negative:
                                neg1e += e
                                b = -b
                            if b is not S.One:
                                pnum_rat.setdefault(b, []).append(e)
                            continue
                        elif b.is_positive or e.is_integer:
                            num_exp.append((b, e))
                            continue

                c_powers.append((b, e))

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
                    b1, e1 = o1.as_base_exp()
                    b2, e2 = o.as_base_exp()
                    new_exp = e1 + e2
                    # Only allow powers to combine if the new exponent is
                    # not an Add. This allow things like a**2*b**3 == a**5
                    # if a.is_commutative == False, but prohibits
                    # a**x*a**y and x**a*x**b from combining (x,y commute).
                    if b1 == b2 and (not new_exp.is_Add):
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
        # We determine if two exponents have the same term by using
        # as_coeff_Mul.
        #
        # Unfortunately, this isn't smart enough to consider combining into
        # exponents that might already be adds, so things like:
        #  z - y    y
        # x      * x  will be left alone.  This is because checking every possible
        # combination can slow things down.

        # gather exponents of common bases...
        def _gather(c_powers):
            common_b = {}  # b:e
            for b, e in c_powers:
                co = e.as_coeff_Mul()
                common_b.setdefault(b, {}).setdefault(
                    co[1], []).append(co[0])
            for b, d in common_b.items():
                for di, li in d.items():
                    d[di] = Add(*li)
            new_c_powers = []
            for b, e in common_b.items():
                new_c_powers.extend([(b, c*t) for t, c in e.items()])
            return new_c_powers

        # in c_powers
        c_powers = _gather(c_powers)

        # and in num_exp
        num_exp = _gather(num_exp)

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
        # - pnum_rat:    {(1/3, [1/3, 2/3, 1/4])}

        #  0             1
        # x  -> 1       x  -> x

        # this should only need to run twice; if it fails because
        # it needs to be run more times, perhaps this should be
        # changed to a "while True" loop -- the only reason it
        # isn't such now is to allow a less-than-perfect result to
        # be obtained rather than raising an error or entering an
        # infinite loop
        for i in range(2):
            new_c_powers = []
            changed = False
            for b, e in c_powers:
                if e.is_zero:
                    # canceling out infinities yields NaN
                    if (b.is_Add or b.is_Mul) and any(infty in b.args
                        for infty in (S.ComplexInfinity, S.Infinity,
                                      S.NegativeInfinity)):
                        return [S.NaN], [], None
                    continue
                if e is S.One:
                    if b.is_Number:
                        coeff *= b
                        continue
                    p = b
                if e is not S.One:
                    p = Pow(b, e)
                    # check to make sure that the base doesn't change
                    # after exponentiation; to allow for unevaluated
                    # Pow, we only do so if b is not already a Pow
                    if p.is_Pow and not b.is_Pow:
                        bi = b
                        b, e = p.as_base_exp()
                        if b != bi:
                            changed = True
                c_part.append(p)
                new_c_powers.append((b, e))
            # there might have been a change, but unless the base
            # matches some other base, there is nothing to do
            if changed and len(set(
                    b for b, e in new_c_powers)) != len(new_c_powers):
                # start over again
                c_part = []
                c_powers = _gather(new_c_powers)
            else:
                break

        #  x    x     x
        # 2  * 3  -> 6
        inv_exp_dict = {}   # exp:Mul(num-bases)     x    x
                            # e.g.  x:6  for  ... * 2  * 3  * ...
        for b, e in num_exp:
            inv_exp_dict.setdefault(e, []).append(b)
        for e, b in inv_exp_dict.items():
            inv_exp_dict[e] = cls(*b)
        c_part.extend([Pow(b, e) for e, b in inv_exp_dict.items() if e])

        # b, e -> e' = sum(e), b
        # {(1/5, [1/3]), (1/2, [1/12, 1/4]} -> {(1/3, [1/5, 1/2])}
        comb_e = {}
        for b, e in pnum_rat.items():
            comb_e.setdefault(Add(*e), []).append(b)
        del pnum_rat
        # process them, reducing exponents to values less than 1
        # and updating coeff if necessary else adding them to
        # num_rat for further processing
        num_rat = []
        for e, b in comb_e.items():
            b = cls(*b)
            if e.q == 1:
                coeff *= Pow(b, e)
                continue
            if e.p > e.q:
                e_i, ep = divmod(e.p, e.q)
                coeff *= Pow(b, e_i)
                e = Rational(ep, e.q)
            num_rat.append((b, e))
        del comb_e

        # extract gcd of bases in num_rat
        # 2**(1/3)*6**(1/4) -> 2**(1/3+1/4)*3**(1/4)
        pnew = defaultdict(list)
        i = 0  # steps through num_rat which may grow
        while i < len(num_rat):
            bi, ei = num_rat[i]
            grow = []
            for j in range(i + 1, len(num_rat)):
                bj, ej = num_rat[j]
                g = bi.gcd(bj)
                if g is not S.One:
                    # 4**r1*6**r2 -> 2**(r1+r2)  *  2**r1 *  3**r2
                    # this might have a gcd with something else
                    e = ei + ej
                    if e.q == 1:
                        coeff *= Pow(g, e)
                    else:
                        if e.p > e.q:
                            e_i, ep = divmod(e.p, e.q)  # change e in place
                            coeff *= Pow(g, e_i)
                            e = Rational(ep, e.q)
                        grow.append((g, e))
                    # update the jth item
                    num_rat[j] = (bj/g, ej)
                    # update bi that we are checking with
                    bi = bi/g
                    if bi is S.One:
                        break
            if bi is not S.One:
                obj = Pow(bi, ei)
                if obj.is_Number:
                    coeff *= obj
                else:
                    # changes like sqrt(12) -> 2*sqrt(3)
                    for obj in Mul.make_args(obj):
                        if obj.is_Number:
                            coeff *= obj
                        else:
                            assert obj.is_Pow
                            bi, ei = obj.args
                            pnew[ei].append(bi)

            num_rat.extend(grow)
            i += 1

        # combine bases of the new powers
        for e, b in pnew.items():
            pnew[e] = cls(*b)

        # handle -1 and I
        if neg1e:
            # treat I as (-1)**(1/2) and compute -1's total exponent
            p, q =  neg1e.as_numer_denom()
            # if the integer part is odd, extract -1
            n, p = divmod(p, q)
            if n % 2:
                coeff = -coeff
            # if it's a multiple of 1/2 extract I
            if q == 2:
                c_part.append(S.ImaginaryUnit)
            elif p:
                # see if there is any positive base this power of
                # -1 can join
                neg1e = Rational(p, q)
                for e, b in pnew.items():
                    if e == neg1e and b.is_positive:
                        pnew[e] = -b
                        break
                else:
                    # keep it separate; we've already evaluated it as
                    # much as possible so evaluate=False
                    c_part.append(Pow(S.NegativeOne, neg1e, evaluate=False))

        # add all the pnew powers
        c_part.extend([Pow(b, e) for e, b in pnew.items()])

        # oo, -oo
        if (coeff is S.Infinity) or (coeff is S.NegativeInfinity):
            def _handle_for_oo(c_part, coeff_sign):
                new_c_part = []
                for t in c_part:
                    if t.is_extended_positive:
                        continue
                    if t.is_extended_negative:
                        coeff_sign *= -1
                        continue
                    new_c_part.append(t)
                return new_c_part, coeff_sign
            c_part, coeff_sign = _handle_for_oo(c_part, 1)
            nc_part, coeff_sign = _handle_for_oo(nc_part, coeff_sign)
            coeff *= coeff_sign

        # zoo
        if coeff is S.ComplexInfinity:
            # zoo might be
            #   infinite_real + bounded_im
            #   bounded_real + infinite_im
            #   infinite_real + infinite_im
            # and non-zero real or imaginary will not change that status.
            c_part = [c for c in c_part if not (fuzzy_not(c.is_zero) and
                                                c.is_extended_real is not None)]
            nc_part = [c for c in nc_part if not (fuzzy_not(c.is_zero) and
                                                  c.is_extended_real is not None)]

        # 0
        elif coeff.is_zero:
            # we know for sure the result will be 0 except the multiplicand
            # is infinity or a matrix
            if any(isinstance(c, MatrixExpr) for c in nc_part):
                return [coeff], nc_part, order_symbols
            if any(c.is_finite == False for c in c_part):
                return [S.NaN], [], order_symbols
            return [coeff], [], order_symbols

        # check for straggling Numbers that were produced
        _new = []
        for i in c_part:
            if i.is_Number:
                coeff *= i
            else:
                _new.append(i)
        c_part = _new

        # order commutative part canonically
        _mulsort(c_part)

        # current code expects coeff to be always in slot-0
        if coeff is not S.One:
            c_part.insert(0, coeff)

        # we are done
        if (global_parameters.distribute and not nc_part and len(c_part) == 2 and
                c_part[0].is_Number and c_part[0].is_finite and c_part[1].is_Add):
            # 2*(1+a) -> 2 + 2 * a
            coeff = c_part[0]
            c_part = [Add(*[coeff*f for f in c_part[1].args])]

        return c_part, nc_part, order_symbols

    def _eval_operation(self, *args, **kwargs):
        args = [a(*args, evaluate=True, **kwargs) if isinstance(a, OpExpr) else a for a in self.args]
        return Mul(*args)

    def _eval_power(b, e):

        # don't break up NC terms: (A*B)**3 != A**3*B**3, it is A*B*A*B*A*B
        cargs, nc = b.args_cnc(split_1=False)

        if e.is_Integer:
            if not isinstance(b, OpExpr):
                return Mul(*[Pow(b, e, evaluate=False) for b in cargs]) * \
                    Pow(Mul._from_args(nc), e, evaluate=False)
            else:
                return OpMul(*[OpPow(b, e, evaluate=False) for b in cargs]) * \
                    OpPow(OpMul._from_args(nc), e, evaluate=False)
        if e.is_Rational and e.q == 2:
            from sympy.core.power import integer_nthroot
            from sympy.functions.elementary.complexes import sign
            if b.is_imaginary:
                a = b.as_real_imag()[1]
                if a.is_Rational:
                    n, d = abs(a/2).as_numer_denom()
                    n, t = integer_nthroot(n, 2)
                    if t:
                        d, t = integer_nthroot(d, 2)
                        if t:
                            r = _sympify(n)/d
                            return _unevaluated_Mul(r**e.p, (1 + sign(a)*S.ImaginaryUnit)**e.p)

        p = Pow(b, e, evaluate=False)

        if e.is_Rational or e.is_Float:
            return p._eval_expand_power_base()

        return p

class OpPow(OpExpr, Pow):
    """Powered operator.

    Examples
    ========

    >>> from sympy import sin, cos
    >>> from sympy.abc import x

    >>> 2**cos
    2**cos
    >>> cos**4
    cos**4
    >>> cos**sin
    cos**sin

    >>> (cos**sin)(x)
    cos(x)**sin(x)
    >>> (sin/cos)(x)
    sin(x)/cos(x)
    >>> (sin**x)(x)
    sin(x)**x
    """
    def __new__(cls, b, e, **kwargs):
        new_b, new_e = cls._process_operator(b, e)
        obj = super(cls, OpPow).__new__(cls, new_b, new_e, **kwargs)
        if not isinstance(obj, OpExpr):
            return Op(obj)
        return obj

    def _eval_operation(self, *args, **kwargs):
        args = [a(*args, evaluate=True, **kwargs) if isinstance(a, OpExpr) else a for a in self.args]
        return Pow(*args)

    def _eval_power(self, other):
        from sympy import arg, exp, floor, im, log, re, sign
        b, e = self.as_base_exp()
        if b is S.NaN:
            return (b**e)**other  # let __new__ handle it

        s = None
        if other.is_integer:
            s = 1
        elif b.is_polar:  # e.g. exp_polar, besselj, var('p', polar=True)...
            s = 1
        elif e.is_extended_real is not None:
            # helper functions ===========================
            def _half(e):
                """Return True if the exponent has a literal 2 as the
                denominator, else None."""
                if getattr(e, 'q', None) == 2:
                    return True
                n, d = e.as_numer_denom()
                if n.is_integer and d == 2:
                    return True
            def _n2(e):
                """Return ``e`` evaluated to a Number with 2 significant
                digits, else None."""
                try:
                    rv = e.evalf(2, strict=True)
                    if rv.is_Number:
                        return rv
                except PrecisionExhausted:
                    pass
            # ===================================================
            if e.is_extended_real:
                # we need _half(other) with constant floor or
                # floor(S.Half - e*arg(b)/2/pi) == 0

                # handle -1 as special case
                if e == -1:
                    # floor arg. is 1/2 + arg(b)/2/pi
                    if _half(other):
                        if b.is_negative is True:
                            return S.NegativeOne**other*OpPow(-b, e*other)
                        elif b.is_negative is False:
                            return OpPow(b, -other)
                elif e.is_even:
                    if b.is_extended_real:
                        b = abs(b)
                    if b.is_imaginary:
                        b = abs(im(b))*S.ImaginaryUnit

                if (abs(e) < 1) == True or e == 1:
                    s = 1  # floor = 0
                elif b.is_extended_nonnegative:
                    s = 1  # floor = 0
                elif re(b).is_extended_nonnegative and (abs(e) < 2) == True:
                    s = 1  # floor = 0
                elif fuzzy_not(im(b).is_zero) and abs(e) == 2:
                    s = 1  # floor = 0
                elif _half(other):
                    s = exp(2*S.Pi*S.ImaginaryUnit*other*floor(
                        S.Half - e*arg(b)/(2*S.Pi)))
                    if s.is_extended_real and _n2(sign(s) - s) == 0:
                        s = sign(s)
                    else:
                        s = None
            else:
                # e.is_extended_real is False requires:
                #     _half(other) with constant floor or
                #     floor(S.Half - im(e*log(b))/2/pi) == 0
                try:
                    s = exp(2*S.ImaginaryUnit*S.Pi*other*
                        floor(S.Half - im(e*log(b))/2/S.Pi))
                    # be careful to test that s is -1 or 1 b/c sign(I) == I:
                    # so check that s is real
                    if s.is_extended_real and _n2(sign(s) - s) == 0:
                        s = sign(s)
                    else:
                        s = None
                except PrecisionExhausted:
                    s = None

        if s is not None:
            return s*OpPow(b, e*other)

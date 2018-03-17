from __future__ import print_function, division

from collections import defaultdict
from functools import cmp_to_key
import operator

from sympy.multipledispatch import dispatch
from sympy.calculus.util import AccumBounds
from sympy import Basic, Expr, Number, Order, Add, Mul, Pow, S, Symbol
from sympy.matrices.expressions import MatrixExpr
from sympy.core.numbers import ComplexInfinity, ImaginaryUnit, Rational
from sympy.core.logic import fuzzy_not, _fuzzy_group
from sympy.core.evaluate import global_distribute


# internal marker to indicate:
#   "there are still non-commutative objects -- don't forget to process them"
NC_Marker = Symbol("NC_Marker", commutative=False, dummy=True)


class MulBuilder(object):
    def __init__(self, seq):
        self.c_part = []         # out: commutative factors
        self.nc_part = []        # out: non-commutative factors

        self.nc_seq = []

        self.coeff = S.One       # standalone term
                            # e.g. 3 * ...

        self.c_powers = []       # (base,exp)      n
                            # e.g. (x,n) for x

        self.num_exp = []        # (num-base, exp)           y
                            # e.g.  (3, y)  for  ... * 3  * ...

        self.neg1e = S.Zero      # exponent on -1 extracted from Number-based Pow and I

        self.pnum_rat = {}       # (num-base, Rat-exp)          1/2
                            # e.g.  (3, 1/2)  for  ... * 3     * ...

        self.order_symbols = None
        self.seq = list(seq)


# O(x)
@dispatch(MulBuilder, Order)
def append_arg_Mul(data, o):
    o, data.order_symbols = o.as_expr_variables(data.order_symbols)
    append_arg_Mul(data, o)

# Mul([...])
@dispatch(MulBuilder, Mul)
def append_arg_Mul(data, o):

    if o.is_commutative:
        data.seq.extend(o.args)    # XXX zerocopy?

    else:
        # NCMul can have commutative parts as well
        for q in o.args:
            if q.is_commutative:
                data.seq.append(q)
            else:
                data.nc_seq.append(q)

        # append non-commutative marker, so we don't forget to
        # process scheduled non-commutative objects
        data.seq.append(NC_Marker)

@dispatch(MulBuilder, Number)
def append_arg_Mul(data, o):
    # 3
    if o is S.NaN or data.coeff is S.ComplexInfinity and o is S.Zero:
        # we know for sure the result will be nan
        return S.NaN
    elif data.coeff.is_Number:  # it could be zoo
        data.coeff *= o
        if data.coeff is S.NaN:
            # we know for sure the result will be nan
            return S.NaN

@dispatch(MulBuilder, AccumBounds)
def append_arg_Mul(data, o):
    data.coeff = o.__mul__(data.coeff)

@dispatch(MulBuilder, MatrixExpr)
def append_arg_Mul(data, o):
    data.coeff = o.__mul__(data.coeff)

@dispatch(MulBuilder, ComplexInfinity)
def append_arg_Mul(data, o):
    if not data.coeff:
        # 0 * zoo = NaN
        return S.NaN
    if data.coeff is S.ComplexInfinity:
        # zoo * zoo = zoo
        return S.ComplexInfinity
    data.coeff = S.ComplexInfinity
    data.coeff = data.coeff

@dispatch(MulBuilder, ImaginaryUnit)
def append_arg_Mul(data, o):
    data.neg1e += S.Half

@dispatch(MulBuilder, Expr)
def append_arg_Mul(data, o):

    if o.is_commutative:
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
                        data.coeff *= Pow(b, e)  # it is an unevaluated power
                        return
                    elif e.is_negative:    # also a sign of an unevaluated power
                        data.seq.append(Pow(b, e))
                        return
                    elif b.is_negative:
                        data.neg1e += e
                        b = -b
                    if b is not S.One:
                        data.pnum_rat.setdefault(b, []).append(e)
                    return
                elif b.is_positive or e.is_integer:
                    data.num_exp.append((b, e))
                    return

            elif b is S.ImaginaryUnit and e.is_Rational:
                data.neg1e += e/2
                return

        data.c_powers.append((b, e))

    # NON-COMMUTATIVE
    # TODO: Make non-commutative exponents not combine automatically
    else:
        if o is not NC_Marker:
            data.nc_seq.append(o)

        # process data.nc_seq (if any)
        while data.nc_seq:
            o = data.nc_seq.pop(0)
            if not data.nc_part:
                data.nc_part.append(o)
                continue

            #                             b    c       b+c
            # try to combine last terms: a  * a   ->  a
            o1 = data.nc_part.pop()
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
                    data.seq.append(o12)
                    continue
                else:
                    data.nc_seq.insert(0, o12)

            else:
                data.nc_part.append(o1)
                data.nc_part.append(o)

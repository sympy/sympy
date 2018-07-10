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
from sympy.core.expr import _imul
from sympy.core.mul import MutableMul


# internal marker to indicate:
#   "there are still non-commutative objects -- don't forget to process them"
NC_Marker = Symbol("NC_Marker", commutative=False, dummy=True)


# O(x)
@dispatch(MutableMul, Order)
def _imul(data, o):
    o, data.order_symbols = o.as_expr_variables(data.order_symbols)
    _imul(data, o)

# Mul([...])
@dispatch(MutableMul, Mul)
def _imul(data, o):

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

@dispatch(MutableMul, Number)
def _imul(data, o):
    # 3
    if o is S.NaN or data.coeff is S.ComplexInfinity and o is S.Zero:
        # we know for sure the result will be nan
        return S.NaN
    elif data.coeff.is_Number:  # it could be zoo
        data.coeff *= o
        if data.coeff is S.NaN:
            # we know for sure the result will be nan
            return S.NaN

@dispatch(MutableMul, AccumBounds)
def _imul(data, o):
    data.coeff = o.__mul__(data.coeff)

@dispatch(MutableMul, MatrixExpr)
def _imul(data, o):
    data.coeff = o.__mul__(data.coeff)

@dispatch(MutableMul, ComplexInfinity)
def _imul(data, o):
    if not data.coeff:
        # 0 * zoo = NaN
        return S.NaN
    if data.coeff is S.ComplexInfinity:
        # zoo * zoo = zoo
        return S.ComplexInfinity
    data.coeff = S.ComplexInfinity
    data.coeff = data.coeff

@dispatch(MutableMul, ImaginaryUnit)
def _imul(data, o):
    data.neg1e += S.Half

@dispatch(MutableMul, Expr)
def _imul(data, o):

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

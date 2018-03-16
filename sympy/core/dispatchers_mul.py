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


# O(x)
@dispatch(dict, Order)
def append_arg_Mul(data, o):
    order_symbols = data["order_symbols"]
    o, order_symbols = o.as_expr_variables(order_symbols)
    data["order_symbols"] = order_symbols
    append_arg_Mul(data, o)

# Mul([...])
@dispatch(dict, Mul)
def append_arg_Mul(data, o):
    seq = data["seq"]
    nc_seq = data["nc_seq"]

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

    data["seq"] = seq
    data["nc_seq"] = nc_seq

@dispatch(dict, Number)
def append_arg_Mul(data, o):
    # 3
    coeff = data["coeff"]
    if o is S.NaN or coeff is S.ComplexInfinity and o is S.Zero:
        # we know for sure the result will be nan
        return S.NaN
    elif coeff.is_Number:  # it could be zoo
        coeff *= o
        if coeff is S.NaN:
            # we know for sure the result will be nan
            return S.NaN
    data["coeff"] = coeff

@dispatch(dict, AccumBounds)
def append_arg_Mul(data, o):
    coeff = data["coeff"]
    coeff = o.__mul__(coeff)
    data["coeff"] = coeff

@dispatch(dict, MatrixExpr)
def append_arg_Mul(data, o):
    coeff = data["coeff"]
    coeff = o.__mul__(coeff)
    data["coeff"] = coeff

@dispatch(dict, ComplexInfinity)
def append_arg_Mul(data, o):
    coeff = data["coeff"]
    if not coeff:
        # 0 * zoo = NaN
        return S.NaN
    if coeff is S.ComplexInfinity:
        # zoo * zoo = zoo
        return S.ComplexInfinity
    coeff = S.ComplexInfinity
    data["coeff"] = coeff

@dispatch(dict, ImaginaryUnit)
def append_arg_Mul(data, o):
    data["neg1e"] += S.Half

@dispatch(dict, Expr)
def append_arg_Mul(data, o):
    coeff = data["coeff"]
    c_powers = data["c_powers"]
    neg1e = data["neg1e"]
    nc_seq = data["nc_seq"]
    nc_part = data["nc_part"]
    seq = data["seq"]

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
                        data["coeff"] *= Pow(b, e)  # it is an unevaluated power
                        return
                    elif e.is_negative:    # also a sign of an unevaluated power
                        data["seq"].append(Pow(b, e))
                        return
                    elif b.is_negative:
                        data["neg1e"] += e
                        b = -b
                    if b is not S.One:
                        data["pnum_rat"].setdefault(b, []).append(e)
                    return
                elif b.is_positive or e.is_integer:
                    data["num_exp"].append((b, e))
                    return

            elif b is S.ImaginaryUnit and e.is_Rational:
                data["neg1e"] += e/2
                return

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

    data["coeff"] = coeff
    data["c_powers"] = c_powers
    data["neg1e"] = neg1e
    data["nc_seq"] = nc_seq
    data["nc_part"] = nc_part
    data["seq"] = seq

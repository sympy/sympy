from sympy.multipledispatch import dispatch
from sympy.calculus.util import AccumBounds
from sympy import Basic, Expr, Number, Order, Add, Mul, Pow, S
from sympy.matrices.expressions import MatrixExpr
from sympy.core.numbers import ComplexInfinity


def append_arg_c_s(klass, data, c, s):
    terms = data["terms"]
    # now we have:
    # o = c*s, where
    #
    # c is a Number
    # s is an expression with number factor extracted
    # let's collect terms with the same s, so e.g.
    # 2*x**2 + 3*x**2  ->  5*x**2
    if s in terms:
        terms[s] += c
        if terms[s] is S.NaN:
            # we know for sure the result will be nan
            return S.NaN
    else:
        terms[s] = c
    data["terms"] = terms
    return klass, data


@dispatch(type, dict, Expr)
def append_arg_Add(klass, data, o):
    # everything else
    c = S.One
    s = o
    return append_arg_c_s(klass, data, c, s)

@dispatch(type, dict, Order)
def append_arg_Add(klass, data, o):
    order_factors = data["order_factors"]
    for o1 in order_factors:
        if o1.contains(o):
            o = None
            break
    if o is None:
        return klass, data
    order_factors = [o] + [
        o1 for o1 in order_factors if not o.contains(o1)]
    data["order_factors"] = order_factors
    return klass, data

@dispatch(type, dict, Number)
def append_arg_Add(klass, data, o):
    coeff = data["coeff"]
    # 3 or NaN
    if (o is S.NaN or coeff is S.ComplexInfinity and
            o.is_finite is False):
        # we know for sure the result will be nan
        return S.NaN
    if coeff.is_Number:
        coeff += o
        if coeff is S.NaN:
            # we know for sure the result will be nan
            return S.NaN
    data["coeff"] = coeff
    return klass, data

@dispatch(type, dict, AccumBounds)
def append_arg_Add(klass, data, o):
    data["coeff"] = o.__add__(data["coeff"])
    return klass, data

@dispatch(type, dict, MatrixExpr)
def append_arg_Add(klass, data, o):
    # can't add 0 to Matrix so make sure coeff is not 0
    coeff = data["coeff"]
    data["coeff"] = o.__add__(coeff) if coeff else o
    return klass, data

@dispatch(type, dict, ComplexInfinity)
def append_arg_Add(klass, data, o):
    coeff = data["coeff"]
    if coeff.is_finite is False:
        # we know for sure the result will be nan
        return S.NaN
    data["coeff"] = S.ComplexInfinity
    return klass, data

# Add([...])
@dispatch(type, dict, Add)
def append_arg_Add(klass, data, o):
    # NB: here we assume Add is always commutative
    # seq.extend(o.args)  # TODO zerocopy?
    for arg in o.args:
        ret = append_arg_Add(klass, data, arg)
        if not isinstance(ret, tuple):
            return ret
        klass, data = ret
    return klass, data

# Mul([...])
@dispatch(type, dict, Mul)
def append_arg_Add(klass, data, o):
    c, s = o.as_coeff_Mul()
    return append_arg_c_s(klass, data, c, s)

@dispatch(type, dict, Pow)
def append_arg_Add(klass, data, o):
    b, e = o.as_base_exp()
    # check for unevaluated Pow, e.g. 2**3 or 2**(-1/2)
    if b.is_Number and (e.is_Integer or
                       (e.is_Rational and e.is_negative)):
        ret = append_arg_Add(klass, data, b**e)
        #f not isinstance(ret, tuple):
        #   return ret
        #klass, data = ret
        return ret
    c, s = S.One, o
    return append_arg_c_s(klass, data, c, s)

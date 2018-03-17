from sympy.multipledispatch import dispatch
from sympy.calculus.util import AccumBounds
from sympy import Basic, Expr, Number, Order, Add, Mul, Pow, S
from sympy.matrices.expressions import MatrixExpr
from sympy.core.numbers import ComplexInfinity


class AddBuilder(object):
    def __init__(self):
        self.terms = {}
        self.coeff = S.Zero
        self.order_factors = []


def append_arg_c_s(data, c, s):
    # now we have:
    # o = c*s, where
    #
    # c is a Number
    # s is an expression with number factor extracted
    # let's collect terms with the same s, so e.g.
    # 2*x**2 + 3*x**2  ->  5*x**2
    if s in data.terms:
        data.terms[s] += c
        if data.terms[s] is S.NaN:
            # we know for sure the result will be nan
            return S.NaN
    else:
        data.terms[s] = c


@dispatch(AddBuilder, Expr)
def append_arg_Add(data, o):
    # everything else
    c = S.One
    s = o
    return append_arg_c_s(data, c, s)

@dispatch(AddBuilder, Order)
def append_arg_Add(data, o):
    for o1 in data.order_factors:
        if o1.contains(o):
            o = None
            break
    if o is None:
        return
    data.order_factors = [o] + [
        o1 for o1 in data.order_factors if not o.contains(o1)]

@dispatch(AddBuilder, Number)
def append_arg_Add(data, o):
    # 3 or NaN
    if (o is S.NaN or data.coeff is S.ComplexInfinity and
            o.is_finite is False):
        # we know for sure the result will be nan
        return S.NaN
    if data.coeff.is_Number:
        data.coeff += o
        if data.coeff is S.NaN:
            # we know for sure the result will be nan
            return S.NaN

@dispatch(AddBuilder, AccumBounds)
def append_arg_Add(data, o):
    data.coeff = o.__add__(data.coeff)

@dispatch(AddBuilder, MatrixExpr)
def append_arg_Add(data, o):
    # can't add 0 to Matrix so make sure coeff is not 0
    data.coeff = o.__add__(data.coeff) if data.coeff else o

@dispatch(AddBuilder, ComplexInfinity)
def append_arg_Add(data, o):
    if data.coeff.is_finite is False:
        # we know for sure the result will be nan
        return S.NaN
    data.coeff = S.ComplexInfinity

# Add([...])
@dispatch(AddBuilder, Add)
def append_arg_Add(data, o):
    # NB: here we assume Add is always commutative
    # seq.extend(o.args)  # TODO zerocopy?
    for arg in o.args:
        ret = append_arg_Add(data, arg)
        if ret is not None:
            return ret

# Mul([...])
@dispatch(AddBuilder, Mul)
def append_arg_Add(data, o):
    c, s = o.as_coeff_Mul()
    return append_arg_c_s(data, c, s)

@dispatch(AddBuilder, Pow)
def append_arg_Add(data, o):
    b, e = o.as_base_exp()
    # check for unevaluated Pow, e.g. 2**3 or 2**(-1/2)
    if b.is_Number and (e.is_Integer or
                       (e.is_Rational and e.is_negative)):
        return append_arg_Add(data, b**e)
    c, s = S.One, o
    return append_arg_c_s(data, c, s)

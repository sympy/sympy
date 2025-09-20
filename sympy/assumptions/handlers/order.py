from sympy.assumptions.assume import recursive_ask
'\nHandlers related to order relations: positive, negative, etc.\n'
from sympy.assumptions import Q
from sympy.core import Add, Basic, Expr, Mul, Pow, S
from sympy.core.logic import fuzzy_not, fuzzy_and, fuzzy_or
from sympy.core.numbers import E, ImaginaryUnit, NaN, I, pi
from sympy.functions import Abs, acos, acot, asin, atan, exp, factorial, log
from sympy.matrices import Determinant, Trace
from sympy.matrices.expressions.matexpr import MatrixElement
from sympy.multipledispatch import MDNotImplementedError
from ..predicates.order import NegativePredicate, NonNegativePredicate, NonZeroPredicate, ZeroPredicate, NonPositivePredicate, PositivePredicate, ExtendedNegativePredicate, ExtendedNonNegativePredicate, ExtendedNonPositivePredicate, ExtendedNonZeroPredicate, ExtendedPositivePredicate

def _NegativePredicate_number(expr, assumptions):
    r, i = expr.as_real_imag()
    if r == S.NaN or i == S.NaN:
        return None
    if not i:
        r = r.evalf(2)
        if r._prec != 1:
            return r < 0
    else:
        i = i.evalf(2)
        if i._prec != 1:
            if i != 0:
                return False
            r = r.evalf(2)
            if r._prec != 1:
                return r < 0

@NegativePredicate.register(Basic)
def _(expr, assumptions, rec):
    if expr.is_number:
        return _NegativePredicate_number(expr, assumptions)

@NegativePredicate.register(Expr)
def _(expr, assumptions, rec):
    ret = expr.is_negative
    if ret is None:
        raise MDNotImplementedError
    return ret

@NegativePredicate.register(Add)
def _(expr, assumptions, rec):
    """
    Positive + Positive -> Positive,
    Negative + Negative -> Negative
    """
    if expr.is_number:
        return _NegativePredicate_number(expr, assumptions)
    r = recursive_ask(Q.real(expr), assumptions=assumptions, rec=rec)
    if r is not True:
        return r
    nonpos = 0
    for arg in expr.args:
        if recursive_ask(Q.negative(arg), assumptions=assumptions, rec=rec) is not True:
            if recursive_ask(Q.positive(arg), assumptions=assumptions, rec=rec) is False:
                nonpos += 1
            else:
                break
    else:
        if nonpos < len(expr.args):
            return True

@NegativePredicate.register(Mul)
def _(expr, assumptions, rec):
    if expr.is_number:
        return _NegativePredicate_number(expr, assumptions)
    result = None
    for arg in expr.args:
        if result is None:
            result = False
        if recursive_ask(Q.negative(arg), assumptions=assumptions, rec=rec):
            result = not result
        elif recursive_ask(Q.positive(arg), assumptions=assumptions, rec=rec):
            pass
        else:
            return
    return result

@NegativePredicate.register(Pow)
def _(expr, assumptions, rec):
    """
    Real ** Even -> NonNegative
    Real ** Odd  -> same_as_base
    NonNegative ** Positive -> NonNegative
    """
    if expr.base == E:
        if recursive_ask(Q.real(expr.exp), assumptions=assumptions, rec=rec):
            return False
        return
    if expr.is_number:
        return _NegativePredicate_number(expr, assumptions)
    if recursive_ask(Q.real(expr.base), assumptions=assumptions, rec=rec):
        if recursive_ask(Q.positive(expr.base), assumptions=assumptions, rec=rec):
            if recursive_ask(Q.real(expr.exp), assumptions=assumptions, rec=rec):
                return False
        if recursive_ask(Q.even(expr.exp), assumptions=assumptions, rec=rec):
            return False
        if recursive_ask(Q.odd(expr.exp), assumptions=assumptions, rec=rec):
            return recursive_ask(Q.negative(expr.base), assumptions=assumptions, rec=rec)

@NegativePredicate.register_many(Abs, ImaginaryUnit)
def _(expr, assumptions, rec):
    return False

@NegativePredicate.register(exp)
def _(expr, assumptions, rec):
    if recursive_ask(Q.real(expr.exp), assumptions=assumptions, rec=rec):
        return False
    raise MDNotImplementedError

@NonNegativePredicate.register(Basic)
def _(expr, assumptions, rec):
    if expr.is_number:
        notnegative = fuzzy_not(_NegativePredicate_number(expr, assumptions))
        if notnegative:
            return recursive_ask(Q.real(expr), assumptions=assumptions, rec=rec)
        else:
            return notnegative

@NonNegativePredicate.register(Expr)
def _(expr, assumptions, rec):
    ret = expr.is_nonnegative
    if ret is None:
        raise MDNotImplementedError
    return ret

@NonZeroPredicate.register(Expr)
def _(expr, assumptions, rec):
    ret = expr.is_nonzero
    if ret is None:
        raise MDNotImplementedError
    return ret

@NonZeroPredicate.register(Basic)
def _(expr, assumptions, rec):
    if recursive_ask(Q.real(expr), assumptions=assumptions, rec=rec) is False:
        return False
    if expr.is_number:
        i = expr.evalf(2)

        def nonz(i):
            if i._prec != 1:
                return i != 0
        return fuzzy_or((nonz(i) for i in i.as_real_imag()))

@NonZeroPredicate.register(Add)
def _(expr, assumptions, rec):
    if all((recursive_ask(Q.positive(x), assumptions=assumptions, rec=rec) for x in expr.args)) or all((recursive_ask(Q.negative(x), assumptions=assumptions, rec=rec) for x in expr.args)):
        return True

@NonZeroPredicate.register(Mul)
def _(expr, assumptions, rec):
    for arg in expr.args:
        result = recursive_ask(Q.nonzero(arg), assumptions=assumptions, rec=rec)
        if result:
            continue
        return result
    return True

@NonZeroPredicate.register(Pow)
def _(expr, assumptions, rec):
    return recursive_ask(Q.nonzero(expr.base), assumptions=assumptions, rec=rec)

@NonZeroPredicate.register(Abs)
def _(expr, assumptions, rec):
    return recursive_ask(Q.nonzero(expr.args[0]), assumptions=assumptions, rec=rec)

@NonZeroPredicate.register(NaN)
def _(expr, assumptions, rec):
    return None

@ZeroPredicate.register(Expr)
def _(expr, assumptions, rec):
    ret = expr.is_zero
    if ret is None:
        raise MDNotImplementedError
    return ret

@ZeroPredicate.register(Basic)
def _(expr, assumptions, rec):
    return fuzzy_and([fuzzy_not(recursive_ask(Q.nonzero(expr), assumptions=assumptions, rec=rec)), recursive_ask(Q.real(expr), assumptions=assumptions, rec=rec)])

@ZeroPredicate.register(Mul)
def _(expr, assumptions, rec):
    return fuzzy_or((recursive_ask(Q.zero(arg), assumptions=assumptions, rec=rec) for arg in expr.args))

@NonPositivePredicate.register(Expr)
def _(expr, assumptions, rec):
    ret = expr.is_nonpositive
    if ret is None:
        raise MDNotImplementedError
    return ret

@NonPositivePredicate.register(Basic)
def _(expr, assumptions, rec):
    if expr.is_number:
        notpositive = fuzzy_not(_PositivePredicate_number(expr, assumptions))
        if notpositive:
            return recursive_ask(Q.real(expr), assumptions=assumptions, rec=rec)
        else:
            return notpositive

def _PositivePredicate_number(expr, assumptions):
    r, i = expr.as_real_imag()
    if not i:
        r = r.evalf(2)
        if r._prec != 1:
            return r > 0
    else:
        i = i.evalf(2)
        if i._prec != 1:
            if i != 0:
                return False
            r = r.evalf(2)
            if r._prec != 1:
                return r > 0

@PositivePredicate.register(Expr)
def _(expr, assumptions, rec):
    ret = expr.is_positive
    if ret is None:
        raise MDNotImplementedError
    return ret

@PositivePredicate.register(Basic)
def _(expr, assumptions, rec):
    if expr.is_number:
        return _PositivePredicate_number(expr, assumptions)

@PositivePredicate.register(Mul)
def _(expr, assumptions, rec):
    if expr.is_number:
        return _PositivePredicate_number(expr, assumptions)
    result = True
    for arg in expr.args:
        if recursive_ask(Q.positive(arg), assumptions=assumptions, rec=rec):
            continue
        elif recursive_ask(Q.negative(arg), assumptions=assumptions, rec=rec):
            result = result ^ True
        else:
            return
    return result

@PositivePredicate.register(Add)
def _(expr, assumptions, rec):
    if expr.is_number:
        return _PositivePredicate_number(expr, assumptions)
    r = recursive_ask(Q.real(expr), assumptions=assumptions, rec=rec)
    if r is not True:
        return r
    nonneg = 0
    for arg in expr.args:
        if recursive_ask(Q.positive(arg), assumptions=assumptions, rec=rec) is not True:
            if recursive_ask(Q.negative(arg), assumptions=assumptions, rec=rec) is False:
                nonneg += 1
            else:
                break
    else:
        if nonneg < len(expr.args):
            return True

@PositivePredicate.register(Pow)
def _(expr, assumptions, rec):
    if expr.base == E:
        if recursive_ask(Q.real(expr.exp), assumptions=assumptions, rec=rec):
            return True
        if recursive_ask(Q.imaginary(expr.exp), assumptions=assumptions, rec=rec):
            return recursive_ask(Q.even(expr.exp / (I * pi)), assumptions=assumptions, rec=rec)
        return
    if expr.is_number:
        return _PositivePredicate_number(expr, assumptions)
    if recursive_ask(Q.positive(expr.base), assumptions=assumptions, rec=rec):
        if recursive_ask(Q.real(expr.exp), assumptions=assumptions, rec=rec):
            return True
    if recursive_ask(Q.negative(expr.base), assumptions=assumptions, rec=rec):
        if recursive_ask(Q.even(expr.exp), assumptions=assumptions, rec=rec):
            return True
        if recursive_ask(Q.odd(expr.exp), assumptions=assumptions, rec=rec):
            return False

@PositivePredicate.register(exp)
def _(expr, assumptions, rec):
    if recursive_ask(Q.real(expr.exp), assumptions=assumptions, rec=rec):
        return True
    if recursive_ask(Q.imaginary(expr.exp), assumptions=assumptions, rec=rec):
        return recursive_ask(Q.even(expr.exp / (I * pi)), assumptions=assumptions, rec=rec)

@PositivePredicate.register(log)
def _(expr, assumptions, rec):
    r = recursive_ask(Q.real(expr.args[0]), assumptions=assumptions, rec=rec)
    if r is not True:
        return r
    if recursive_ask(Q.positive(expr.args[0] - 1), assumptions=assumptions, rec=rec):
        return True
    if recursive_ask(Q.negative(expr.args[0] - 1), assumptions=assumptions, rec=rec):
        return False

@PositivePredicate.register(factorial)
def _(expr, assumptions, rec):
    x = expr.args[0]
    if recursive_ask(Q.integer(x) & Q.positive(x), assumptions=assumptions, rec=rec):
        return True

@PositivePredicate.register(ImaginaryUnit)
def _(expr, assumptions, rec):
    return False

@PositivePredicate.register(Abs)
def _(expr, assumptions, rec):
    return recursive_ask(Q.nonzero(expr), assumptions=assumptions, rec=rec)

@PositivePredicate.register(Trace)
def _(expr, assumptions, rec):
    if recursive_ask(Q.positive_definite(expr.arg), assumptions=assumptions, rec=rec):
        return True

@PositivePredicate.register(Determinant)
def _(expr, assumptions, rec):
    if recursive_ask(Q.positive_definite(expr.arg), assumptions=assumptions, rec=rec):
        return True

@PositivePredicate.register(MatrixElement)
def _(expr, assumptions, rec):
    if expr.i == expr.j and recursive_ask(Q.positive_definite(expr.parent), assumptions=assumptions, rec=rec):
        return True

@PositivePredicate.register(atan)
def _(expr, assumptions, rec):
    return recursive_ask(Q.positive(expr.args[0]), assumptions=assumptions, rec=rec)

@PositivePredicate.register(asin)
def _(expr, assumptions, rec):
    x = expr.args[0]
    if recursive_ask(Q.positive(x) & Q.nonpositive(x - 1), assumptions=assumptions, rec=rec):
        return True
    if recursive_ask(Q.negative(x) & Q.nonnegative(x + 1), assumptions=assumptions, rec=rec):
        return False

@PositivePredicate.register(acos)
def _(expr, assumptions, rec):
    x = expr.args[0]
    if recursive_ask(Q.nonpositive(x - 1) & Q.nonnegative(x + 1), assumptions=assumptions, rec=rec):
        return True

@PositivePredicate.register(acot)
def _(expr, assumptions, rec):
    return recursive_ask(Q.real(expr.args[0]), assumptions=assumptions, rec=rec)

@PositivePredicate.register(NaN)
def _(expr, assumptions, rec):
    return None

@ExtendedNegativePredicate.register(object)
def _(expr, assumptions, rec):
    return recursive_ask(Q.negative(expr) | Q.negative_infinite(expr), assumptions=assumptions, rec=rec)

@ExtendedPositivePredicate.register(object)
def _(expr, assumptions, rec):
    return recursive_ask(Q.positive(expr) | Q.positive_infinite(expr), assumptions=assumptions, rec=rec)

@ExtendedNonZeroPredicate.register(object)
def _(expr, assumptions, rec):
    return recursive_ask(Q.negative_infinite(expr) | Q.negative(expr) | Q.positive(expr) | Q.positive_infinite(expr), assumptions=assumptions, rec=rec)

@ExtendedNonPositivePredicate.register(object)
def _(expr, assumptions, rec):
    return recursive_ask(Q.negative_infinite(expr) | Q.negative(expr) | Q.zero(expr), assumptions=assumptions, rec=rec)

@ExtendedNonNegativePredicate.register(object)
def _(expr, assumptions, rec):
    return recursive_ask(Q.zero(expr) | Q.positive(expr) | Q.positive_infinite(expr), assumptions=assumptions, rec=rec)
"""
AskHandlers related to order relations: positive, negative, etc.
"""

from sympy.assumptions import Q, ask
from sympy.assumptions.handlers import AskHandlerClass, CommonHandler
from sympy.core.logic import fuzzy_not, fuzzy_and, fuzzy_or
from sympy.core import (
    Expr, Basic, Add, Mul, Pow
)
from sympy.core.numbers import ImaginaryUnit, NaN
from sympy.functions import (
    Abs, exp, log, factorial, atan, asin, acos, acot
)
from sympy.matrices import Trace, Determinant
from sympy.matrices.expressions.matexpr import MatrixElement


### AskNegativeHandler ###

AskNegativeHandler = AskHandlerClass(
    'AskNegativeHandler',
    doc="""
    This is called by ask() when key='negative'

    Test that an expression is less (strict) than zero.

    Examples
    ========

    >>> from sympy import ask, Q, pi
    >>> ask(Q.negative(pi+1)) # this calls AskNegativeHandler.Add
    False
    >>> ask(Q.negative(pi**2)) # this calls AskNegativeHandler.Pow
    False

    """,
    base=CommonHandler
)

for sig in (ImaginaryUnit, Abs):
    AskNegativeHandler.register(sig)(AskNegativeHandler.AlwaysFalse)

@AskNegativeHandler.register(Expr)
def _(expr, assumptions):
    return expr.is_negative

def _number_negative(expr, assumptions):
    r, i = expr.as_real_imag()
    # If the imaginary part can symbolically be shown to be zero then
    # we just evaluate the real part; otherwise we evaluate the imaginary
    # part to see if it actually evaluates to zero and if it does then
    # we make the comparison between the real part and zero.
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

@AskNegativeHandler.register(Basic)
def _(expr, assumptions):
    if expr.is_number:
        return _number_negative(expr, assumptions)

@AskNegativeHandler.register(Add)
def _(expr, assumptions):
    """
    Positive + Positive -> Positive,
    Negative + Negative -> Negative
    """
    if expr.is_number:
        return _number_negative(expr, assumptions)

    r = ask(Q.real(expr), assumptions)
    if r is not True:
        return r

    nonpos = 0
    for arg in expr.args:
        if ask(Q.negative(arg), assumptions) is not True:
            if ask(Q.positive(arg), assumptions) is False:
                nonpos += 1
            else:
                break
    else:
        if nonpos < len(expr.args):
            return True

@AskNegativeHandler.register(Mul)
def _(expr, assumptions):
    if expr.is_number:
        return _number_negative(expr, assumptions)
    result = None
    for arg in expr.args:
        if result is None:
            result = False
        if ask(Q.negative(arg), assumptions):
            result = not result
        elif ask(Q.positive(arg), assumptions):
            pass
        else:
            return
    return result

@AskNegativeHandler.register(Pow)
def _(expr, assumptions):
    """
    Real ** Even -> NonNegative
    Real ** Odd  -> same_as_base
    NonNegative ** Positive -> NonNegative
    """
    if expr.is_number:
        return _number_negative(expr, assumptions)
    if ask(Q.real(expr.base), assumptions):
        if ask(Q.positive(expr.base), assumptions):
            if ask(Q.real(expr.exp), assumptions):
                return False
        if ask(Q.even(expr.exp), assumptions):
            return False
        if ask(Q.odd(expr.exp), assumptions):
            return ask(Q.negative(expr.base), assumptions)

@AskNegativeHandler.register(exp)
def _(expr, assumptions):
    if ask(Q.real(expr.args[0]), assumptions):
        return False


### AskNonNegativeHandler ###

AskNonNegativeHandler = AskHandlerClass(
    'AskNonNegativeHandler',
    doc="""
    Handler for key 'nonnegative'
    """,
    base=CommonHandler
)

@AskNonNegativeHandler.register(Expr)
def _(expr, assumptions):
    return expr.is_nonnegative

@AskNonNegativeHandler.register(Basic)
def _(expr, assumptions):
    if expr.is_number:
        notnegative = fuzzy_not(_number_negative(expr, assumptions))
        if notnegative:
            return ask(Q.real(expr), assumptions)
        else:
            return notnegative


### AskNonZeroHandler ###

AskNonZeroHandler = AskHandlerClass(
    'AskNonZeroHandler',
    doc="""
    Handler for key 'nonzero'
    Test that an expression is not identically zero
    """,
    base=CommonHandler
)

@AskNonZeroHandler.register(Expr)
def _(expr, assumptions):
    return expr.is_nonzero

@AskNonZeroHandler.register(Basic)
def _(expr, assumptions):
    if ask(Q.real(expr)) is False:
        return False
    if expr.is_number:
        # if there are no symbols just evalf
        i = expr.evalf(2)
        def nonz(i):
            if i._prec != 1:
                return i != 0
        return fuzzy_or(nonz(i) for i in i.as_real_imag())

@AskNonZeroHandler.register(Add)
def _(expr, assumptions):
    if all(ask(Q.positive(x), assumptions) for x in expr.args) \
            or all(ask(Q.negative(x), assumptions) for x in expr.args):
        return True

@AskNonZeroHandler.register(Mul)
def _(expr, assumptions):
    for arg in expr.args:
        result = ask(Q.nonzero(arg), assumptions)
        if result:
            continue
        return result
    return True

@AskNonZeroHandler.register(Pow)
def _(expr, assumptions):
    return ask(Q.nonzero(expr.base), assumptions)

@AskNonZeroHandler.register(NaN)
def _(expr, assumptions):
    return True

@AskNonZeroHandler.register(Abs)
def _(expr, assumptions):
    return ask(Q.nonzero(expr.args[0]), assumptions)


### AskZeroHandler ###

AskZeroHandler = AskHandlerClass(
    'AskZeroHandler',
    doc="""
    Handler for key 'zero'
    Test that an expression is identically zero
    """,
    base=CommonHandler
)

@AskZeroHandler.register(Expr)
def _(expr, assumptions):
    return expr.is_zero

@AskZeroHandler.register(Basic)
def _(expr, assumptions):
    return fuzzy_and([fuzzy_not(ask(Q.nonzero(expr), assumptions)),
        ask(Q.real(expr), assumptions)])

@AskZeroHandler.register(Mul)
def _(expr, assumptions):
    # TODO: This should be deducible from the nonzero handler
    return fuzzy_or(ask(Q.zero(arg), assumptions) for arg in expr.args)


### AskNonPositiveHandler ###

AskNonPositiveHandler = AskHandlerClass(
    'AskNonPositiveHandler',
    doc="""
    Handler for key 'nonpositive'
    """,
    base=CommonHandler
)

@AskNonPositiveHandler.register(Expr)
def _(expr, assumptions):
    return expr.is_nonpositive

@AskNonPositiveHandler.register(Basic)
def _(expr, assumptions):
    if expr.is_number:
        notpositive = fuzzy_not(_number_positive(expr, assumptions))
        if notpositive:
            return ask(Q.real(expr), assumptions)
        else:
            return notpositive


### AskPositiveHandler ###

AskPositiveHandler = AskHandlerClass(
    'AskPositiveHandler',
    doc="""
    Handler for key 'positive'
    Test that an expression is greater (strict) than zero
    """,
    base=CommonHandler
)

AskPositiveHandler.register(ImaginaryUnit)(AskPositiveHandler.AlwaysFalse)

@AskPositiveHandler.register(Expr)
def _(expr, assumptions):
    return expr.is_positive

def _number_positive(expr, assumptions):
    r, i = expr.as_real_imag()
    # If the imaginary part can symbolically be shown to be zero then
    # we just evaluate the real part; otherwise we evaluate the imaginary
    # part to see if it actually evaluates to zero and if it does then
    # we make the comparison between the real part and zero.
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

@AskPositiveHandler.register(Basic)
def _(expr, assumptions):
    if expr.is_number:
        return _number_positive(expr, assumptions)

@AskPositiveHandler.register(Mul)
def _(expr, assumptions):
    if expr.is_number:
        return _number_positive(expr, assumptions)
    result = True
    for arg in expr.args:
        if ask(Q.positive(arg), assumptions):
            continue
        elif ask(Q.negative(arg), assumptions):
            result = result ^ True
        else:
            return
    return result

@AskPositiveHandler.register(Add)
def _(expr, assumptions):
    if expr.is_number:
        return _number_positive(expr, assumptions)

    r = ask(Q.real(expr), assumptions)
    if r is not True:
        return r

    nonneg = 0
    for arg in expr.args:
        if ask(Q.positive(arg), assumptions) is not True:
            if ask(Q.negative(arg), assumptions) is False:
                nonneg += 1
            else:
                break
    else:
        if nonneg < len(expr.args):
            return True

@AskPositiveHandler.register(Pow)
def _(expr, assumptions):
    if expr.is_number:
        return _number_positive(expr, assumptions)
    if ask(Q.positive(expr.base), assumptions):
        if ask(Q.real(expr.exp), assumptions):
            return True
    if ask(Q.negative(expr.base), assumptions):
        if ask(Q.even(expr.exp), assumptions):
            return True
        if ask(Q.odd(expr.exp), assumptions):
            return False

@AskPositiveHandler.register(exp)
def _(expr, assumptions):
    if ask(Q.real(expr.args[0]), assumptions):
        return True
    if ask(Q.imaginary(expr.args[0]), assumptions):
        from sympy import pi, I
        return ask(Q.even(expr.args[0]/(I*pi)), assumptions)

@AskPositiveHandler.register(log)
def _(expr, assumptions):
    r = ask(Q.real(expr.args[0]), assumptions)
    if r is not True:
        return r
    if ask(Q.positive(expr.args[0] - 1), assumptions):
        return True
    if ask(Q.negative(expr.args[0] - 1), assumptions):
        return False

@AskPositiveHandler.register(factorial)
def _(expr, assumptions):
    x = expr.args[0]
    if ask(Q.integer(x) & Q.positive(x), assumptions):
        return True

@AskPositiveHandler.register(Abs)
def _(expr, assumptions):
    return ask(Q.nonzero(expr), assumptions)

@AskPositiveHandler.register(Trace)
def _(expr, assumptions):
    if ask(Q.positive_definite(expr.arg), assumptions):
        return True

@AskPositiveHandler.register(Determinant)
def _(expr, assumptions):
    if ask(Q.positive_definite(expr.arg), assumptions):
        return True

@AskPositiveHandler.register(MatrixElement)
def _(expr, assumptions):
    if (expr.i == expr.j
            and ask(Q.positive_definite(expr.parent), assumptions)):
        return True

@AskPositiveHandler.register(atan)
def _(expr, assumptions):
    return ask(Q.positive(expr.args[0]), assumptions)

@AskPositiveHandler.register(asin)
def _(expr, assumptions):
    x = expr.args[0]
    if ask(Q.positive(x) & Q.nonpositive(x - 1), assumptions):
        return True
    if ask(Q.negative(x) & Q.nonnegative(x + 1), assumptions):
        return False

@AskPositiveHandler.register(acos)
def _(expr, assumptions):
    x = expr.args[0]
    if ask(Q.nonpositive(x - 1) & Q.nonnegative(x + 1), assumptions):
        return True

@AskPositiveHandler.register(acot)
def _(expr, assumptions):
    return ask(Q.real(expr.args[0]), assumptions)

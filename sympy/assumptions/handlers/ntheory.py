"""
Handlers for keys related to number theory: prime, even, odd, etc.
"""

from sympy.assumptions import Q, ask
from sympy.assumptions.handlers import CommonHandler
from sympy.ntheory import isprime
from sympy.core import (
    S, Expr, Basic, Mul, Pow, Add
)
from sympy.core.numbers import (
    Float, Integer, Rational, Infinity, NegativeInfinity, ImaginaryUnit, NumberSymbol
)
from sympy.functions import Abs, re, im


### AskPrimeHandler ###

AskPrimeHandler = CommonHandler.copy(
    'AskPrimeHandler',
    doc="""
    Handler for key 'prime'
    Test that an expression represents a prime number. When the
    expression is an exact number, the result (when True) is subject to
    the limitations of isprime() which is used to return the result.
    """
)

for sig in (Rational, Infinity, NegativeInfinity, ImaginaryUnit):
    AskPrimeHandler.register(sig)(AskPrimeHandler.AlwaysFalse)

@AskPrimeHandler.register(Expr)
def _(expr, assumptions):
    return expr.is_prime

def _number_prime(expr, assumptions):
    # helper method
    exact = not expr.atoms(Float)
    try:
        i = int(expr.round())
        if (expr - i).equals(0) is False:
            raise TypeError
    except TypeError:
        return False
    if exact:
        return isprime(i)
    # when not exact, we won't give a True or False
    # since the number represents an approximate value

@AskPrimeHandler.register(Basic)
def _(expr, assumptions):
    if expr.is_number:
        return _number_prime(expr, assumptions)

@AskPrimeHandler.register(Mul)
def _(expr, assumptions):
    if expr.is_number:
        return _number_prime(expr, assumptions)
    for arg in expr.args:
        if not ask(Q.integer(arg), assumptions):
            return None
    for arg in expr.args:
        if arg.is_number and arg.is_composite:
            return False

@AskPrimeHandler.register(Pow)
def _(expr, assumptions):
    """
    Integer**Integer     -> !Prime
    """
    if expr.is_number:
        return _number_prime(expr, assumptions)
    if ask(Q.integer(expr.exp), assumptions) and \
            ask(Q.integer(expr.base), assumptions):
        return False

@AskPrimeHandler.register(Integer)
def _(expr, assumptions):
    return isprime(expr)

@AskPrimeHandler.register(Float)
def _(expr, assumptions):
    return _number_prime(expr, assumptions)

@AskPrimeHandler.register(NumberSymbol)
def _(expr, assumptions):
    return _number_prime(expr, assumptions)


### AskCompositeHandler ###

AskCompositeHandler = CommonHandler.copy(
    'AskCompositeHandler',
    doc="""
    Handler for key 'composite'
    """
)

@AskCompositeHandler.register(Expr)
def _(expr, assumptions):
    return expr.is_composite

@AskCompositeHandler.register(Basic)
def _(expr, assumptions):
    _positive = ask(Q.positive(expr), assumptions)
    if _positive:
        _integer = ask(Q.integer(expr), assumptions)
        if _integer:
            _prime = ask(Q.prime(expr), assumptions)
            if _prime is None:
                return
            # Positive integer which is not prime is not
            # necessarily composite
            if expr.equals(1):
                return False
            return not _prime
        else:
            return _integer
    else:
        return _positive


### AskEvenHandler ###

AskEvenHandler = CommonHandler.copy(
    'AskEvenHandler',
    doc="""
    Handler for key 'even'
    """
)

for sig in (Rational, Infinity, NegativeInfinity, ImaginaryUnit):
    AskEvenHandler.register(sig)(AskEvenHandler.AlwaysFalse)

@AskEvenHandler.register(Expr)
def _(expr, assumptions):
    return expr.is_even

def _number_even(expr, assumptions):
    # helper method
    try:
        i = int(expr.round())
        if not (expr - i).equals(0):
            raise TypeError
    except TypeError:
        return False
    if isinstance(expr, (float, Float)):
        return False
    return i % 2 == 0

@AskEvenHandler.register(Basic)
def _(expr, assumptions):
    if expr.is_number:
        return _number_even(expr, assumptions)

@AskEvenHandler.register(Mul)
def _(expr, assumptions):
    """
    Even * Integer    -> Even
    Even * Odd        -> Even
    Integer * Odd     -> ?
    Odd * Odd         -> Odd
    Even * Even       -> Even
    Integer * Integer -> Even if Integer + Integer = Odd
    otherwise         -> ?
    """
    if expr.is_number:
        return _number_even(expr, assumptions)
    even, odd, irrational, acc = False, 0, False, 1
    for arg in expr.args:
        # check for all integers and at least one even
        if ask(Q.integer(arg), assumptions):
            if ask(Q.even(arg), assumptions):
                even = True
            elif ask(Q.odd(arg), assumptions):
                odd += 1
            elif not even and acc != 1:
                if ask(Q.odd(acc + arg), assumptions):
                    even = True
        elif ask(Q.irrational(arg), assumptions):
            # one irrational makes the result False
            # two makes it undefined
            if irrational:
                break
            irrational = True
        else:
            break
        acc = arg
    else:
        if irrational:
            return False
        if even:
            return True
        if odd == len(expr.args):
            return False

@AskEvenHandler.register(Add)
def _(expr, assumptions):
    """
    Even + Odd  -> Odd
    Even + Even -> Even
    Odd  + Odd  -> Even

    """
    if expr.is_number:
        return _number_even(expr, assumptions)
    _result = True
    for arg in expr.args:
        if ask(Q.even(arg), assumptions):
            pass
        elif ask(Q.odd(arg), assumptions):
            _result = not _result
        else:
            break
    else:
        return _result

@AskEvenHandler.register(Pow)
def _(expr, assumptions):
    if expr.is_number:
        return _number_even(expr, assumptions)
    if ask(Q.integer(expr.exp), assumptions):
        if ask(Q.positive(expr.exp), assumptions):
            return ask(Q.even(expr.base), assumptions)
        elif ask(~Q.negative(expr.exp) & Q.odd(expr.base), assumptions):
            return False
        elif expr.base is S.NegativeOne:
            return False

@AskEvenHandler.register(Integer)
def _(expr, assumptions):
    return not bool(expr.p & 1)

@AskEvenHandler.register(NumberSymbol)
def _(expr, assumptions):
    return _number_even(expr, assumptions)

@AskEvenHandler.register(Abs)
def _(expr, assumptions):
    if ask(Q.real(expr.args[0]), assumptions):
        return ask(Q.even(expr.args[0]), assumptions)

@AskEvenHandler.register(re)
def _(expr, assumptions):
    if ask(Q.real(expr.args[0]), assumptions):
        return ask(Q.even(expr.args[0]), assumptions)

@AskEvenHandler.register(im)
def _(expr, assumptions):
    if ask(Q.real(expr.args[0]), assumptions):
        return True


### AskOddHandler ###

AskOddHandler = CommonHandler.copy(
    'AskOddHandler',
    doc="""
    Handler for key 'odd'
    Test that an expression represents an odd number
    """
)

@AskOddHandler.register(Expr)
def _(expr, assumptions):
    return expr.is_odd

@AskOddHandler.register(Basic)
def _(expr, assumptions):
    _integer = ask(Q.integer(expr), assumptions)
    if _integer:
        _even = ask(Q.even(expr), assumptions)
        if _even is None:
            return None
        return not _even
    return _integer

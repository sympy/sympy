"""
Handlers for keys related to number theory: prime, even, odd, etc.
"""

from sympy.assumptions import Q, ask
from sympy.core import Add, Basic, Expr, Float, Mul, Pow, S
from sympy.core.numbers import (ImaginaryUnit, Infinity, Integer, NaN,
    NegativeInfinity, NumberSymbol, Rational)
from sympy.functions import Abs, im, re
from sympy.ntheory import isprime

from sympy.multipledispatch import MDNotImplementedError

from ..predicates.ntheory import (PrimePredicate, CompositePredicate,
    EvenPredicate, OddPredicate)


# PrimePredicate

def _PrimePredicate_number(expr, assumptions):
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

@PrimePredicate.register(Expr) # type: ignore
def _(expr, assumptions):
    ret = expr.is_prime
    if ret is None:
        raise MDNotImplementedError
    return ret

@PrimePredicate.register(Basic) # type: ignore
def _(expr, assumptions):
    if expr.is_number:
        return _PrimePredicate_number(expr, assumptions)

@PrimePredicate.register(Mul) # type: ignore
def _(expr, assumptions):
    if expr.is_number:
        return _PrimePredicate_number(expr, assumptions)
    for arg in expr.args:
        if not ask(Q.integer(arg), assumptions):
            return None
    for arg in expr.args:
        if arg.is_number and arg.is_composite:
            return False

@PrimePredicate.register(Pow) # type: ignore
def _(expr, assumptions):
    """
    Integer**Integer     -> !Prime
    """
    if expr.is_number:
        return _PrimePredicate_number(expr, assumptions)
    if ask(Q.integer(expr.exp), assumptions) and \
            ask(Q.integer(expr.base), assumptions):
        return False

@PrimePredicate.register(Integer) # type: ignore
def _(expr, assumptions):
    return isprime(expr)

@PrimePredicate.register_many(Rational, Infinity, NegativeInfinity, ImaginaryUnit) # type: ignore
def _(expr, assumptions):
    return False

@PrimePredicate.register(Float) # type: ignore
def _(expr, assumptions):
    return _PrimePredicate_number(expr, assumptions)

@PrimePredicate.register(NumberSymbol) # type: ignore
def _(expr, assumptions):
    return _PrimePredicate_number(expr, assumptions)

@PrimePredicate.register(NaN) # type: ignore
def _(expr, assumptions):
    return None


# CompositePredicate

@CompositePredicate.register(Expr) # type: ignore
def _(expr, assumptions):
    ret = expr.is_composite
    if ret is None:
        raise MDNotImplementedError
    return ret

@CompositePredicate.register(Basic) # type: ignore
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


# EvenPredicate

def _EvenPredicate_number(expr, assumptions):
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

@EvenPredicate.register(Expr) # type: ignore
def _(expr, assumptions):
    ret = expr.is_even
    if ret is None:
        raise MDNotImplementedError
    return ret

@EvenPredicate.register(Basic) # type: ignore
def _(expr, assumptions):
    if expr.is_number:
        return _EvenPredicate_number(expr, assumptions)

@EvenPredicate.register(Mul) # type: ignore
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
        return _EvenPredicate_number(expr, assumptions)
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

@EvenPredicate.register(Add) # type: ignore
def _(expr, assumptions):
    """
    Even + Odd  -> Odd
    Even + Even -> Even
    Odd  + Odd  -> Even

    """
    if expr.is_number:
        return _EvenPredicate_number(expr, assumptions)
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

@EvenPredicate.register(Pow) # type: ignore
def _(expr, assumptions):
    if expr.is_number:
        return _EvenPredicate_number(expr, assumptions)
    if ask(Q.integer(expr.exp), assumptions):
        if ask(Q.positive(expr.exp), assumptions):
            return ask(Q.even(expr.base), assumptions)
        elif ask(~Q.negative(expr.exp) & Q.odd(expr.base), assumptions):
            return False
        elif expr.base is S.NegativeOne:
            return False

@EvenPredicate.register(Integer) # type: ignore
def _(expr, assumptions):
    return not bool(expr.p & 1)

@EvenPredicate.register_many(Rational, Infinity, NegativeInfinity, ImaginaryUnit) # type: ignore
def _(expr, assumptions):
    return False

@EvenPredicate.register(NumberSymbol) # type: ignore
def _(expr, assumptions):
    return _EvenPredicate_number(expr, assumptions)

@EvenPredicate.register(Abs) # type: ignore
def _(expr, assumptions):
    if ask(Q.real(expr.args[0]), assumptions):
        return ask(Q.even(expr.args[0]), assumptions)

@EvenPredicate.register(re) # type: ignore
def _(expr, assumptions):
    if ask(Q.real(expr.args[0]), assumptions):
        return ask(Q.even(expr.args[0]), assumptions)

@EvenPredicate.register(im) # type: ignore
def _(expr, assumptions):
    if ask(Q.real(expr.args[0]), assumptions):
        return True

@EvenPredicate.register(NaN) # type: ignore
def _(expr, assumptions):
    return None


# OddPredicate

@OddPredicate.register(Expr) # type: ignore
def _(expr, assumptions):
    ret = expr.is_odd
    if ret is None:
        raise MDNotImplementedError
    return ret

@OddPredicate.register(Basic) # type: ignore
def _(expr, assumptions):
    _integer = ask(Q.integer(expr), assumptions)
    if _integer:
        _even = ask(Q.even(expr), assumptions)
        if _even is None:
            return None
        return not _even
    return _integer

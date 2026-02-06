from sympy.assumptions.assume import recursive_ask
'\nThis module contains query handlers responsible for calculus queries:\ninfinitesimal, finite, etc.\n'
from sympy.assumptions import Q
from sympy.core import Expr, Add, Mul, Pow, Symbol
from sympy.core.numbers import NegativeInfinity, GoldenRatio, Infinity, Exp1, ComplexInfinity, ImaginaryUnit, NaN, Number, Pi, E, TribonacciConstant
from sympy.functions import cos, exp, log, sign, sin
from sympy.logic.boolalg import conjuncts
from ..predicates.calculus import FinitePredicate, InfinitePredicate, PositiveInfinitePredicate, NegativeInfinitePredicate

@FinitePredicate.register(Symbol)
def _(expr, assumptions):
    """
    Handles Symbol.
    """
    if expr.is_finite is not None:
        return expr.is_finite
    if Q.finite(expr) in conjuncts(assumptions):
        return True
    return None

@FinitePredicate.register(Add)
def _(expr, assumptions):
    """
    Return True if expr is bounded, False if not and None if unknown.

    Truth Table:

    +-------+-----+-----------+-----------+
    |       |     |           |           |
    |       |  B  |     U     |     ?     |
    |       |     |           |           |
    +-------+-----+---+---+---+---+---+---+
    |       |     |   |   |   |   |   |   |
    |       |     |'+'|'-'|'x'|'+'|'-'|'x'|
    |       |     |   |   |   |   |   |   |
    +-------+-----+---+---+---+---+---+---+
    |       |     |           |           |
    |   B   |  B  |     U     |     ?     |
    |       |     |           |           |
    +---+---+-----+---+---+---+---+---+---+
    |   |   |     |   |   |   |   |   |   |
    |   |'+'|     | U | ? | ? | U | ? | ? |
    |   |   |     |   |   |   |   |   |   |
    |   +---+-----+---+---+---+---+---+---+
    |   |   |     |   |   |   |   |   |   |
    | U |'-'|     | ? | U | ? | ? | U | ? |
    |   |   |     |   |   |   |   |   |   |
    |   +---+-----+---+---+---+---+---+---+
    |   |   |     |           |           |
    |   |'x'|     |     ?     |     ?     |
    |   |   |     |           |           |
    +---+---+-----+---+---+---+---+---+---+
    |       |     |           |           |
    |   ?   |     |           |     ?     |
    |       |     |           |           |
    +-------+-----+-----------+---+---+---+

        * 'B' = Bounded

        * 'U' = Unbounded

        * '?' = unknown boundedness

        * '+' = positive sign

        * '-' = negative sign

        * 'x' = sign unknown

        * All Bounded -> True

        * 1 Unbounded and the rest Bounded -> False

        * >1 Unbounded, all with same known sign -> False

        * Any Unknown and unknown sign -> None

        * Else -> None

    When the signs are not the same you can have an undefined
    result as in oo - oo, hence 'bounded' is also undefined.
    """
    sign = -1
    result = True
    for arg in expr.args:
        _bounded = recursive_ask(Q.finite(arg), assumptions=assumptions)
        if _bounded:
            continue
        s = recursive_ask(Q.extended_positive(arg), assumptions=assumptions)
        if sign != -1 and s != sign or (s is None and None in (_bounded, sign)):
            return None
        else:
            sign = s
        if result is not False:
            result = _bounded
    return result

@FinitePredicate.register(Mul)
def _(expr, assumptions):
    """
    Return True if expr is bounded, False if not and None if unknown.

    Truth Table:

    +---+---+---+--------+
    |   |   |   |        |
    |   | B | U |   ?    |
    |   |   |   |        |
    +---+---+---+---+----+
    |   |   |   |   |    |
    |   |   |   | s | /s |
    |   |   |   |   |    |
    +---+---+---+---+----+
    |   |   |   |        |
    | B | B | U |   ?    |
    |   |   |   |        |
    +---+---+---+---+----+
    |   |   |   |   |    |
    | U |   | U | U | ?  |
    |   |   |   |   |    |
    +---+---+---+---+----+
    |   |   |   |        |
    | ? |   |   |   ?    |
    |   |   |   |        |
    +---+---+---+---+----+

        * B = Bounded

        * U = Unbounded

        * ? = unknown boundedness

        * s = signed (hence nonzero)

        * /s = not signed
    """
    result = True
    possible_zero = False
    for arg in expr.args:
        _bounded = recursive_ask(Q.finite(arg), assumptions=assumptions)
        if _bounded:
            if recursive_ask(Q.zero(arg), assumptions=assumptions) is not False:
                if result is False:
                    return None
                possible_zero = True
        elif _bounded is None:
            if result is None:
                return None
            if recursive_ask(Q.extended_nonzero(arg), assumptions=assumptions) is None:
                return None
            if result is not False:
                result = None
        else:
            if possible_zero:
                return None
            result = False
    return result

@FinitePredicate.register(Pow)
def _(expr, assumptions):
    """
    * Unbounded ** NonZero -> Unbounded

    * Bounded ** Bounded -> Bounded

    * Abs()<=1 ** Positive -> Bounded

    * Abs()>=1 ** Negative -> Bounded

    * Otherwise unknown
    """
    if expr.base == E:
        return recursive_ask(Q.finite(expr.exp), assumptions=assumptions)
    base_bounded = recursive_ask(Q.finite(expr.base), assumptions=assumptions)
    exp_bounded = recursive_ask(Q.finite(expr.exp), assumptions=assumptions)
    if base_bounded is None and exp_bounded is None:
        return None
    if base_bounded is False and recursive_ask(Q.extended_nonzero(expr.exp), assumptions=assumptions):
        return False
    if base_bounded and exp_bounded:
        is_base_zero = recursive_ask(Q.zero(expr.base), assumptions=assumptions)
        is_exp_negative = recursive_ask(Q.negative(expr.exp), assumptions=assumptions)
        if is_base_zero is True and is_exp_negative is True:
            return False
        if is_base_zero is not False and is_exp_negative is not False:
            return None
        return True
    if (abs(expr.base) <= 1) == True and recursive_ask(Q.extended_positive(expr.exp), assumptions=assumptions):
        return True
    if (abs(expr.base) >= 1) == True and recursive_ask(Q.extended_negative(expr.exp), assumptions=assumptions):
        return True
    if (abs(expr.base) >= 1) == True and exp_bounded is False:
        return False
    return None

@FinitePredicate.register(exp)
def _(expr, assumptions):
    return recursive_ask(Q.finite(expr.exp), assumptions=assumptions)

@FinitePredicate.register(log)
def _(expr, assumptions):
    if recursive_ask(Q.infinite(expr.args[0]), assumptions=assumptions):
        return False
    return recursive_ask(~Q.zero(expr.args[0]), assumptions=assumptions)

@FinitePredicate.register_many(cos, sin, Number, Pi, Exp1, GoldenRatio, TribonacciConstant, ImaginaryUnit, sign)
def _(expr, assumptions):
    return True

@FinitePredicate.register_many(ComplexInfinity, Infinity, NegativeInfinity)
def _(expr, assumptions):
    return False

@FinitePredicate.register(NaN)
def _(expr, assumptions):
    return None

@InfinitePredicate.register(Expr)
def _(expr, assumptions):
    is_finite = Q.finite(expr)._eval_ask(assumptions)
    if is_finite is None:
        return None
    return not is_finite

@PositiveInfinitePredicate.register(Infinity)
def _(expr, assumptions):
    return True

@PositiveInfinitePredicate.register_many(NegativeInfinity, ComplexInfinity)
def _(expr, assumptions):
    return False

@NegativeInfinitePredicate.register(NegativeInfinity)
def _(expr, assumptions):
    return True

@NegativeInfinitePredicate.register_many(Infinity, ComplexInfinity)
def _(expr, assumptions):
    return False

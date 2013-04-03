"""
AskHandlers related to order relations: positive, negative, etc.
"""
from sympy.assumptions import Q, ask
from sympy.assumptions.handlers import CommonHandler
from sympy.assumptions.assume import AppliedPredicate
from sympy.logic.boolalg import And
from sympy.core.relational import (StrictGreaterThan, StrictLessThan,
                                   _Greater, _Less)

class AskNegativeHandler(CommonHandler):
    """
    This is called by ask() when key='negative'

    Test that an expression is less (strict) than zero.

    Examples:

    >>> from sympy import ask, Q, pi
    >>> ask(Q.negative(pi+1)) # this calls AskNegativeHandler.Add
    False
    >>> ask(Q.negative(pi**2)) # this calls AskNegativeHandler.Pow
    False

    """

    @staticmethod
    def _number(expr, assumptions):
        if not expr.as_real_imag()[1]:
            return expr.evalf() < 0
        else:
            return False

    @staticmethod
    def Basic(expr, assumptions):
        if expr.is_number:
            return AskNegativeHandler._number(expr, assumptions)

    @staticmethod
    def Expr(expr, assumptions):
        return AskPositiveHandler.Expr(-expr, assumptions)

    @staticmethod
    def Add(expr, assumptions):
        """
        Positive + Positive -> Positive,
        Negative + Negative -> Negative
        """
        return AskPositiveHandler.Add(-expr, assumptions)

    @staticmethod
    def Mul(expr, assumptions):
        return ask(Q.positive(-expr), assumptions)

    @staticmethod
    def Pow(expr, assumptions):
        """
        Real ** Even -> NonNegative
        Real ** Odd  -> same_as_base
        NonNegative ** Positive -> NonNegative
        """
        if expr.is_number:
            return AskNegativeHandler._number(expr, assumptions)
        if ask(Q.real(expr.base), assumptions):
            if ask(Q.positive(expr.base), assumptions):
                return False
            if ask(Q.even(expr.exp), assumptions):
                return False
            if ask(Q.odd(expr.exp), assumptions):
                return ask(Q.negative(expr.base), assumptions)

    ImaginaryUnit, Abs = [staticmethod(CommonHandler.AlwaysFalse)]*2


class AskNonZeroHandler(CommonHandler):
    """
    Handler for key 'zero'
    Test that an expression is not identically zero
    """

    @staticmethod
    def Basic(expr, assumptions):
        if expr.is_number:
            # if there are no symbols just evalf
            return expr.evalf() != 0

    @staticmethod
    def Add(expr, assumptions):
        if all(ask(Q.positive(x), assumptions) for x in expr.args) \
                or all(ask(Q.negative(x), assumptions) for x in expr.args):
            return True

    @staticmethod
    def Mul(expr, assumptions):
        for arg in expr.args:
            result = ask(Q.nonzero(arg), assumptions)
            if result:
                continue
            return result
        return True

    @staticmethod
    def Pow(expr, assumptions):
        return ask(Q.nonzero(expr.base), assumptions)

    NaN = staticmethod(CommonHandler.AlwaysTrue)

    @staticmethod
    def Abs(expr, assumptions):
        return ask(Q.nonzero(expr.args[0]), assumptions)


class AskPositiveHandler(CommonHandler):
    """
    Handler for key 'positive'
    Test that an expression is greater (strict) than zero
    """

    @staticmethod
    def _number(expr, assumptions):
        if not expr.as_real_imag()[1]:
            return expr.evalf() > 0
        else:
            return False

    @staticmethod
    def Basic(expr, assumptions):
        if expr.is_number:
            return AskPositiveHandler._number(expr, assumptions)


    @staticmethod
    def Expr(expr, assumptions):
        if expr.is_number:
            return AskPositiveHandler._number(expr, assumptions)
        if assumptions is True:
            return

        # See if assumtion a implies Q.positive(expr)
        # TODO replace expr > xxx
        # by ask(Q.is_true(expr > xxx), assumptions)
        # This must be done carefully to avoid recursive call.
        def transitivity(relation):
            # Prefer not to use the sign function,
            # so that sign function could call ask(Q.positive)
            positive = expr - relation.gts + relation.lts
            if (positive > 0) is True:
                return True
            negative = expr + relation.gts - relation.lts
            if (negative < 0) is True:
                return False
            if isinstance(relation, (StrictGreaterThan, StrictLessThan)):
                if ask(Q.nonzero(positive)) is False:
                    return True
                if ask(Q.nonzero(negative)) is False:
                    return False
        # Iterate assumptions.
        assumptions_list = assumptions.args if assumptions.func is And \
                else [assumptions]
        for a in assumptions_list:
            if a.func is Q.positive:
                rtn = transitivity(a.arg > 0)
            elif a.func is Q.negative:
                rtn = transitivity(a.arg < 0)
            elif a.func is Q.is_true and isinstance(a.arg, (_Greater, _Less)):
                rtn = transitivity(a.arg)
            else:
                rtn = None
            if rtn != None:
                return rtn

        # TODO here we could use solve or _solve_inequality
        # with assumptions that contains expr.
        # However, we must be careful with that. For example:
        # Suppose x**2 - 1 < 0 is it true that x < 1 ?
        # solve(x**2 - 1 < 0) says
        # And(-1 < re(x), im(x) == 0, re(x) < 1)
        # so we think that for all solutions x < 1,
        # but I**2 - 1 < 0  and I is not < 1.
        # It would help if solve could have a mode find_all
        # were it fails if it doesn't find all solutions.


    @staticmethod
    def Mul(expr, assumptions):
        if expr.is_number:
            return AskPositiveHandler._number(expr, assumptions)
        result = True
        for arg in expr.args:
            if ask(Q.positive(arg), assumptions):
                continue
            elif ask(Q.negative(arg), assumptions):
                result = result ^ True
            else:
                return
        return result

    @staticmethod
    def Add(expr, assumptions):
        """
        Positive + Positive -> Positive,
        Negative + Negative -> Negative
        """
        if expr.is_number:
            return AskPositiveHandler._number(expr, assumptions)
        for arg in expr.args:
            if ask(Q.positive(arg), assumptions) is not True:
                break
        else:
            # if all argument's are positive
            return True
        for arg in expr.args:
            if ask(Q.negative(arg), assumptions) is not True:
                break
        else:
            # if all argument's are negative
            return False
        return AskPositiveHandler.Expr(expr, assumptions)

    @staticmethod
    def Pow(expr, assumptions):
        if expr.is_number:
            return expr.evalf() > 0
        if ask(Q.positive(expr.base), assumptions):
            return True
        if ask(Q.negative(expr.base), assumptions):
            if ask(Q.even(expr.exp), assumptions):
                return True
            if ask(Q.even(expr.exp), assumptions):
                return False

    @staticmethod
    def exp(expr, assumptions):
        if ask(Q.real(expr.args[0]), assumptions):
            return True

    ImaginaryUnit = staticmethod(CommonHandler.AlwaysFalse)

    @staticmethod
    def Abs(expr, assumptions):
        return ask(Q.nonzero(expr), assumptions)

"""
AskHandlers related to order relations: positive, negative, etc.
"""
from sympy.assumptions import Q, ask
from sympy.assumptions.handlers import CommonHandler


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
        else: return False

    @staticmethod
    def Basic(expr, assumptions):
        if expr.is_number:
            return AskNegativeHandler._number(expr, assumptions)

    @staticmethod
    def Add(expr, assumptions):
        """
        Positive + Positive -> Positive,
        Negative + Negative -> Negative
        """
        if expr.is_number:
            return AskNegativeHandler._number(expr, assumptions)
        for arg in expr.args:
            if not ask(Q.negative(arg), assumptions):
                break
        else:
            # if all argument's are negative
            return True

    @staticmethod
    def Mul(expr, assumptions):
        if expr.is_number:
            return AskNegativeHandler._number(expr, assumptions)
        result = None
        for arg in expr.args:
            if result is None: result = False
            if ask(Q.negative(arg), assumptions):
                result = not result
            elif ask(Q.positive(arg), assumptions):
                pass
            else: return
        return result

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

    @staticmethod
    def ImaginaryUnit(expr, assumptions):
        return False

    @staticmethod
    def Abs(expr, assumptions):
        return False

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
            if result: continue
            return result
        return True

    @staticmethod
    def Pow(expr, assumptions):
        return ask(Q.nonzero(expr.base), assumptions)

    @staticmethod
    def NaN(expr, assumptions):
        return True

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
        else: return False

    @staticmethod
    def Basic(expr, assumptions):
        if expr.is_number:
            return AskPositiveHandler._number(expr, assumptions)

    @staticmethod
    def Mul(expr, assumptions):
        if expr.is_number:
            return AskPositiveHandler._number(expr, assumptions)
        result = True
        for arg in expr.args:
            if ask(Q.positive(arg), assumptions): continue
            elif ask(Q.negative(arg), assumptions):
                result = result ^ True
            else: return
        return result

    @staticmethod
    def Add(expr, assumptions):
        if expr.is_number:
            return AskPositiveHandler._number(expr, assumptions)
        for arg in expr.args:
            if ask(Q.positive(arg), assumptions) is not True:
                break
        else:
            # if all argument's are positive
            return True

    @staticmethod
    def Pow(expr, assumptions):
        if expr.is_number: return expr.evalf() > 0
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

    @staticmethod
    def ImaginaryUnit(expr, assumptions):
        return False

    @staticmethod
    def Abs(expr, assumptions):
        return ask(Q.nonzero(expr), assumptions)

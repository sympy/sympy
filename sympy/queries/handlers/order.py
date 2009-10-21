"""
AskHandlers related to order relations: positive, negative, etc.
"""
from sympy.utilities import all # python2.4 compatibility
from sympy.queries import Q, ask
from sympy.queries.handlers import CommonHandler


class AskNegativeHandler(CommonHandler):
    """
    This is called by ask() when key='negative'

    Test that an expression is less (strict) than zero.

    Examples:

    >>> from sympy import *
    >>> ask(pi+1, Q.negative) # this calls AskNegativeHandler.Add
    False
    >>> ask(pi**2, Q.negative) # this calls AskNegativeHandler.Pow
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
            if not ask(arg, Q.negative, assumptions):
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
            if ask(arg, Q.negative, assumptions):
                result = not result
            elif ask(arg, Q.positive, assumptions):
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
        if ask(expr.base, Q.real, assumptions):
            if ask(expr.base, Q.positive, assumptions):
                return False
            if ask(expr.exp, Q.even, assumptions):
                return False
            if ask(expr.exp, Q.odd, assumptions):
                return ask(expr.base, Q.negative, assumptions)

    @staticmethod
    def ImaginaryUnit(expr, assumptions):
        return False

    @staticmethod
    def abs(expr, assumptions):
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
        if all([ask(x, Q.positive, assumptions) for x in expr.args]) \
            or all([ask(x, Q.negative, assumptions) for x in expr.args]):
            return True

    @staticmethod
    def Mul(expr, assumptions):
        for arg in expr.args:
            result = ask(arg, Q.nonzero, assumptions)
            if result: continue
            return result
        return True

    @staticmethod
    def Pow(expr, assumptions):
        return ask(expr.base, Q.nonzero, assumptions)

    @staticmethod
    def NaN(expr, assumptions):
        return True

    @staticmethod
    def abs(expr, assumptions):
        return ask(expr.args[0], Q.nonzero, assumptions)

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
            if ask(arg, Q.positive, assumptions): continue
            elif ask(arg, Q.negative, assumptions):
                result = result ^ True
            else: return
        return result

    @staticmethod
    def Add(expr, assumptions):
        if expr.is_number:
            return AskPositiveHandler._number(expr, assumptions)
        for arg in expr.args:
            if ask(arg, Q.positive, assumptions) is not True:
                break
        else:
            # if all argument's are positive
            return True

    @staticmethod
    def Pow(expr, assumptions):
        if expr.is_number: return expr.evalf() > 0
        if ask(expr.base, Q.positive, assumptions):
            return True
        if ask(expr.base, Q.negative, assumptions):
            if ask(expr.exp, Q.even, assumptions):
                return True
            if ask(expr.exp, Q.even, assumptions):
                return False

    @staticmethod
    def exp(expr, assumptions):
        if ask(expr.args[0], Q.real, assumptions):
            return True

    @staticmethod
    def ImaginaryUnit(expr, assumptions):
        return False

    @staticmethod
    def abs(expr, assumptions):
        return ask(expr, Q.nonzero, assumptions)

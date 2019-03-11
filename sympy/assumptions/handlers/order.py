"""
AskHandlers related to order relations: positive, negative, etc.
"""
from __future__ import print_function, division

from sympy.assumptions import Q, ask
from sympy.assumptions.handlers import CommonHandler
from sympy.core.logic import fuzzy_not, fuzzy_and, fuzzy_or


class AskNegativeHandler(CommonHandler):
    """
    This is called by ask() when key='negative'

    Test that an expression is less (strict) than zero.

    Examples
    ========

    >>> from sympy import ask, Q, pi
    >>> ask(Q.negative(pi+1)) # this calls AskNegativeHandler.Add
    False
    >>> ask(Q.negative(pi**2)) # this calls AskNegativeHandler.Pow
    False

    """

    @staticmethod
    def Expr(expr, assumptions):
        return expr.is_negative

    @staticmethod
    def _number(expr, assumptions):
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

        r = ask(Q.real(expr), assumptions)
        if r is not True:
            return r

        # Helper method to check positive or negative
        def neg_pos(expr):
            if ask(Q.positive(expr), assumptions):
                res = 1  # positive
            elif ask(Q.negative(expr), assumptions):
                res = 2  # negative
            else:
                nn = ask(Q.nonnegative(expr), assumptions)
                np = ask(Q.nonpositive(expr), assumptions)
                if nn and np:
                    res = 3  # zero
                elif nn:
                    res = 4  # non-negative
                elif np:
                    res = 5  # non-positive
                elif nn is None and np is None:
                    return None
            return res

        val = 3  # zero
        contras = [
            (1, 2), (2, 1),  # positive + negative := unknown
            (1, 5), (5, 1),  # positive + non-positive := unknown
            (2, 4), (4, 2),  # negative + non-negative := unknown
            (4, 5), (5, 4),  # non-positive + non-negative := unknown
        ]

        for arg in expr.args:
            check = neg_pos(arg)
            if check is None:
                return check
            elif check == val or check == 3 or val == 3:  # added with 0 or same type.
                val = check
            elif (val, check) in contras:
                return None
            elif (val == 1 and check == 4) or (check == 1 and val == 4):
                val = 1  # positive + non-negative := positive
            elif (val == 2 and check == 5) or (check == 2 and val == 5):
                val = 2

        # True if negative, False if positive, zero, non-negative
        # None if non-positive
        result = True if val == 2 else (None if val == 5 else False)
        return result

    @staticmethod
    def Mul(expr, assumptions):
        if expr.is_number:
            return AskNegativeHandler._number(expr, assumptions)
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
                if ask(Q.real(expr.exp), assumptions):
                    return False
            if ask(Q.even(expr.exp), assumptions):
                return False
            if ask(Q.odd(expr.exp), assumptions):
                return ask(Q.negative(expr.base), assumptions)

    ImaginaryUnit, Abs = [staticmethod(CommonHandler.AlwaysFalse)]*2

    @staticmethod
    def exp(expr, assumptions):
        if ask(Q.real(expr.args[0]), assumptions):
            return False


class AskNonNegativeHandler(CommonHandler):

    @staticmethod
    def Expr(expr, assumptions):
        return expr.is_nonnegative

    @staticmethod
    def Basic(expr, assumptions):
        if expr.is_number:
            notnegative = fuzzy_not(AskNegativeHandler._number(expr, assumptions))
            if notnegative:
                return ask(Q.real(expr), assumptions)
            else:
                return notnegative


class AskNonZeroHandler(CommonHandler):
    """
    Handler for key 'zero'
    Test that an expression is not identically zero
    """

    @staticmethod
    def Expr(expr, assumptions):
        return expr.is_nonzero

    @staticmethod
    def Basic(expr, assumptions):
        if ask(Q.real(expr)) is False:
            return False
        if expr.is_number:
            # if there are no symbols just evalf
            i = expr.evalf(2)
            def nonz(i):
                if i._prec != 1:
                    return i != 0
            return fuzzy_or(nonz(i) for i in i.as_real_imag())

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

class AskZeroHandler(CommonHandler):

    @staticmethod
    def Expr(expr, assumptions):
        return expr.is_zero

    @staticmethod
    def Basic(expr, assumptions):
        return fuzzy_and([fuzzy_not(ask(Q.nonzero(expr), assumptions)),
            ask(Q.real(expr), assumptions)])

    @staticmethod
    def Mul(expr, assumptions):
        # TODO: This should be deducible from the nonzero handler
        return fuzzy_or(ask(Q.zero(arg), assumptions) for arg in expr.args)

class AskNonPositiveHandler(CommonHandler):

    @staticmethod
    def Expr(expr, assumptions):
        return expr.is_nonpositive

    @staticmethod
    def Basic(expr, assumptions):
        if expr.is_number:
            notpositive = fuzzy_not(AskPositiveHandler._number(expr, assumptions))
            if notpositive:
                return ask(Q.real(expr), assumptions)
            else:
                return notpositive

class AskPositiveHandler(CommonHandler):
    """
    Handler for key 'positive'
    Test that an expression is greater (strict) than zero
    """

    @staticmethod
    def Expr(expr, assumptions):
        return expr.is_positive

    @staticmethod
    def _number(expr, assumptions):
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
            if ask(Q.positive(arg), assumptions):
                continue
            elif ask(Q.negative(arg), assumptions):
                result = result ^ True
            else:
                return
        return result

    @staticmethod
    def Add(expr, assumptions):
        if expr.is_number:
            return AskPositiveHandler._number(expr, assumptions)

        r = ask(Q.real(expr), assumptions)
        if r is not True:
            return r

        # Helper method to check positive or negative
        def neg_pos(expr):
            if ask(Q.positive(expr), assumptions):
                res = 1  # positive
            elif ask(Q.negative(expr), assumptions):
                res = 2  # negative
            else:
                nn = ask(Q.nonnegative(expr), assumptions)
                np = ask(Q.nonpositive(expr), assumptions)
                if nn and np:
                    res = 3  # zero
                elif nn:
                    res = 4  # non-negative
                elif np:
                    res = 5  # non-positive
                elif nn is None and np is None:
                    return None
            return res

        val = 3  # zero
        contras = [
            (1, 2), (2, 1),  # positive + negative := unknown
            (1, 5), (5, 1),  # positive + non-positive := unknown
            (2, 4), (4, 2),  # negative + non-negative := unknown
            (4, 5), (5, 4),  # non-positive + non-negative := unknown
        ]

        for arg in expr.args:
            check = neg_pos(arg)
            if check is None:
                return check
            elif check == val or check == 3 or val == 3:  # added with 0 or same type.
                val = check
            elif (val, check) in contras:
                return None
            elif (val == 1 and check == 4) or (check == 1 and val == 4):
                val = 1   # positive + non-negative := positive
            elif (val == 2 and check == 5) or (check == 2 and val == 5):
                val = 2

        # True if positive, False if negative, zero, non-positive
        # None if non-negative
        result = True if val == 1 else (None if val == 4 else False)
        return result

    @staticmethod
    def Pow(expr, assumptions):
        if expr.is_number:
            return AskPositiveHandler._number(expr, assumptions)
        if ask(Q.positive(expr.base), assumptions):
            if ask(Q.real(expr.exp), assumptions):
                return True
        if ask(Q.negative(expr.base), assumptions):
            if ask(Q.even(expr.exp), assumptions):
                return True
            if ask(Q.odd(expr.exp), assumptions):
                return False

    @staticmethod
    def exp(expr, assumptions):
        if ask(Q.real(expr.args[0]), assumptions):
            return True
        if ask(Q.imaginary(expr.args[0]), assumptions):
            from sympy import pi, I
            return ask(Q.even(expr.args[0]/(I*pi)), assumptions)

    @staticmethod
    def log(expr, assumptions):
        r = ask(Q.real(expr.args[0]), assumptions)
        if r is not True:
            return r
        if ask(Q.positive(expr.args[0] - 1), assumptions):
            return True
        if ask(Q.negative(expr.args[0] - 1), assumptions):
            return False

    @staticmethod
    def factorial(expr, assumptions):
        x = expr.args[0]
        if ask(Q.integer(x) & Q.positive(x), assumptions):
            return True

    ImaginaryUnit = staticmethod(CommonHandler.AlwaysFalse)

    @staticmethod
    def Abs(expr, assumptions):
        return ask(Q.nonzero(expr), assumptions)

    @staticmethod
    def Trace(expr, assumptions):
        if ask(Q.positive_definite(expr.arg), assumptions):
            return True

    @staticmethod
    def Determinant(expr, assumptions):
        if ask(Q.positive_definite(expr.arg), assumptions):
            return True

    @staticmethod
    def MatrixElement(expr, assumptions):
        if (expr.i == expr.j
                and ask(Q.positive_definite(expr.parent), assumptions)):
            return True

    @staticmethod
    def atan(expr, assumptions):
        return ask(Q.positive(expr.args[0]), assumptions)

    @staticmethod
    def asin(expr, assumptions):
        x = expr.args[0]
        if ask(Q.positive(x) & Q.nonpositive(x - 1), assumptions):
            return True
        if ask(Q.negative(x) & Q.nonnegative(x + 1), assumptions):
            return False

    @staticmethod
    def acos(expr, assumptions):
        x = expr.args[0]
        if ask(Q.nonpositive(x - 1) & Q.nonnegative(x + 1), assumptions):
            return True

    @staticmethod
    def acot(expr, assumptions):
        return ask(Q.real(expr.args[0]), assumptions)

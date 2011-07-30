"""
This module contains query handlers responsible for calculus queries:
infinitesimal, bounded, etc.
"""
from sympy.logic.boolalg import conjuncts
from sympy.assumptions import Q, ask
from sympy.assumptions.handlers import CommonHandler

class AskInfinitesimalHandler(CommonHandler):
    """
    Handler for key 'infinitesimal'
    Test that a given expression is equivalent to an infinitesimal
    number
    """

    @staticmethod
    def _number(expr, assumptions):
        # helper method
        return expr.evalf() == 0

    @staticmethod
    def Basic(expr, assumptions):
        if expr.is_number:
            return AskInfinitesimalHandler._number(expr, assumptions)

    @staticmethod
    def Mul(expr, assumptions):
        """
        Infinitesimal*Bounded -> Infinitesimal
        """
        if expr.is_number:
            return AskInfinitesimalHandler._number(expr, assumptions)
        result = False
        for arg in expr.args:
            if ask(Q.infinitesimal(arg), assumptions):
                result = True
            elif ask(Q.bounded(arg), assumptions):
                continue
            else: break
        else:
            return result

    Add, Pow = Mul, Mul

    @staticmethod
    def Number(expr, assumptions):
        return expr == 0

    NumberSymbol = Number

    @staticmethod
    def ImaginaryUnit(expr, assumptions):
        return False


class AskBoundedHandler(CommonHandler):
    """
    Handler for key 'bounded'.

    Test that an expression is bounded respect to all its variables.

    Example of usage:

    >>> from sympy import Symbol, Q
    >>> from sympy.assumptions.handlers.calculus import AskBoundedHandler
    >>> from sympy.abc import x
    >>> a = AskBoundedHandler()
    >>> a.Symbol(x, Q.positive(x))
    False
    >>> a.Symbol(x, Q.bounded(x))
    True

    """

    @staticmethod
    def Symbol(expr, assumptions):
        """
        Handles Symbol.

        Example:

        >>> from sympy import Symbol, Q
        >>> from sympy.assumptions.handlers.calculus import AskBoundedHandler
        >>> from sympy.abc import x
        >>> a = AskBoundedHandler()
        >>> a.Symbol(x, Q.positive(x))
        False
        >>> a.Symbol(x, Q.bounded(x))
        True

        """
        if Q.bounded(expr) in conjuncts(assumptions):
            return True
        return False

    @staticmethod
    def Add(expr, assumptions):
        """
        Bounded + Bounded     -> Bounded
        Unbounded + Bounded   -> Unbounded
        Unbounded + Unbounded -> ?
        """
        result = True
        for arg in expr.args:
            _bounded = ask(Q.bounded(arg), assumptions)
            if _bounded: continue
            elif _bounded is None: return
            elif _bounded is False:
                if result: result = False
                else: return
        return result

    Mul = Add

    @staticmethod
    def Pow(expr, assumptions):
        """
        Unbounded ** Whatever -> Unbounded
        Bounded ** Unbounded -> Unbounded if base > 1
        Bounded ** Unbounded -> Unbounded if base < 1
        """
        base_bounded = ask(Q.bounded(expr.base), assumptions)
        if not base_bounded:
            return base_bounded
        if ask(Q.bounded(expr.exp), assumptions) and base_bounded:
            return True
        if base_bounded and expr.base.is_number:
            # We need to implement relations for this
            if abs(expr.base) > 1:
                return False
            return True

    @staticmethod
    def log(expr, assumptions):
        return ask(Q.bounded(expr.args[0]), assumptions)

    exp = log

    @staticmethod
    def sin(expr, assumptions):
        return True

    cos = sin

    @staticmethod
    def Number(expr, assumptions):
        return True

    @staticmethod
    def Infinity(expr, assumptions):
        return False

    @staticmethod
    def NegativeInfinity(expr, assumptions):
        return False

    @staticmethod
    def Pi(expr, assumptions):
        return True

    @staticmethod
    def Exp1(expr, assumptions):
        return True

    @staticmethod
    def ImaginaryUnit(expr, assumptions):
        return True

    @staticmethod
    def sign(expr, assumptions):
        return True

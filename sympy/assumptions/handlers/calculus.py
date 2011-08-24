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
    >>> a.Symbol(x, Q.positive(x)) == None
    True
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
        >>> a.Symbol(x, Q.positive(x)) == None
        True
        >>> a.Symbol(x, Q.bounded(x))
        True

        """
        if Q.bounded(expr) in conjuncts(assumptions):
            return True
        return None

    @staticmethod
    def Add(expr, assumptions):
        """
        Return True if expr is bounded, False if not and None if unknown.

               TRUTH TABLE

              B    U     ?
                 + - x + - x
            +---+-----+-----+
        B   | B |  U  |? ? ?|  legend:
            +---+-----+-----+    B  = Bounded
          +     |U ? ?|U ? ?|    U  = Unbounded
        U -     |? U ?|? U ?|    ?  = unknown boundedness
          x     |? ? ?|? ? ?|    +  = positive sign
                +-----+--+--+    -  = negative sign
        ?             |? ? ?|    x  = sign unknown
                      +--+--+


        All Bounded -> True
        Any Unbounded and all same sign -> False
        Any Unknown and unknown sign -> None
        Else -> None

        When the signs are not the same you can have an undefined
        (hence bounded undefined) result as in oo - oo
        """

        result = True
        sign = -1 # not assigned yet
        for arg in expr.args:
            _bounded = ask(Q.bounded(arg), assumptions)
            if _bounded:
                continue
            if result is None and _bounded is None and sign is None:
                return None
            if result is not False:
                result = _bounded
            pos = ask(Q.positive(arg), assumptions)
            if sign == -1:
                sign = pos
                continue
            if sign != pos:
                return None
            if sign is None and pos is None:
                return None
        return result

    @staticmethod
    def Mul(expr, assumptions):
        """
        Return True if expr is bounded, False if not and None if unknown.

               TRUTH TABLE

              B   U     ?
                      s   /s
            +---+---+---+---+
         B  | B | U |   ?   |  legend:
            +---+---+---+---+    B  = Bounded
         U      | U | U | ? |    U  = Unbounded
                +---+---+---+    ?  = unknown boundedness
         ?          |   ?   |    s  = signed (hence nonzero)
                    +---+---+    /s = not signed

        """
        result = True
        for arg in expr.args:
            _bounded = ask(Q.bounded(arg), assumptions)
            if _bounded:
                continue
            elif _bounded is None:
                if result is None:
                    return None
                if ask(Q.nonzero(arg), assumptions) is None:
                    return None
                if result is not False:
                    result = None
            else:
                result = False
        return result

    @staticmethod
    def Pow(expr, assumptions):
        """
        Unbounded ** NonZero -> Unbounded
        Bounded ** Bounded -> Bounded
        Abs()<=1 ** Positive -> Bounded
        Abs()>=1 ** Negative -> Bounded
        Otherwise unknown
        """
        base_bounded = ask(Q.bounded(expr.base), assumptions)
        exp_bounded = ask(Q.bounded(expr.exp), assumptions)
        if base_bounded==None and exp_bounded==None: # Common Case
            return None
        if base_bounded==False and ask(Q.nonzero(expr.exp), assumptions):
            return False
        if base_bounded and exp_bounded:
            return True
        if abs(expr.base)<=1 and ask(Q.positive(expr.exp), assumptions):
            return True
        if abs(expr.base)>=1 and ask(Q.negative(expr.exp), assumptions):
            return True
        if abs(expr.base)>=1 and exp_bounded==False:
            return False
        return None

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

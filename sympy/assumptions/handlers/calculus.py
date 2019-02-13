"""
This module contains query handlers responsible for calculus queries:
infinitesimal, finite, etc.
"""
from __future__ import print_function, division

from sympy.logic.boolalg import conjuncts
from sympy.assumptions import Q, ask
from sympy.assumptions.handlers import CommonHandler


class AskFiniteHandler(CommonHandler):
    """
    Handler for key 'finite'.

    Test that an expression is bounded respect to all its variables.

    Examples of usage:

    >>> from sympy import Symbol, Q
    >>> from sympy.assumptions.handlers.calculus import AskFiniteHandler
    >>> from sympy.abc import x
    >>> a = AskFiniteHandler()
    >>> a.Symbol(x, Q.positive(x)) == None
    True
    >>> a.Symbol(x, Q.finite(x))
    True

    """

    @staticmethod
    def Symbol(expr, assumptions):
        """
        Handles Symbol.

        Examples
        ========

        >>> from sympy import Symbol, Q
        >>> from sympy.assumptions.handlers.calculus import AskFiniteHandler
        >>> from sympy.abc import x
        >>> a = AskFiniteHandler()
        >>> a.Symbol(x, Q.positive(x)) == None
        True
        >>> a.Symbol(x, Q.finite(x))
        True

        """
        if expr.is_finite is not None:
            return expr.is_finite
        if Q.finite(expr) in conjuncts(assumptions):
            return True
        return None

    @staticmethod
    def Add(expr, assumptions):
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

|

            * All Bounded -> True

            * 1 Unbounded and the rest Bounded -> False

            * >1 Unbounded, all with same known sign -> False

            * Any Unknown and unknown sign -> None

            * Else -> None

        When the signs are not the same you can have an undefined
        result as in oo - oo, hence 'bounded' is also undefined.

        """

        sign = -1  # sign of unknown or infinite
        result = True
        for arg in expr.args:
            _bounded = ask(Q.finite(arg), assumptions)
            if _bounded:
                continue
            s = ask(Q.positive(arg), assumptions)
            # if there has been more than one sign or if the sign of this arg
            # is None and Bounded is None or there was already
            # an unknown sign, return None
            if sign != -1 and s != sign or \
                    s is None and (s == _bounded or s == sign):
                return None
            else:
                sign = s
            # once False, do not change
            if result is not False:
                result = _bounded
        return result

    @staticmethod
    def Mul(expr, assumptions):
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
        for arg in expr.args:
            _bounded = ask(Q.finite(arg), assumptions)
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
        base_bounded = ask(Q.finite(expr.base), assumptions)
        exp_bounded = ask(Q.finite(expr.exp), assumptions)
        if base_bounded is None and exp_bounded is None:  # Common Case
            return None
        if base_bounded is False and ask(Q.nonzero(expr.exp), assumptions):
            return False
        if base_bounded and exp_bounded:
            return True
        if (abs(expr.base) <= 1) == True and ask(Q.positive(expr.exp), assumptions):
            return True
        if (abs(expr.base) >= 1) == True and ask(Q.negative(expr.exp), assumptions):
            return True
        if (abs(expr.base) >= 1) == True and exp_bounded is False:
            return False
        return None

    @staticmethod
    def log(expr, assumptions):
        return ask(Q.finite(expr.args[0]), assumptions)

    exp = log

    cos, sin, Number, Pi, Exp1, GoldenRatio, TribonacciConstant, ImaginaryUnit, sign = \
        [staticmethod(CommonHandler.AlwaysTrue)]*9

    Infinity, NegativeInfinity = [staticmethod(CommonHandler.AlwaysFalse)]*2

class AskTranscendentalHandler(CommonHandler):
    """
    Handler for key 'transcendental'.

    Test that an expression is transcendental.

    Examples of usage:

    >>> from sympy import Symbol, Q, E, Pi
    >>> from sympy.assumptions.handlers.calculus import AskTranscendentalHandler
    >>> from sympy.abc import x
    >>> a = AskTranscendentalHandler()
    >>> a.Add(E + Pi)
    None
    >>> a.Mul(2*E)
    True

    """

    @staticmethod
    def Symbol(expr, assumptions):
        """
        Handles Symbol.

        Examples
        ========

        >>> from sympy import Symbol, Q
        >>> from sympy.assumptions.handlers.calculus import AskFiniteHandler
        >>> from sympy.abc import x
        >>> a = AskTranscendentalHandler()
        >>> a.Symbol(x, Q.transcendental(x))
        True
        >>> a.Symbol(x)
        None

        """

        if expr.is_transcendental is not None:
            return expr.is_transcendental
        if Q.transcendental(expr) in conjuncts(assumptions):
            return True
        return None

    @staticmethod
    def Add(expr, assumptions):
        """
        Returns True if expr can be proven transcendental.

        If only one of the arguments is transcendental then
        `True` is returned.
        If all of the arguments are non-transcendental then
        `False` is returned.
        In all other cases, `None` is returned.

        """

        _num_transcendental = 0
        for arg in expr.args:
            if ask(Q.transcendental(arg), assumptions):
                _num_transcendental += 1
        if _num_transcendental == 0:
            return False
        if _num_transcendental == 1:
            return True
        return None

    @staticmethod
    def Mul(expr, assumptions):
        """
        Returns True if expr can be proven transcendental.

        If only one of the arguments is transcendental then
        True is returned.
        If all of the arguments are non-transcendental then
        False is returned.
        In all other cases, None is returned.

        """

        expr = sum([arg for arg in expr.args])
        return Add(expr, assumptions)

    @staticmethod
    def Pow(expr, assumptions):
        """
        Returns True if expr can be proven transcendental.

        If base is transcendental and exponent is
        non-transcendental then True is returned.
        If base is non-transcendental and exponent is
        irrational then also True is returned.
        In all other cases, None is returned.

        """

        base_transcendental = ask(Q.transcendental(expr.base), assumptions)
        exp_transcendental = ask(Q.transcendental(exp.exp), assumptions)
        exp_irrational = ask(Q.rational(exp.exp), assumptions)
        if base_transcendental and not exp_transcendental:
            return True
        if not base_transcendental and exp_irrational:
            return True
        return None

    @staticmethod
    def log(expr, assumptions):
        if expr.args[0] not in (0, 1):
            return not ask(Q.transcendental(expr.args[0]), assumptions)
        return None

    @staticmethod
    def cos(expr, assumptions):
        arg_transcendental = ask(Q.transcendental(expr.args[0]), assumptions)
        if not arg_transcendental and expr.args[0] != 0:
            return True
        return None

    @staticmethod
    def sin(expr, assumptions):
        return cos(expr, assumptions)

    Pi, Exp1 = [staticmethod(CommonHandler.AlwaysTrue)]*2
    GoldenRatio, TribonacciConstant, ImaginaryUnit = [staticmethod(CommonHandler.AlwaysFalse)]*3

    Infinity, NegativeInfinity = [staticmethod(CommonHandler.AlwaysNone)]*2

"""
Handlers for keys related to number theory: prime, even, odd, etc.
"""
from sympy.assumptions import Q, ask
from sympy.assumptions.handlers import CommonHandler
from sympy.ntheory import isprime

class AskPrimeHandler(CommonHandler):
    """
    Handler for key 'prime'
    Test that an expression represents a prime number. When the
    expression is a number the result, when True, is subject to
    the limitations of isprime() which is used to return the result.
    """

    @staticmethod
    def _number(expr, assumptions):
        # helper method
        try:
            i = int(expr.round())
            if not (expr - i).equals(0):
                raise TypeError
        except TypeError:
            return False
        return isprime(i)

    @staticmethod
    def Basic(expr, assumptions):
        # Just use int(expr) once
        # http://code.google.com/p/sympy/issues/detail?id=1462
        # is solved
        if expr.is_number:
            return AskPrimeHandler._number(expr, assumptions)

    @staticmethod
    def Mul(expr, assumptions):
        if expr.is_number:
            return AskPrimeHandler._number(expr, assumptions)
        for arg in expr.args:
            if ask(Q.integer(arg), assumptions):
                pass
            else: break
        else:
            # a product of integers can't be a prime
            return False

    @staticmethod
    def Pow(expr, assumptions):
        """
        Integer**Integer     -> !Prime
        """
        if expr.is_number:
            return AskPrimeHandler._number(expr, assumptions)
        if ask(Q.integer(expr.exp), assumptions) and \
                ask(Q.integer(expr.base), assumptions):
            return False

    @staticmethod
    def Integer(expr, assumptions):
        return isprime(expr)

    @staticmethod
    def Rational(expr, assumptions):
        return False

    @staticmethod
    def Float(expr, assumptions):
        return AskPrimeHandler._number(expr, assumptions)

    @staticmethod
    def Infinity(expr, assumptions):
        return False

    @staticmethod
    def NegativeInfinity(expr, assumptions):
        return False

    @staticmethod
    def ImaginaryUnit(expr, assumptions):
        return False

    @staticmethod
    def NumberSymbol(expr, assumptions):
        return AskPrimeHandler._number(expr, assumptions)

class AskCompositeHandler(CommonHandler):

    @staticmethod
    def Basic(expr, assumptions):
        _positive = ask(Q.positive(expr), assumptions)
        if _positive:
            _integer = ask(Q.integer(expr), assumptions)
            if _integer:
                _prime = ask(Q.prime(expr), assumptions)
                if _prime is None: return
                return not _prime
            else: return _integer
        else: return _positive

class AskEvenHandler(CommonHandler):

    @staticmethod
    def _number(expr, assumptions):
        # helper method
        try:
            i = int(expr.round())
            if not (expr - i).equals(0):
                raise TypeError
        except TypeError:
            return False
        return i % 2 == 0

    @staticmethod
    def Basic(expr, assumptions):
        if expr.is_number:
            return AskEvenHandler._number(expr, assumptions)

    @staticmethod
    def Mul(expr, assumptions):
        """
        Even * Integer -> Even
        Even * Odd     -> Even
        Integer * Odd  -> ?
        Odd * Odd      -> Odd
        """
        if expr.is_number:
            return AskEvenHandler._number(expr, assumptions)
        even, odd, irrational = False, 0, False
        for arg in expr.args:
            # check for all integers and at least one even
            if ask(Q.integer(arg), assumptions):
                if ask(Q.even(arg), assumptions):
                    even = True
                elif ask(Q.odd(arg), assumptions):
                    odd += 1
            elif ask(Q.irrational(arg), assumptions):
                # one irrational makes the result False
                # two makes it undefined
                if irrational:
                    break
                irrational = True
            else: break
        else:
            if irrational: return False
            if even: return True
            if odd == len(expr.args): return False

    @staticmethod
    def Add(expr, assumptions):
        """
        Even + Odd  -> Odd
        Even + Even -> Even
        Odd  + Odd  -> Even

        """
        if expr.is_number:
            return AskEvenHandler._number(expr, assumptions)
        _result = True
        for arg in expr.args:
            if ask(Q.even(arg), assumptions):
                pass
            elif ask(Q.odd(arg), assumptions):
                _result = not _result
            else: break
        else:
            return _result

    @staticmethod
    def Integer(expr, assumptions):
        return not bool(expr.p & 1)

    @staticmethod
    def Rational(expr, assumptions):
        return False

    @staticmethod
    def Float(expr, assumptions):
        return expr % 2 == 0

    @staticmethod
    def Infinity(expr, assumptions):
        return False

    @staticmethod
    def NegativeInfinity(expr, assumptions):
        return False

    @staticmethod
    def NumberSymbol(expr, assumptions):
        return AskEvenHandler._number(expr, assumptions)

    @staticmethod
    def ImaginaryUnit(expr, assumptions):
        return False

    @staticmethod
    def Abs(expr, assumptions):
        if ask(Q.real(expr.args[0]), assumptions):
            return ask(Q.even(expr.args[0]), assumptions)

    @staticmethod
    def re(expr, assumptions):
        if ask(Q.real(expr.args[0]), assumptions):
            return ask(Q.even(expr.args[0]), assumptions)

    @staticmethod
    def im(expr, assumptions):
        if ask(Q.real(expr.args[0]), assumptions):
            return True

class AskOddHandler(CommonHandler):
    """
    Handler for key 'odd'
    Test that an expression represents an odd number
    """

    @staticmethod
    def Basic(expr, assumptions):
        _integer = ask(Q.integer(expr), assumptions)
        if _integer:
            _even = ask(Q.even(expr), assumptions)
            if _even is None: return None
            return not _even
        return _integer

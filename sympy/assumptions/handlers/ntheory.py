"""
Handlers for keys related to number theory: prime, even, odd, etc.
"""
from sympy.assumptions import Q, ask
from sympy.assumptions.handlers import CommonHandler
from sympy.ntheory import isprime

class AskPrimeHandler(CommonHandler):
    """
    Handler for key 'prime'
    Test that an expression represents a prime number
    """

    @staticmethod
    def _number(expr, assumptions):
        # helper method
        if (expr.as_real_imag()[1] == 0) and int(expr.evalf()) == expr:
            return isprime(expr.evalf(1))
        return False

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
            if ask(arg, Q.integer, assumptions):
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
        if ask(expr.exp, Q.integer, assumptions) and \
                ask(expr.base, Q.integer, assumptions):
            return False

    @staticmethod
    def Integer(expr, assumptions):
        return isprime(expr)

    @staticmethod
    def Rational(expr, assumptions):
        return False

    @staticmethod
    def Real(expr, assumptions):
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
        _positive = ask(expr, Q.positive, assumptions)
        if _positive:
            _integer = ask(expr, Q.integer, assumptions)
            if _integer:
                _prime = ask(expr, Q.prime, assumptions)
                if _prime is None: return
                return not _prime
            else: return _integer
        else: return _positive

class AskEvenHandler(CommonHandler):

    @staticmethod
    def _number(expr, assumptions):
        # helper method
        if (expr.as_real_imag()[1] == 0) and expr.evalf(1) == expr:
            return float(expr.evalf()) % 2 == 0
        else: return False

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
            if ask(arg, Q.integer, assumptions):
                if ask(arg, Q.even, assumptions):
                    even = True
                elif ask(arg, Q.odd, assumptions):
                    odd += 1
            elif ask(arg, Q.irrational, assumptions):
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

        TODO: remove float() when issue
        http://code.google.com/p/sympy/issues/detail?id=1473
        is solved
        """
        if expr.is_number:
            return AskEvenHandler._number(expr, assumptions)
        _result = True
        for arg in expr.args:
            if ask(arg, Q.even, assumptions):
                pass
            elif ask(arg, Q.odd, assumptions):
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
    def Real(expr, assumptions):
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
    def abs(expr, assumptions):
        if ask(expr.args[0], Q.real, assumptions):
            return ask(expr.args[0], Q.even, assumptions)

    @staticmethod
    def re(expr, assumptions):
        if ask(expr.args[0], Q.real, assumptions):
            return ask(expr.args[0], Q.even, assumptions)

    @staticmethod
    def im(expr, assumptions):
        if ask(expr.args[0], Q.real, assumptions):
            return True

class AskOddHandler(CommonHandler):
    """
    Handler for key 'odd'
    Test that an expression represents an odd number
    """

    @staticmethod
    def Basic(expr, assumptions):
        _integer = ask(expr, Q.integer, assumptions)
        if _integer:
            _even = ask(expr, Q.even, assumptions)
            if _even is None: return None
            return not _even
        return _integer


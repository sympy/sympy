"""
Handlers for keys related to set membership: integer, rational, etc.
"""
from sympy.assumptions import Q, ask
from sympy.assumptions.handlers import CommonHandler

class AskIntegerHandler(CommonHandler):
    """
    Handler for key 'integer'
    Test that an expression belongs to the field of integer numbers
    """

    @staticmethod
    def _number(expr, assumptions):
        # helper method
        if expr.as_real_imag()[1] == 0:
            return expr.evalf(1) == expr
        return False

    @staticmethod
    def Add(expr, assumptions):
        """
        Integer + Integer       -> Integer
        Integer + !Integer      -> !Integer
        !Integer + !Integer -> ?
        """
        if expr.is_number:
            return AskIntegerHandler._number(expr, assumptions)
        return test_closed_group(expr, assumptions, 'integer')

    @staticmethod
    def Mul(expr, assumptions):
        """
        Integer*Integer      -> Integer
        Integer*Irrational   -> !Integer
        Odd/Even             -> !Integer
        Integer*Rational     -> ?
        """
        if expr.is_number:
            return AskIntegerHandler._number(expr, assumptions)
        _output = True
        for arg in expr.args:
            if not ask(arg, Q.integer, assumptions):
                if arg.is_Rational:
                    if arg.q == 2:
                        return ask(2*expr, Q.even, assumptions)
                    if ~(arg.q & 1):
                        return None
                elif ask(arg, Q.irrational, assumptions):
                    if _output:
                        _output = False
                    else:
                        return
                else:
                    return
        else:
            return _output

    Pow = Add

    @staticmethod
    def int(expr, assumptions):
        return True

    @staticmethod
    def Integer(expr, assumptions):
        return True

    @staticmethod
    def Rational(expr, assumptions):
        # rationals with denominator one get
        # evaluated to Integers
        return False

    @staticmethod
    def Real(expr, assumptions):
        return int(expr) == expr

    @staticmethod
    def Pi(expr, assumptions):
        return False

    @staticmethod
    def Exp1(expr, assumptions):
        return False

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
    def abs(expr, assumptions):
        return ask(expr.args[0], Q.integer, assumptions)

class AskRationalHandler(CommonHandler):
    """
    Handler for key 'rational'
    Test that an expression belongs to the field of rational numbers
    """

    @staticmethod
    def Add(expr, assumptions):
        """
        Rational + Rational     -> Rational
        Rational + !Rational    -> !Rational
        !Rational + !Rational   -> ?
        """
        if expr.is_number:
            if expr.as_real_imag()[1]:
                return False
        return test_closed_group(expr, assumptions, 'rational')

    Mul = Add

    @staticmethod
    def Pow(expr, assumptions):
        """
        Rational ** Integer      -> Rational
        Irrational ** Rational   -> Irrational
        Rational ** Irrational   -> ?
        """
        if ask(expr.exp, Q.integer, assumptions):
            return ask(expr.base, Q.rational, assumptions)
        elif ask(expr.exp, Q.rational, assumptions):
            if ask(expr.base, Q.prime, assumptions):
                return False

    @staticmethod
    def Rational(expr, assumptions):
        return True

    @staticmethod
    def Real(expr, assumptions):
        # it's finite-precission
        return True

    @staticmethod
    def ImaginaryUnit(expr, assumptions):
        return False

    @staticmethod
    def Infinity(expr, assumptions):
        return False

    @staticmethod
    def NegativeInfinity(expr, assumptions):
        return False

    @staticmethod
    def Pi(expr, assumptions):
        return False

    @staticmethod
    def Exp1(expr, assumptions):
        return False

class AskIrrationalHandler(CommonHandler):

    @staticmethod
    def Basic(expr, assumptions):
        _real = ask(expr, Q.real, assumptions)
        if _real:
            _rational = ask(expr, Q.rational, assumptions)
            if _rational is None: return None
            return not _rational
        else: return _real

class AskRealHandler(CommonHandler):
    """
    Handler for key 'real'
    Test that an expression belongs to the field of real numbers
    """

    @staticmethod
    def _number(expr, assumptions):
        return not expr.as_real_imag()[1]

    @staticmethod
    def Add(expr, assumptions):
        """
        Real + Real              -> Real
        Real + (Complex & !Real) -> !Real
        """
        if expr.is_number:
            return AskRealHandler._number(expr, assumptions)
        return test_closed_group(expr, assumptions, 'real')

    @staticmethod
    def Mul(expr, assumptions):
        """
        Real*Real               -> Real
        Real*Imaginary          -> !Real
        Imaginary*Imaginary     -> Real
        """
        if expr.is_number:
            return AskRealHandler._number(expr, assumptions)
        result = True
        for arg in expr.args:
            if ask(arg, Q.real, assumptions):
                pass
            elif ask(arg, Q.imaginary, assumptions):
                result = result ^ True
            else:
                break
        else:
            return result

    @staticmethod
    def Pow(expr, assumptions):
        """
        Real**Integer         -> Real
        Positive**Real        -> Real
        Real**(Integer/Even)  -> Real if base is nonnegative
        Real**(Integer/Odd)   -> Real
        """
        if expr.is_number:
            return AskRealHandler._number(expr, assumptions)
        if ask(expr.base, Q.real, assumptions):
            if ask(expr.exp, Q.integer, assumptions):
                return True
            elif expr.exp.is_Rational:
                if (expr.exp.q % 2 == 0):
                    return ask(expr.base, Q.real, assumptions) and \
                       not ask(expr.base, Q.negative, assumptions)
                else: return True
            elif ask(expr.exp, Q.real, assumptions):
                if ask(expr.base, Q.positive, assumptions):
                    return True

    @staticmethod
    def Rational(expr, assumptions):
        return True

    @staticmethod
    def Real(expr, assumptions):
        return True

    @staticmethod
    def Pi(expr, assumptions):
        return True

    @staticmethod
    def Exp1(expr, assumptions):
        return True

    @staticmethod
    def abs(expr, assumptions):
        return True

    @staticmethod
    def re(expr, assumptions):
        return True

    im = re

    @staticmethod
    def ImaginaryUnit(expr, assumptions):
        return False

    @staticmethod
    def Infinity(expr, assumptions):
        return False

    @staticmethod
    def NegativeInfinity(expr, assumptions):
        return False

    @staticmethod
    def sin(expr, assumptions):
        if ask(expr.args[0], Q.real, assumptions):
            return True

    cos, exp = sin, sin

class AskExtendedRealHandler(AskRealHandler):
    """
    Handler for key 'extended_real'
    Test that an expression belongs to the field of extended real numbers,
    that is real numbers union {Infinity, -Infinity}
    """

    @staticmethod
    def Add(expr, assumptions):
        return test_closed_group(expr, assumptions, 'extended_real')

    Mul, Pow = Add, Add

    @staticmethod
    def Infinity(expr, assumptions):
        return True

    @staticmethod
    def NegativeInfinity(expr, assumptions):
        return True

class AskComplexHandler(CommonHandler):
    """
    Handler for key 'complex'
    Test that an expression belongs to the field of complex numbers
    """

    @staticmethod
    def Add(expr, assumptions):
        return test_closed_group(expr, assumptions, 'complex')

    Mul, Pow = Add, Add

    @staticmethod
    def Number(expr, assumptions):
        return True

    @staticmethod
    def NumberSymbol(expr, assumptions):
        return True

    @staticmethod
    def abs(expr, assumptions):
        return True

    @staticmethod
    def ImaginaryUnit(expr, assumptions):
        return True

    @staticmethod
    def Infinity(expr, assumptions):
        return False

    @staticmethod
    def NegativeInfinity(expr, assumptions):
        return False

    sin, cos, exp, re, im = [abs]*5 # they are all complex functions

class AskImaginaryHandler(CommonHandler):
    """
    Handler for key 'imaginary'
    Test that an expression belongs to the field of imaginary numbers,
    that is, numbers in the form x*I, where x is real
    """

    @staticmethod
    def _number(expr, assumptions):
        # helper method
        return not expr.as_real_imag()[0]

    @staticmethod
    def Add(expr, assumptions):
        """
        Imaginary + Imaginary -> Imaginary
        Imaginary + Complex   -> ?
        Imaginary + Real      -> !Imaginary
        """
        if expr.is_number:
            return AskImaginaryHandler._number(expr, assumptions)
        reals = 0
        for arg in expr.args:
            if ask(arg, Q.imaginary, assumptions):
                pass
            elif ask(arg, Q.real, assumptions):
                reals += 1
            else:
                break
        else:
            if reals == 0:
                return True
            if reals == 1 or (len(expr.args) == reals):
                # two reals could sum 0 thus giving an imaginary
                return False

    @staticmethod
    def Mul(expr, assumptions):
        """
        Real*Imaginary      -> Imaginary
        Imaginary*Imaginary -> Real
        """
        if expr.is_number:
            return AskImaginaryHandler._number(expr, assumptions)
        result = False
        reals = 0
        for arg in expr.args:
            if ask(arg, Q.imaginary, assumptions):
                result = result ^ True
            elif not ask(arg, Q.real, assumptions):
                break
        else:
            if reals == len(expr.args):
                return False
            return result

    Pow = Add

    @staticmethod
    def Number(expr, assumptions):
        return not (expr.as_real_imag()[1] == 0)

    NumberSymbol = Number

    @staticmethod
    def ImaginaryUnit(expr, assumptions):
        return True

class AskAlgebraicHandler(CommonHandler):
    """Handler for 'algebraic' key. """

    @staticmethod
    def Add(expr, assumptions):
        return test_closed_group(expr, assumptions, 'algebraic')

    @staticmethod
    def Mul(expr, assumptions):
        return test_closed_group(expr, assumptions, 'algebraic')

    @staticmethod
    def Pow(expr, assumptions):
        return expr.exp.is_Rational and ask(expr.base, 'algebraic', assumptions)

    @staticmethod
    def Number(expr, assumptions):
        return False

    @staticmethod
    def Rational(expr, assumptions):
        return expr.q != 0

    @staticmethod
    def ImaginaryUnit(expr, assumptions):
        return True

    @staticmethod
    def AlgebraicNumber(expr, assumptions):
        return True

#### Helper methods

def test_closed_group(expr, assumptions, key):
    """
    Test for membership in a group with respect
    to the current operation
    """
    result = True
    for arg in expr.args:
        _out = ask(arg, key, assumptions)
        if _out is None: break
        elif _out is False:
            if result: result = False
            else: break
    else:
        return result


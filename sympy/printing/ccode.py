"""
C code printer
"""

from str import StrPrinter
from sympy.printing.precedence import precedence, PRECEDENCE
from sympy.core.basic import S

class CCodePrinter(StrPrinter):
    """A printer to convert python expressions to stings of c code"""

    def _print_Pow(self, expr):
        PREC = precedence(expr)
        if expr.exp is S.NegativeOne:
            return '1/%s'%(self.parenthesize(expr.base, PREC))
        else:
            return 'pow(%s,%s)'%(self.parenthesize(expr.base, PREC),
                                 self.parenthesize(expr.exp, PREC))

    def _print_Exp1(self, expr):
        return "exp(1)"

def ccode(expr):
    r"""Converts an expr to a string of c code

        Works for simple expressions using math.h functions.

        >>> from sympy import *
        >>> from sympy.abc import *

        >>> ccode((2*tau)**Rational(7,2))
        '8*pow(2,(1/2))*pow(tau,(7/2))'
    """
    return CCodePrinter().doprint(expr)

def print_ccode(expr):
    """Prints C representation of the given expression."""
    print ccode(expr)

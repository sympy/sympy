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

    def _print_Piecewise(self, expr):
        ecpairs = ["(%s) {\n%s\n}\n" % (self._print(c), self._print(e)) \
                       for e, c in expr.args[:-1]]
        last_line = ""
        if expr.args[-1].cond is S.One:
            last_line = "else {\n%s\n}" % self._print(expr.args[-1].expr)
        else:
            ecpairs.append("(%s) {\n%s\n" % \
                           (self._print(expr.args[-1].cond),
                            self._print(expr.args[-1].expr)))
        code = "if %s" + last_line
        return code % "else if ".join(ecpairs)


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

"""
C code printer

The CCodePrinter converts single sympy expressions into single C expressions,
using the functions defined in math.h where possible.

A complete code generator, which uses ccode extensively, can be found in
sympy.utilities.codegen. The codegen module can be used to generate complete
source code files that are compilable without further modifications.
"""

from str import StrPrinter
from sympy.printing.precedence import precedence, PRECEDENCE
from sympy.core.basic import S

class CCodePrinter(StrPrinter):
    """A printer to convert python expressions to strings of c code"""
    printmethod = "_ccode_"

    def _print_Pow(self, expr):
        PREC = precedence(expr)
        if expr.exp is S.NegativeOne:
            return '1.0/%s'%(self.parenthesize(expr.base, PREC))
        else:
            return 'pow(%s,%s)'%(self.parenthesize(expr.base, PREC),
                                 self.parenthesize(expr.exp, PREC))

    def _print_Rational(self, expr):
        p, q = int(expr.p), int(expr.q)
        return '%d.0/%d.0' % (p, q)

    def _print_Exp1(self, expr):
        return "M_E"

    def _print_Pi(self, expr):
        return 'M_PI'

    def _print_Infinity(self, expr):
        return 'HUGE_VAL'

    def _print_NegativeInfinity(self, expr):
        return '-HUGE_VAL'

    def _print_Piecewise(self, expr):
        ecpairs = ["(%s) {\n%s\n}\n" % (self._print(c), self._print(e)) \
                       for e, c in expr.args[:-1]]
        last_line = ""
        if expr.args[-1].cond == True:
            last_line = "else {\n%s\n}" % self._print(expr.args[-1].expr)
        else:
            ecpairs.append("(%s) {\n%s\n" % \
                           (self._print(expr.args[-1].cond),
                            self._print(expr.args[-1].expr)))
        code = "if %s" + last_line
        return code % "else if ".join(ecpairs)

    def _print_Function(self, expr):
        if expr.func.__name__ == "ceiling":
            return "ceil(%s)" % self.stringify(expr.args, ", ")
        elif expr.func.__name__ == "abs" and not expr.args[0].is_integer:
            return "fabs(%s)" % self.stringify(expr.args, ", ")
        else:
            return StrPrinter._print_Function(self, expr)


def ccode(expr):
    r"""Converts an expr to a string of c code

        Works for simple expressions using math.h functions.

        >>> from sympy import *
        >>> from sympy.abc import *

        >>> ccode((2*tau)**Rational(7,2))
        '8*pow(2,(1.0/2.0))*pow(tau,(7.0/2.0))'

    """
    return CCodePrinter().doprint(expr)

def print_ccode(expr):
    """Prints C representation of the given expression."""
    print ccode(expr)

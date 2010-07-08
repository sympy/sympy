"""
C code printer

The CCodePrinter converts single sympy expressions into single C expressions,
using the functions defined in math.h where possible.

A complete code generator, which uses ccode extensively, can be found in
sympy.utilities.codegen. The codegen module can be used to generate complete
source code files that are compilable without further modifications.
"""

from str import StrPrinter
from sympy.printing.precedence import precedence
from sympy.core.basic import S


# dictionary mapping sympy function to (argument_conditions, C_function).
# Used in CCodePrinter._print_Function(self)
known_functions = {
        "ceiling": [(lambda x: True, "ceil")],
        "abs": [(lambda x: not x.is_integer, "fabs")],
        }

class CCodePrinter(StrPrinter):
    """A printer to convert python expressions to strings of c code"""
    printmethod = "_ccode"

    _default_settings = {
        'order': None,
        'full_prec': 'auto',
        'precision': 15,
        'user_functions': {},
        'human': True,
    }

    def __init__(self, settings={}):
        """Register function mappings supplied by user"""
        StrPrinter.__init__(self, settings)
        self.known_functions = dict(known_functions)
        userfuncs = settings.get('user_functions', {})
        for k,v in userfuncs.items():
            if not isinstance(v, tuple):
                userfuncs[k] = (lambda *x: True, v)
        self.known_functions.update(userfuncs)

    def _print_Pow(self, expr):
        PREC = precedence(expr)
        if expr.exp is S.NegativeOne:
            return '1.0/%s'%(self.parenthesize(expr.base, PREC))
        elif expr.exp == 0.5:
            return 'sqrt(%s)' % self._print(expr.base)
        else:
            return 'pow(%s, %s)'%(self._print(expr.base),
                                 self._print(expr.exp))

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

    def _print_And(self, expr):
        PREC = precedence(expr)
        return '&&'.join(self.parenthesize(a, PREC) for a in expr.args)

    def _print_Or(self, expr):
        PREC = precedence(expr)
        return '||'.join(self.parenthesize(a, PREC) for a in expr.args)

    def _print_Not(self, expr):
        PREC = precedence(expr)
        return '!'+self.parenthesize(expr.args[0], PREC)

    def _print_Function(self, expr):
        if expr.func.__name__ in self.known_functions:
            cond_cfunc = self.known_functions[expr.func.__name__]
            for cond, cfunc in cond_cfunc:
                if cond(*expr.args):
                    return "%s(%s)" % (cfunc, self.stringify(expr.args, ", "))
        return StrPrinter._print_Function(self, expr)


def ccode(expr, **settings):
    r"""Converts an expr to a string of c code

        Arguments:
          expr  --  a sympy expression to be converted

        Optional arguments:
          precision  --  the precision for numbers such as pi [default=15]
          user_functions  --  A dictionary where keys are FunctionClass instances
                              and values are there string representations.
                              Alternatively, the dictionary value can be a list
                              of tuples i.e. [(argument_test, cfunction_string)].
                              See below for examples.
          human  --  If True, the result is a single string that may contain
                     some constant declarations for the number symbols. If
                     False, the same information is returned in a more
                     programmer-friendly data structure.

        >>> from sympy import ccode, symbols, Rational, sin
        >>> x, tau = symbols(["x", "tau"])
        >>> ccode((2*tau)**Rational(7,2))
        '8*sqrt(2)*pow(tau, 7.0/2.0)'
        >>> ccode(sin(x), assign_to="s")
        's = sin(x);'


    """
    return CCodePrinter(settings).doprint(expr)

def print_ccode(expr, **settings):
    """Prints C representation of the given expression."""
    print ccode(expr, **settings)

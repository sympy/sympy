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
from sympy.core.numbers import NumberSymbol
from sympy.functions import Piecewise, piecewise_fold
from sympy.tensor import Idx


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

    def doprint(self, expr, assign_to=None):

        # keep a set of expressions that are not strictly translatable to C
        # and number constants that must be declared and initialized
        not_c = self._not_c = set()
        self._number_symbols = set()

        # We treat Piecewise here to ensure it is top-level and outside loops
        expr = piecewise_fold(expr)
        lines = []
        if isinstance(expr, Piecewise):
            for i, (e, c) in enumerate(expr.args):
                if i == 0:
                    lines.append("if (%s) {" % self._print(c))
                elif i == len(expr.args)-1 and c == True:
                    lines.append("else {")
                else:
                    lines.append("else if (%s) {" % self._print(c))
                code0 = self._doprint_a_piece(e, assign_to)
                lines.extend(code0)
                lines.append("}")
        else:
            code0 = self._doprint_a_piece(expr, assign_to)
            lines.extend(code0)

        # format the output
        if self._settings["human"]:
            frontlines = []
            if len(not_c) > 0:
                frontlines.append("// Not C:")
                for expr in sorted(not_c, key=str):
                    frontlines.append("// %s" % expr)
            for name, value in sorted(self._number_symbols, key=str):
                frontlines.append("double const %s = %s" % (name, value))
            lines = frontlines + lines
            # lines = self.indent_code(lines)
            result = "\n".join(lines)
        else:
            # lines = self.indent_code(lines)
            result = self._number_symbols, not_c, "\n".join(lines)
        del self._not_c
        del self._number_symbols
        return result

    def _doprint_a_piece(self, expr, assign_to=None):

        # Setup loops if expression contain Indexed objects
        openloop, closeloop, local_ints = self._get_loop_opening_ending_ints(expr)

        # the lhs may contain loops that are not in the rhs
        lhs = assign_to
        if lhs:
            open_lhs, close_lhs, lhs_ints = self._get_loop_opening_ending_ints(lhs)
            for n,ind in enumerate(lhs_ints):
                if ind not in local_ints:
                    openloop.insert(0, open_lhs[n])
                    closeloop.append(close_lhs[n])
            lhs_printed = self._print(lhs)

        lines = openloop
        line = StrPrinter.doprint(self, expr)
        if assign_to is None:
            text = "%s" % line
        else:
            text = "%s = %s" % (lhs_printed, line)
        lines.append(text)
        lines.extend(closeloop)
        return lines

    def _get_loop_opening_ending_ints(self, expr):
        """Returns a tuple (open_lines, close_lines) containing lists of codelines
        """
        indices = expr.atoms(Idx)
        # FIXME: sort indices in an optimized way
        open_lines = []
        close_lines = []
        local_ints = []

        loopstart = "for (int %(var)s = %(start)s; %(var)s < %(end)s; %(var)s++){"
        for i in indices:
            # C arrays start at 0 and end at dimension-1
            open_lines.append(loopstart % {
                'var': i.label,
                'start': i.lower,
                'end': i.upper})
            close_lines.append("}")
            local_ints.append(i)
        return open_lines, close_lines, local_ints

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

    def _print_IndexedElement(self, expr):
        # calculate index for 1d array
        dims = expr.dimensions
        inds = [ i.label for i in expr.indices ]
        elem = S.Zero
        offset = S.One
        for i in reversed(range(expr.rank)):
            elem += offset*inds[i]
            offset *= dims[i]
        return "%s[%s]" % (self._print(expr.stem), self._print(elem))

    def _print_Exp1(self, expr):
        return "M_E"

    def _print_Pi(self, expr):
        return 'M_PI'

    def _print_Infinity(self, expr):
        return 'HUGE_VAL'

    def _print_NegativeInfinity(self, expr):
        return '-HUGE_VAL'

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


def ccode(expr, assign_to=None, **settings):
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
    return CCodePrinter(settings).doprint(expr, assign_to)

def print_ccode(expr, **settings):
    """Prints C representation of the given expression."""
    print ccode(expr, **settings)

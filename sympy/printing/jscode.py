"""
Javascript code printer

The JavascriptCodePrinter converts single sympy expressions into single
Javascript expressions, using the functions defined in the Javascript
Math object where possible.

"""

from __future__ import print_function, division

from sympy.core import S, C
from sympy.printing.codeprinter import CodePrinter, Assignment
from sympy.printing.precedence import precedence
from sympy.core.compatibility import string_types


# dictionary mapping sympy function to (argument_conditions, Javascript_function).
# Used in JavascriptCodePrinter._print_Function(self)
known_functions = {
}

function_translations = {
    'Abs': 'Math.abs',
    'acos': 'Math.acos',
    'asin': 'Math.asin',
    'atan': 'Math.atan',
    'ceiling': 'Math.ceil',
    'cos': 'Math.cos',
    'exp': 'Math.exp',
    'floor': 'Math.floor',
    'log': 'Math.log',
    'sin': 'Math.sin',
    'tan': 'Math.tan',
}


class JavascriptCodePrinter(CodePrinter):
    """"A Printer to convert python expressions to strings of javascript code
    """
    printmethod = '_javascript'
    language = 'Javascript'

    _default_settings = {
        'order': None,
        'full_prec': 'auto',
        'precision': 15,
        'user_functions': {},
        'human': True,
    }

    def __init__(self, settings={}):
        CodePrinter.__init__(self, settings)
        self.known_functions = dict(known_functions)
        userfuncs = settings.get('user_functions', {})
        for k, v in userfuncs.items():
            if not isinstance(v, tuple):
                userfuncs[k] = (lambda *x: True, v)
        self.known_functions.update(userfuncs)

    def _rate_index_position(self, p):
        return p*5

    def _get_statement(self, codestring):
        return "%s;" % codestring

    def _get_comment(self, text):
        return "// {:}".format(text)

    def _declare_number_const(self, name, value):
        return "var {:} = {:};".format(name, value)

    def _format_code(self, lines):
        return self.indent_code(lines)

    def _traverse_matrix_indices(self, mat):
        rows, cols = mat.shape
        return ((i, j) for i in range(rows) for j in range(cols))

    def _get_loop_opening_ending(self, indices):
        open_lines = []
        close_lines = []
        loopstart = "for (var %(varble)s=%(start)s; %(varble)s<%(end)s; %(varble)s++){"
        for i in indices:
            # Javascript arrays start at 0 and end at dimension-1
            open_lines.append(loopstart % {
                'varble': self._print(i.label),
                'start': self._print(i.lower),
                'end': self._print(i.upper + 1)})
            close_lines.append("}")
        return open_lines, close_lines

    def _print_Pow(self, expr):
        PREC = precedence(expr)
        if expr.exp == -1:
            return '1/%s' % (self.parenthesize(expr.base, PREC))
        elif expr.exp == 0.5:
            return 'Math.sqrt(%s)' % self._print(expr.base)
        else:
            return 'Math.pow(%s, %s)' % (self._print(expr.base),
                                 self._print(expr.exp))

    def _print_Rational(self, expr):
        p, q = int(expr.p), int(expr.q)
        return '%d/%d' % (p, q)

    def _print_Indexed(self, expr):
        # calculate index for 1d array
        dims = expr.shape
        inds = [ i.label for i in expr.indices ]
        elem = S.Zero
        offset = S.One
        for i in reversed(range(expr.rank)):
            elem += offset*inds[i]
            offset *= dims[i]
        return "%s[%s]" % (self._print(expr.base.label), self._print(elem))

    def _print_Exp1(self, expr):
        return "Math.E"

    def _print_Pi(self, expr):
        return 'Math.PI'

    def _print_Infinity(self, expr):
        return 'Number.POSITIVE_INFINITY'

    def _print_NegativeInfinity(self, expr):
        return 'Number.NEGATIVE_INFINITY'

    def _print_Piecewise(self, expr):
        lines = []
        if expr.has(Assignment):
            for i, (e, c) in enumerate(expr.args):
                if i == 0:
                    lines.append("if (%s) {" % self._print(c))
                elif i == len(expr.args) - 1 and c == True:
                    lines.append("else {")
                else:
                    lines.append("else if (%s) {" % self._print(c))
                code0 = self._print(e)
                lines.append(code0)
                lines.append("}")
            return "\n".join(lines)
        else:
            # The piecewise was used in an expression, need to do inline
            # operators. This has the downside that if none of the conditions
            # are true, the last expression will still be returned. Also, these
            # inline operators will not work for statements that span multiple
            # lines (Matrix or Indexed expressions).
            ecpairs = ["((%s) ? (\n%s\n)\n" % (self._print(c), self._print(e))
                    for e, c in expr.args[:-1]]
            last_line = ""
            if expr.args[-1].cond == True:
                last_line = ": (\n%s\n)" % self._print(expr.args[-1].expr)
            else:
                ecpairs.append("(%s) ? (\n%s\n)" %
                (self._print(expr.args[-1].cond),
                    self._print(expr.args[-1].expr)))
            code = "%s" + last_line
            return code % ": ".join(ecpairs) + " ".join([")"*len(ecpairs)])

    def _print_Function(self, expr):
        if expr.func.__name__ in self.known_functions:
            cond_cfunc = self.known_functions[expr.func.__name__]
            for cond, cfunc in cond_cfunc:
                if cond(*expr.args):
                    return "%s(%s)" % (cfunc, self.stringify(expr.args, ", "))
        if expr.func.__name__ in function_translations:
            tr = function_translations[expr.func.__name__]
            return "%s(%s)" % (tr, self.stringify(expr.args, ", "))
        if hasattr(expr, '_imp_') and isinstance(expr._imp_, C.Lambda):
            # inlined function
            return self._print(expr._imp_(*expr.args))
        return CodePrinter._print_Function(self, expr)

    def _print_MatrixElement(self, expr):
        return "{:}[{:}][{:}]".format(expr.parent, expr.i, expr.j)

    def indent_code(self, code):
        """Accepts a string of code or a list of code lines"""

        if isinstance(code, string_types):
            code_lines = self.indent_code(code.splitlines(True))
            return ''.join(code_lines)

        tab = "   "
        inc_token = ('{', '(', '{\n', '(\n')
        dec_token = ('}', ')')

        code = [ line.lstrip(' \t') for line in code ]

        increase = [ int(any(map(line.endswith, inc_token))) for line in code ]
        decrease = [ int(any(map(line.startswith, dec_token)))
                     for line in code ]

        pretty = []
        level = 0
        for n, line in enumerate(code):
            if line == '' or line == '\n':
                pretty.append(line)
                continue
            level -= decrease[n]
            pretty.append("%s%s" % (tab*level, line))
            level += increase[n]
        return pretty


def jscode(expr, assign_to=None, **settings):
    """Converts an expr to a string of javascript code

       Parameters
       ==========

       expr : sympy.core.Expr
           a sympy expression to be converted
       precision : optional
           the precision for numbers such as pi [default=15]
       user_functions : optional
           A dictionary where keys are FunctionClass instances and values
           are their string representations. Alternatively the
           dictionary values can be a list of tuples i.e. [(argument_test,
           jsfunction_string)].
       human : optional
           If True, the result is a single string that may contain some
           constant declarations for the number symbols. If False, the
           same information is returned in a more programmer-friendly
           data structure.

       Examples
       ========

       >>> from sympy import jscode, symbols, Rational, sin
       >>> x, tau = symbols(["x", "tau"])
       >>> jscode((2*tau)**Rational(7,2))
       '8*Math.sqrt(2)*Math.pow(tau, 7/2)'
       >>> jscode(sin(x), assign_to="s")
       's = Math.sin(x);'

    """
    return JavascriptCodePrinter(settings).doprint(expr, assign_to)


def print_jscode(expr, **settings):
    """Prints the Javascript representation of the given expression.

       See jscode for the meaning of the optional arguments.
    """
    print(jscode(expr, **settings))

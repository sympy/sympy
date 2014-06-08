"""
C code printer

The CCodePrinter converts single sympy expressions into single C expressions,
using the functions defined in math.h where possible.

A complete code generator, which uses ccode extensively, can be found in
sympy.utilities.codegen. The codegen module can be used to generate complete
source code files that are compilable without further modifications.


"""

from __future__ import print_function, division

from sympy.core import S, C
from sympy.core.compatibility import string_types
from sympy.printing.codeprinter import CodePrinter
from sympy.printing.precedence import precedence

# dictionary mapping sympy function to (argument_conditions, C_function).
# Used in CCodePrinter._print_Function(self)
known_functions = {
    "ceiling": [(lambda x: True, "ceil")],
    "Abs": [(lambda x: not x.is_integer, "fabs")],
    "gamma": [(lambda x: True, "tgamma")],
}


class CCodePrinter(CodePrinter):
    """A printer to convert python expressions to strings of c code"""
    printmethod = "_ccode"

    _default_settings = {
        'order': None,
        'full_prec': 'auto',
        'precision': 15,
        'user_functions': {},
        'human': True,
        'contract': True,
    }

    def __init__(self, settings={}):
        """Register function mappings supplied by user"""
        CodePrinter.__init__(self, settings)
        self.known_functions = dict(known_functions)
        userfuncs = settings.get('user_functions', {})
        for k, v in userfuncs.items():
            if not isinstance(v, list):
                userfuncs[k] = [(lambda *x: True, v)]
        self.known_functions.update(userfuncs)

    def _rate_index_position(self, p):
        """function to calculate score based on position among indices

        This method is used to sort loops in an optimized order, see
        CodePrinter._sort_optimized()
        """
        return p*5

    def _get_statement(self, codestring):
        return "%s;" % codestring

    def doprint(self, expr, assign_to=None):
        """
        Actually format the expression as C code.
        """

        if isinstance(assign_to, string_types):
            assign_to = C.Symbol(assign_to)
        elif not isinstance(assign_to, (C.Basic, type(None))):
            raise TypeError("CCodePrinter cannot assign to object of type %s" %
                    type(assign_to))

        # keep a set of expressions that are not strictly translatable to C
        # and number constants that must be declared and initialized
        not_c = self._not_supported = set()
        self._number_symbols = set()

        # We treat top level Piecewise here to get if tests outside loops
        lines = []
        if isinstance(expr, C.Piecewise):
            for i, (e, c) in enumerate(expr.args):
                if i == 0:
                    lines.append("if (%s) {" % self._print(c))
                elif i == len(expr.args) - 1 and c == True:
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
                    frontlines.append("// %s" % repr(expr))
            for name, value in sorted(self._number_symbols, key=str):
                frontlines.append("double const %s = %s;" % (name, value))
            lines = frontlines + lines
            lines = "\n".join(lines)
            result = self.indent_code(lines)
        else:
            lines = self.indent_code("\n".join(lines))
            result = self._number_symbols, not_c, lines
        del self._not_supported
        del self._number_symbols
        return result

    def _get_loop_opening_ending(self, indices):
        """Returns a tuple (open_lines, close_lines) containing lists of codelines
        """
        open_lines = []
        close_lines = []
        loopstart = "for (int %(var)s=%(start)s; %(var)s<%(end)s; %(var)s++){"
        for i in indices:
            # C arrays start at 0 and end at dimension-1
            open_lines.append(loopstart % {
                'var': self._print(i.label),
                'start': self._print(i.lower),
                'end': self._print(i.upper + 1)})
            close_lines.append("}")
        return open_lines, close_lines

    def _print_Pow(self, expr):
        if "Pow" in self.known_functions:
            return self._print_Function(expr)
        PREC = precedence(expr)
        if expr.exp == -1:
            return '1.0/%s' % (self.parenthesize(expr.base, PREC))
        elif expr.exp == 0.5:
            return 'sqrt(%s)' % self._print(expr.base)
        else:
            return 'pow(%s, %s)' % (self._print(expr.base),
                                 self._print(expr.exp))

    def _print_Rational(self, expr):
        p, q = int(expr.p), int(expr.q)
        return '%d.0L/%d.0L' % (p, q)

    def _print_Indexed(self, expr):
        # calculate index for 1d array
        dims = expr.shape
        elem = S.Zero
        offset = S.One
        for i in reversed(range(expr.rank)):
            elem += expr.indices[i]*offset
            offset *= dims[i]
        return "%s[%s]" % (self._print(expr.base.label), self._print(elem))

    def _print_Idx(self, expr):
        return self._print(expr.label)

    def _print_Exp1(self, expr):
        return "M_E"

    def _print_Pi(self, expr):
        return 'M_PI'

    def _print_Infinity(self, expr):
        return 'HUGE_VAL'

    def _print_NegativeInfinity(self, expr):
        return '-HUGE_VAL'

    def _print_Piecewise(self, expr):
        # This method is called only for inline if constructs
        # Top level piecewise is handled in doprint()
        ecpairs = ["((%s) ? (\n%s\n)\n" % (self._print(c), self._print(e))
                   for e, c in expr.args[:-1]]
        last_line = ""
        if expr.args[-1].cond == True:
            last_line = ": (\n%s\n)" % self._print(expr.args[-1].expr)
        else:
            ecpairs.append("(%s) ? (\n%s\n" %
                           (self._print(expr.args[-1].cond),
                            self._print(expr.args[-1].expr)))
        code = "%s" + last_line
        return code % ": ".join(ecpairs) + " )"

    def _print_Function(self, expr):
        if expr.func.__name__ in self.known_functions:
            cond_cfunc = self.known_functions[expr.func.__name__]
            for cond, cfunc in cond_cfunc:
                if cond(*expr.args):
                    return "%s(%s)" % (cfunc, self.stringify(expr.args, ", "))
        if hasattr(expr, '_imp_') and isinstance(expr._imp_, C.Lambda):
            # inlined function
            return self._print(expr._imp_(*expr.args))
        return CodePrinter._print_Function(self, expr)

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


def ccode(expr, assign_to=None, **settings):
    r"""Converts an expr to a string of c code

        Parameters
        ==========

        expr : sympy.core.Expr
            a sympy expression to be converted
        precision : optional
            the precision for numbers such as pi [default=15]
        user_functions : optional
            A dictionary where keys are FunctionClass instances and values
            are their string representations.  Alternatively, the
            dictionary value can be a list of tuples i.e. [(argument_test,
            cfunction_string)].  See below for examples.
        human : optional
            If True, the result is a single string that may contain some
            constant declarations for the number symbols. If False, the
            same information is returned in a more programmer-friendly
            data structure.
        contract: optional
            If True, `Indexed` instances are assumed to obey
            tensor contraction rules and the corresponding nested
            loops over indices are generated. Setting contract = False
            will not generate loops, instead the user is responsible
            to provide values for the indices in the code. [default=True]


        Examples
        ========

        >>> from sympy import ccode, symbols, Rational, sin, ceiling, Abs
        >>> x, tau = symbols(["x", "tau"])
        >>> ccode((2*tau)**Rational(7,2))
        '8*sqrt(2)*pow(tau, 7.0L/2.0L)'
        >>> ccode(sin(x), assign_to="s")
        's = sin(x);'
        >>> custom_functions = {
        ...   "ceiling": "CEIL",
        ...   "Abs": [(lambda x: not x.is_integer, "fabs"),
        ...           (lambda x: x.is_integer, "ABS")]
        ... }
        >>> ccode(Abs(x) + ceiling(x), user_functions=custom_functions)
        'fabs(x) + CEIL(x)'
        >>> from sympy import Eq, IndexedBase, Idx
        >>> len_y = 5
        >>> y = IndexedBase('y', shape=(len_y,))
        >>> t = IndexedBase('t', shape=(len_y,))
        >>> Dy = IndexedBase('Dy', shape=(len_y-1,))
        >>> i = Idx('i', len_y-1)
        >>> e=Eq(Dy[i], (y[i+1]-y[i])/(t[i+1]-t[i]))
        >>> ccode(e.rhs, assign_to=e.lhs, contract=False)
        'Dy[i] = (y[i + 1] - y[i])/(t[i + 1] - t[i]);'

    """
    return CCodePrinter(settings).doprint(expr, assign_to)


def print_ccode(expr, **settings):
    """Prints C representation of the given expression."""
    print(ccode(expr, **settings))

"""
Octave/Matlab code printer
"""

from __future__ import print_function, division
from sympy.core import C, Add, Mul, Pow, S, Rational
from sympy.core.compatibility import string_types
from sympy.core.mul import _keep_coeff
from sympy.printing.codeprinter import CodePrinter, Assignment
from sympy.printing.str import StrPrinter
from sympy.printing.precedence import precedence

# dictionary mapping sympy function to (argument_conditions, C_function).
# Used in OctaveCodePrinter._print_Function(self)
# FIXME: almost certainly incomplete, perhaps better to default yet!
known_functions = {
    "sin": "sin",
    "cos": "cos",
    "tan": "tan",
    "asin": "asin",
    "acos": "acos",
    "atan": "atan",
    "atan2": "atan2",
    "sinh": "sinh",
    "cosh": "cosh",
    "tanh": "tanh",
    "asinh": "asinh",
    "acosh": "acosh",
    "atanh": "atanh",
    "log": "log",
    "exp": "exp",
    "erf": "erf",
    "gamma": "gamma",
    "sign": "sign",
    "floor": "floor",
    "Abs": "abs",       # different names after this point
    "ceiling": "ceil",
    "conjugate": "conj",
}
# FIXME: test and/or merge these in
# csc
# sec
# cot
# coth
# acot
# acoth
# factorial
# erfc
# erfinv
# erfcinv
# erfi
# heaviside|Heaviside
# dirac|DiracDelta


class OctaveCodePrinter(CodePrinter):
    """
    A printer to convert expressions to strings of Octave/Matlab code.
    """
    printmethod = "_octave"
    language = "Octave"

    _operators = {
        'and': '&&',
        'or': '||',
        'not': '~',
    }

    _default_settings = {
        'order': None,
        'full_prec': 'auto',
        'precision': 15,
        'user_functions': {},
        'human': True,
        'contract': False
    }
    # FIXME: contract is for expressing tensors as loops (if True), or
    # just assignment (if False).  Needs tests for tensors, borrow
    # some from C?

    def __init__(self, settings={}):
        CodePrinter.__init__(self, settings)
        self.known_functions = dict(known_functions)
        userfuncs = settings.get('user_functions', {})
        self.known_functions.update(userfuncs)

    def _rate_index_position(self, p):
        return p*5

    def _get_statement(self, codestring):
        return "%s;" % codestring

    def _get_comment(self, text):
        return "% {0}".format(text)

    def _declare_number_const(self, name, value):
        return "{0} = {1};".format(name, value)

    def _format_code(self, lines):
        return self.indent_code(lines)

    def _traverse_matrix_indices(self, mat):
        # Octave uses Fortran order (column-major)
        # FIXME: need tests for matrix traverse stuff
        rows, cols = mat.shape
        return ((i, j) for j in range(cols) for i in range(rows))

    def _get_loop_opening_ending(self, indices):
        open_lines = []
        close_lines = []
        for i in indices:
            # Octave arrays start at 1 and end at dimension
            var, start, stop = map(self._print,
                    [i.label, i.lower + 1, i.upper + 1])
            open_lines.append("for %s = %s:%s" % (var, start, stop))
            close_lines.append("end")
        return open_lines, close_lines


    def _print_Mul(self, expr):
        # print complex numbers nicely in Octave
        if (expr.is_number and expr.is_imaginary and
                expr.as_coeff_Mul()[0].is_integer):
            return "%si" % self._print(-S.ImaginaryUnit*expr)

        # cribbed from str.py
        prec = precedence(expr)

        c, e = expr.as_coeff_Mul()
        if c < 0:
            expr = _keep_coeff(-c, e)
            sign = "-"
        else:
            sign = ""

        a = []  # items in the numerator
        b = []  # items that are in the denominator (if any)

        if self.order not in ('old', 'none'):
            args = expr.as_ordered_factors()
        else:
            # use make_args in case expr was something like -x -> x
            args = Mul.make_args(expr)

        # Gather args for numerator/denominator
        for item in args:
            if (item.is_commutative and item.is_Pow and item.exp.is_Rational
                    and item.exp.is_negative):
                if item.exp != -1:
                    b.append(Pow(item.base, -item.exp, evaluate=False))
                else:
                    b.append(Pow(item.base, -item.exp))
            elif item.is_Rational and item is not S.Infinity:
                if item.p != 1:
                    a.append(Rational(item.p))
                if item.q != 1:
                    b.append(Rational(item.q))
            else:
                a.append(item)

        a = a or [S.One]

        a_str = list(map(lambda x: self.parenthesize(x, prec), a))
        b_str = list(map(lambda x: self.parenthesize(x, prec), b))

        # from here it differs from str.py to deal with "*" and ".*"
        def multjoin(a, a_str):
            # here we probably are assuming the constants will come first
            r = a_str[0]
            #for (ai, ai_str) in zip(a, a_str)[1:]:
            #    if ai.is_constant()
            for i in range(1, len(a)):
                if a[i-1].is_constant():
                    mulsym = '*'
                else:
                    mulsym = '.*'
                r = r + mulsym + a_str[i]
            return r

        if len(b) == 0:
            return sign + multjoin(a, a_str)
        elif len(b) == 1:
            if b[0].is_constant():
                divsym = '/'
            else:
                divsym = './'
            return sign + multjoin(a, a_str) + divsym + b_str[0]
        else:
            if all([bi.is_constant() for bi in b]):
                divsym = '/'
            else:
                divsym = './'
            return (sign + multjoin(a, a_str) +
                    divsym + "(%s)" % multjoin(b, b_str))

    # FIXME: need tests that use HadamardProduct: user needs this to
    # mix matrix symbols and .* products?
    #def _print_MatMul(self, expr):
    #def _print_HadamardProduct(self, expr):


    def _print_Pow(self, expr):
        return super(OctaveCodePrinter, self)._print_Pow(expr, powsymbol='.^')

    def _print_MatPow(self, expr):
        return super(OctaveCodePrinter, self)._print_MatPow(expr, powsymbol='^')

    def _print_Pi(self, expr):
        return 'pi'

    def _print_ImaginaryUnit(self, expr):
        return "1i"

    def _print_Exp1(self, expr):
        return "exp(1)"

    def _print_GoldenRatio(self, expr):
        # FIXME: how to do better, e.g., for octave_code(2*GoldenRatio)?
        return "(1+sqrt(5))/2"

    def _print_NumberSymbol(self, expr):
        # FIXME: perhaps this should be based on Assignment?  Like in
        # Piecewise.  This form is better for inline code, but the
        # CodePrinter implementation is nicer for longer programs.
        return "%.18g" % float(expr)

    def _print_Infinity(self, expr):
        return 'inf'

    def _print_NegativeInfinity(self, expr):
        return '-inf'

    def _print_NaN(self, expr):
        return 'NaN'

    def _print_list(self, expr):
        return '{' + ', '.join(self._print(a) for a in expr) + '}'
    _print_tuple = _print_list
    _print_Tuple = _print_list

    def _print_MatrixBase(self, A):
        # Handle zero dimensions:
        if (A.rows, A.cols) == (0, 0):
            return '[]'
        elif A.rows == 0 or A.cols == 0:
            return 'zeros(%s, %s)' % (A.rows, A.cols)
        elif (A.rows, A.cols) == (1, 1):
            # Octave does not distinguish between scalars and 1x1 matrices
            return self._print(A[0, 0])
        elif A.rows == 1 :
            return "[%s]" % A.table(self, rowstart='', rowend='', colsep=' ')
        elif A.cols == 1 :
            return "[%s]" % A.table(self, rowstart='', rowend='',
                                    rowsep='; ', colsep=' ')
        return "[%s]" % A.table(self, rowstart='', rowend='',
                                rowsep='; ...\n', colsep=' ')
    # FIXME: see my prosposed change for _print_NumberSymbol, same here
    _print_SparseMatrix = \
        _print_MutableSparseMatrix = \
        _print_ImmutableSparseMatrix = \
        _print_Matrix = \
        _print_DenseMatrix = \
        _print_MutableDenseMatrix = \
        _print_ImmutableMatrix = \
        _print_ImmutableDenseMatrix = \
        _print_MatrixBase

    def _print_Identity(self, expr):
        return "eye(%s)" % self._print(expr.shape[0])

    def _print_Piecewise(self, expr):
        if expr.args[-1].cond != True:
            # We need the last conditional to be a True, otherwise the resulting
            # function may not return a result.
            raise ValueError("All Piecewise expressions must contain an "
                             "(expr, True) statement to be used as a default "
                             "condition. Without one, the generated "
                             "expression may not evaluate to anything under "
                             "some condition.")
        lines = []
        # FIXME: for Octave, user might want to force the inline mode (its
        # better for vector operations for example).  How to expose this?
        # Choosing just based on Assignment seems too coarse.
        if expr.has(Assignment):
            for i, (e, c) in enumerate(expr.args):
                if i == 0:
                    lines.append("if (%s)" % self._print(c))
                elif i == len(expr.args) - 1 and c == True:
                    lines.append("else")
                else:
                    lines.append("elseif (%s)" % self._print(c))
                code0 = self._print(e)
                lines.append(code0)
                if i == len(expr.args) - 1:
                    lines.append("end")
            return "\n".join(lines)
        else:
            # This Piecewise was used in an expression, so do inline
            # where each cond, expr pair is like a nested Horner form:
            #   (condition) .* (expr) + (not cond) .* (<others>)
            # FIXME: ccode says some things won't work inline, true here?
            ecpairs = ["({0}).*({1}) + (~({0})).*(".format
                       (self._print(c), self._print(e))
                       for e, c in expr.args[:-1]]
            elast = "%s" % self._print(expr.args[-1].expr)
            return " ...\n".join(ecpairs) + elast + ")"*len(ecpairs)
            # FIXME: need brackets for 2*pw, see XFAIL test

    def indent_code(self, code):
        """Accepts a string of code or a list of code lines"""

        # code mostly copied from ccode
        if isinstance(code, string_types):
            code_lines = self.indent_code(code.splitlines(True))
            return ''.join(code_lines)

        tab = "  "
        inc_token = ('if ', 'function ', 'else', 'elseif ')
        dec_token = ('end')

        code = [ line.lstrip(' \t') for line in code ]

        increase = [ int(any(map(line.startswith, inc_token)))
                     for line in code ]
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


def octave_code(expr, assign_to=None, **settings):
    r"""Converts an expr to a string of Octave/Matlab code

    Examples
    ========

    FIXME, DUDE WRITE SOME DOCS ;-)
    >>> from sympy import octave_code as mcode, symbols, sin
    >>> x = symbols('x')
    >>> mcode(sin(x).series(x).removeO())
    '(1/120)*x^5 - 1/6*x^3 + x'
    """
    return OctaveCodePrinter(settings).doprint(expr, assign_to)

def print_octave_code(expr, **settings):
    """Prints the Octave/Matlab representation of the given expression.

       See `octave_code` for the meaning of the optional arguments.
    """
    print(octave_code(expr, **settings))

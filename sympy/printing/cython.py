from sympy.core import S
from sympy.core.compatibility import string_types, range
from sympy.printing.codeprinter import CodePrinter

from sympy.core import Add, Mul, Pow, S, Eq
from sympy.core.basic import Basic
from sympy.core.relational import Relational
from sympy.core.sympify import _sympify, sympify
from sympy.core.symbol import Symbol
from sympy.printing.str import StrPrinter

from sympy.printing.precedence import precedence
from sympy.sets.fancysets import Range

# dictionary mapping sympy function to (argument_conditions, C_function).
# Used in CCodePrinter._print_Function(self)
known_functions = {
    "Abs": [(lambda x: not x.is_integer, "fabs")],
    "gamma": "tgamma",
    "sin": "sin",
    "cos": "cos",
    "tan": "tan",
    "asin": "asin",
    "acos": "acos",
    "atan": "atan",
    "atan2": "atan2",
    "exp": "exp",
    "log": "log",
    "erf": "erf",
    "sinh": "sinh",
    "cosh": "cosh",
    "tanh": "tanh",
    "asinh": "asinh",
    "acosh": "acosh",
    "atanh": "atanh",
    "floor": "floor",
    "ceiling": "ceil",
}

reserved_words = ['auto',
                  'if',
                  'break',
                  'int',
                  'case',
                  'long',
                  'char',
                  'register',
                  'continue',
                  'return',
                  'default',
                  'short',
                  'do',
                  'sizeof',
                  'double',
                  'static',
                  'else',
                  'struct',
                  'entry',
                  'switch',
                  'extern',
                  'typedef',
                  'float',
                  'union',
                  'for',
                  'unsigned',
                  'goto',
                  'while',
                  'enum',
                  'void',
                  'const',
                  'signed',
                  'volatile']

class CythonCodePrinter(CodePrinter):
    """A printer to convert python expressions to strings of cython code"""
    printmethod = "_cythoncode"
    language = "Cython"

    _default_settings = {
        'order': None,
        'full_prec': 'auto',
        'precision': 16,
        'user_functions': {},
        'human': True,
        'contract': True,
        'inline': True,
        'dereference': set(),
        'error_on_reserved': False,
        'reserved_word_suffix': '_',
    }

    def __init__(self, settings={}):
        super(CythonCodePrinter, self).__init__(settings)
        self.known_functions = dict(known_functions)
        userfuncs = settings.get('user_functions', {})
        self.known_functions.update(userfuncs)
        self._dereference = set(settings.get('dereference', []))
        self.reserved_words = set(reserved_words)

    def doprint(self, expr, assign_to=None):
        """
        Print the expression as code.

        Parameters
        ----------
        expr : Expression
            The expression to be printed.

        assign_to : Symbol, MatrixSymbol, or string (optional)
            If provided, the printed code will set the expression to a
            variable with name ``assign_to``.
        """
        from sympy.matrices.expressions.matexpr import MatrixSymbol

        if isinstance(assign_to, string_types):
            if expr.is_Matrix:
                assign_to = MatrixSymbol(assign_to, *expr.shape)
            else:
                assign_to = Symbol(assign_to)
        elif not isinstance(assign_to, (Basic, type(None))):
            raise TypeError("{0} cannot assign to object of type {1}".format(
                    type(self).__name__, type(assign_to)))

        if assign_to:
            expr = Assignment(assign_to, expr)
        else:
            # _sympify is not enough b/c it errors on iterables
            expr = sympify(expr)

        # keep a set of expressions that are not strictly translatable to Code
        # and number constants that must be declared and initialized
        self._not_supported = set()
        self._number_symbols = set()

        if isinstance(expr, Eq):
            expr = Assignment(expr.lhs, expr.rhs)

        lines = self._print(expr).splitlines()

        # format the output
        if self._settings["human"]:
            frontlines = []
            if len(self._not_supported) > 0:
                frontlines.append(self._get_comment(
                        "Not supported in {0}:".format(self.language)))
                for expr in sorted(self._not_supported, key=str):
                    frontlines.append(self._get_comment(type(expr).__name__))
            for name, value in sorted(self._number_symbols, key=str):
                frontlines.append(self._declare_number_const(name, value))
            lines = frontlines + lines
            lines = self._format_code(lines)
            result = "\n".join(lines)
        else:
            lines = self._format_code(lines)
            result = (self._number_symbols, self._not_supported,
                    "\n".join(lines))
        return result

    def _get_expression_indices(self, expr, assign_to):
        from sympy.tensor import get_indices
        from sympy.tensor.indexed import Idx

        # need to remove not Idx indices !!!
        rinds = expr.atoms(Idx)
        linds = assign_to.atoms(Idx)

        # support broadcast of scalar
        if linds and not rinds:
            rinds = linds
        return self._sort_optimized(rinds, assign_to)

    def _rate_index_position(self, p):
        return p*5

    def _get_statement(self, codestring):
        return "%s" % codestring

    def _get_comment(self, text):
        return "# {0}".format(text)

    def _declare_number_const(self, name, value):
        return "cdef double const {0} = {1}".format(name, value)

    def _format_code(self, lines):
        return self.indent_code(lines)

    def _traverse_matrix_indices(self, mat):
        rows, cols = mat.shape
        return ((i, j) for i in range(rows) for j in range(cols))

    def _print_Pow(self, expr):
        if "Pow" in self.known_functions:
            return self._print_Function(expr)
        PREC = precedence(expr)
        if expr.exp == -1:
            return '1.0/%s' % (self.parenthesize(expr.base, PREC))
        elif expr.exp == 0.5:
            return 'sqrt(%s)' % self._print(expr.base)
        elif expr.exp > 0 and expr.exp.is_integer:
            line = "%s"%(self.parenthesize(expr.base, PREC))
            for i in range(1, expr.exp):
                line += '*' + "%s"%(self.parenthesize(expr.base, PREC))
            return '(%s)' % line
        else:
            return 'pow(%s, %s)' % (self._print(expr.base),
                                 self._print(expr.exp))

    def _print_Rational(self, expr):
        return self._print(expr.evalf(self._settings["precision"]))

    def _print_Indexed(self, expr):
        elem = []
        for i in range(expr.rank):
            elem.append(self._print(expr.indices[i]))
        return "%s[%s]" % (self._print(expr.base.label), ', '.join(elem))

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

    def _print_If(self, expr):
        lines = []
        for c, e in expr.statement:
            #temp1, temp2, cond = self.doprint(c)
            c = c.replace(Eq, AssignmentIf)
            lines.append("if %s:"%self._print(c))
            for ee in e:
                temp1, temp2, output = self.doprint(ee)
                lines.append(output)
        lines.append("#end")
        return "\n".join(lines)

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
        if expr.has(Assignment):
            for i, (e, c) in enumerate(expr.args):
                if i == 0:
                    lines.append("if (%s):" % self._print(c))
                elif i == len(expr.args) - 1 and c == True:
                    lines.append("else:")
                else:
                    lines.append("else if (%s):" % self._print(c))
                code0 = self._print(e)
                lines.append(code0)
                lines.append("}")
            return "\n".join(lines)
        else:
            # The piecewise was used in an expression, need to do inline
            # operators. This has the downside that inline operators will
            # not work for statements that span multiple lines (Matrix or
            # Indexed expressions).
            ecpairs = ["((%s) ? (\n%s\n)\n" % (self._print(c), self._print(e))
                    for e, c in expr.args[:-1]]
            last_line = ": (\n%s\n)" % self._print(expr.args[-1].expr)
            return ": ".join(ecpairs) + last_line + " ".join([")"*len(ecpairs)])

    def _print_ITE(self, expr):
        from sympy.functions import Piecewise
        _piecewise = Piecewise((expr.args[1], expr.args[0]), (expr.args[2], True))
        return self._print(_piecewise)

    def _print_MatrixElement(self, expr):
        return "{0}[{1}]".format(expr.parent, expr.j +
                expr.i*expr.parent.shape[1])

    def _print_Symbol(self, expr):

        name = super(CythonCodePrinter, self)._print_Symbol(expr)

        if expr in self._dereference:
            return '(*{0})'.format(name)
        else:
            return name

    def _print_For(self, expr):
        target = self._print(expr.target)
        if isinstance(expr.iterable, Range):
            start, stop, step = expr.iterable.args
        else:
            raise NotImplementedError("Only iterable currently supported is Range")
        body = self._print(expr.body)
        return ('for {target} in range({start}, {stop}, {step}):\n{body}\n#end\n').format(target=target, start=start,
                stop=stop, step=step, body=body)

    def _print_Assignment(self, expr):
        from sympy.functions.elementary.piecewise import Piecewise
        from sympy.matrices.expressions.matexpr import MatrixSymbol
        from sympy.matrices import MatrixBase, MatrixSlice
        from sympy.tensor.indexed import IndexedBase
        lhs = expr.lhs
        rhs = expr.rhs
        # We special case assignments that take multiple lines
        if isinstance(expr.rhs, Piecewise):
            # Here we modify Piecewise so each expression is now
            # an Assignment, and then continue on the print.
            expressions = []
            conditions = []
            for (e, c) in rhs.args:
                expressions.append(Assignment(lhs, e))
                conditions.append(c)
            temp = Piecewise(*zip(expressions, conditions))
            return self._print(temp)
        #elif isinstance(lhs, MatrixSymbol):
        elif isinstance(lhs, (MatrixBase, MatrixSymbol, MatrixSlice)):
            # Here we form an Assignment for each element in the array,
            # printing each one.
            lines = []
            for (i, j) in self._traverse_matrix_indices(lhs):
                temp = Assignment(lhs[i, j], rhs[i, j])
                code0 = self._print(temp)
                lines.append(code0)
            return "\n".join(lines)
        else:
            lhs_code = self._print(lhs)
            rhs_code = self._print(rhs)
            # hack to avoid the printing of m[i] = m[i]
            if lhs_code == rhs_code:
                return ""
            return self._get_statement("%s = %s" % (lhs_code, rhs_code))

    def _print_AssignmentIf(self, expr):
        from sympy.functions.elementary.piecewise import Piecewise
        from sympy.matrices.expressions.matexpr import MatrixSymbol
        from sympy.matrices import MatrixBase
        from sympy.tensor.indexed import IndexedBase
        lhs = expr.lhs
        rhs = expr.rhs
        # We special case assignments that take multiple lines
        if isinstance(expr.rhs, Piecewise):
            # Here we modify Piecewise so each expression is now
            # an Assignment, and then continue on the print.
            expressions = []
            conditions = []
            for (e, c) in rhs.args:
                expressions.append(Assignment(lhs, e))
                conditions.append(c)
            temp = Piecewise(*zip(expressions, conditions))
            return self._print(temp)
        elif isinstance(lhs, (MatrixBase, MatrixSymbol)):
            # Here we form an Assignment for each element in the array,
            # printing each one.
            lines = []
            for (i, j) in self._traverse_matrix_indices(lhs):
                temp = Assignment(lhs[i, j], rhs[i, j])
                code0 = self._print(temp)
                lines.append(code0)
            return "\n".join(lines)
        else:
            lhs_code = self._print(lhs)
            rhs_code = self._print(rhs)
            # hack to avoid the printing of m[i] = m[i]
            if lhs_code == rhs_code:
                return ""
            return self._get_statement("%s == %s" % (lhs_code, rhs_code))

    def _print_sign(self, func):
        return '((({0}) > 0) - (({0}) < 0))'.format(self._print(func.args[0]))

    def indent_code(self, code):
        """Accepts a string of code or a list of code lines"""

        if isinstance(code, string_types):
            code_lines = self.indent_code(code.splitlines(True))
            return ''.join(code_lines)

        tab = "    "
        inc_token = (':', ":\n")
        dec_token = ('end\n',)

        code = [ line.lstrip(' \t') for line in code ]


        increase = [ int(any(map(line.endswith, inc_token))) for line in code ]
        decrease = [ int(any(map(line.endswith, dec_token)))
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

    def _print_MatrixBase(self, A):
        # Handle zero dimensions:
        if A.rows == 0 or A.cols == 0:
            return 'zeros(%s, %s)' % (A.rows, A.cols)
        elif (A.rows, A.cols) == (1, 1):
            return "[%s]" % A[0, 0]
        elif A.rows == 1:
            return "[%s]" % A.table(self, rowstart='', rowend='', colsep=' ')
        elif A.cols == 1:
            # note .table would unnecessarily equispace the rows
            return "[%s]" % ", ".join([self._print(a) for a in A])
        return "[%s]" % A.table(self, rowstart='', rowend='',
                                rowsep=';\n', colsep=' ')

    _print_Matrix = \
        _print_DenseMatrix = \
        _print_MutableDenseMatrix = \
        _print_ImmutableMatrix = \
        _print_ImmutableDenseMatrix = \
        _print_MatrixBase

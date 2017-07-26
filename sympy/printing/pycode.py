from .precedence import precedence
from .codeprinter import CodePrinter

_kw_py2and3 = {
    'and', 'as', 'assert', 'break', 'class', 'continue', 'def', 'del', 'elif',
    'else', 'except', 'finally', 'for', 'from', 'global', 'if', 'import', 'in',
    'is', 'lambda', 'not', 'or', 'pass', 'raise', 'return', 'try', 'while',
    'with', 'yield', 'None'  # 'None' is actually not in Python 2's keyword.kwlist
}
_kw_only_py2 = {'exec', 'print'}
_kw_only_py3 = {'False', 'nonlocal', 'True'}

_known_functions = {
    'Abs': 'abs'
}

class PythonCodePrinter(CodePrinter):
    printmethod = "_pythoncode"
    language = "Python"
    standard = "python3"
    reserved_words = _kw_py2and3.union(_kw_only_py3)
    tab = '    '
    _kf = _known_functions
    _operators = {'and': 'and', 'or': 'or', 'not': 'not'}
    _default_settings = dict(CodePrinter._default_settings, precision=15)

    def __init__(self, settings=None):
        super(PythonCodePrinter, self).__init__(settings)
        self.known_functions = dict(self._kf, **(settings or {}).get(
            'user_functions', {}))

    def _format_code(self, lines):
        return lines

    def _get_comment(self, text):
        return "  # {0}".format(text)

    def _print_Mod(self, expr):
        PREC = precedence(expr)
        return ('{0} % {1}'.format(*map(lambda x: self.parenthesize(x, PREC), expr.args)))

    def _print_Piecewise(self, expr):
        lines = []
        for i, (e, c) in enumerate(expr.args):
            if i == 0:
                lines.append("if %s:" % self._print(c))
            elif i == len(expr.args) - 1 and c == True:
                lines.append('else:')
            else:
                lines.append('elif %s:' % self._print(c))
            lines.append(self.tab + 'return ' + self._print(e))
            if i == len(expr.args) - 1 and c != True:
                lines.append('else:')
                lines.append('%sraise NotImplementedError("Unhandled condition in: %s")' % (
                    self.tab, expr))
        return '\n'.join(lines)

    def _print_Sum(self, expr):
        loops = (
            'for {i} in range({a}, {b}+1)'.format(
                i=self._print(i),
                a=self._print(a),
                b=self._print(b))
            for i, a, b in expr.limits)
        return '(builtins.sum({function} {loops}))'.format(
            function=self._print(expr.function),
            loops=' '.join(loops))

    def _print_ImaginaryUnit(self, expr):
        return '1j'


def pycode(expr, **settings):
    return PythonCodePrinter(settings).doprint(expr)


class MpmathPrinter(PythonCodePrinter):
    """
    Lambda printer for mpmath which maintains precision for floats
    """
    def _print_Integer(self, e):
        return 'mpf(%d)' % e

    def _print_Float(self, e):
        # XXX: This does not handle setting mpmath.mp.dps. It is assumed that
        # the caller of the lambdified function will have set it to sufficient
        # precision to match the Floats in the expression.

        # Remove 'mpz' if gmpy is installed.
        args = str(tuple(map(int, e._mpf_)))
        return 'mpf(%s)' % args

    def _print_uppergamma(self,e): #printer for the uppergamma function
        return "gammainc({0}, {1}, inf)".format(self._print(e.args[0]), self._print(e.args[1]))

    def _print_lowergamma(self,e): #printer for the lowergamma functioin
        return "gammainc({0}, 0, {1})".format(self._print(e.args[0]), self._print(e.args[1]))


class NumPyPrinter(PythonCodePrinter):
    """
    Numpy printer which handles vectorized piecewise functions,
    logical operators, etc.
    """

    def _print_seq(self, seq, delimiter=', '):
        "General sequence printer: converts to tuple"
        # Print tuples here instead of lists because numba supports
        #     tuples in nopython mode.
        return '({},)'.format(delimiter.join(self._print(item) for item in seq))

    def _print_MatMul(self, expr):
        "Matrix multiplication printer"
        return '({0})'.format(').dot('.join(self._print(i) for i in expr.args))

    def _print_DotProduct(self, expr):
        # DotProduct allows any shape order, but numpy.dot does matrix
        # multiplication, so we have to make sure it gets 1 x n by n x 1.
        arg1, arg2 = expr.args
        if arg1.shape[0] != 1:
            arg1 = arg1.T
        if arg2.shape[1] != 1:
            arg2 = arg2.T

        return "numpy.dot(%s, %s)" % (self._print(arg1), self._print(arg2))

    def _print_Piecewise(self, expr):
        "Piecewise function printer"
        exprs = '[{0}]'.format(','.join(self._print(arg.expr) for arg in expr.args))
        conds = '[{0}]'.format(','.join(self._print(arg.cond) for arg in expr.args))
        # If [default_value, True] is a (expr, cond) sequence in a Piecewise object
        #     it will behave the same as passing the 'default' kwarg to select()
        #     *as long as* it is the last element in expr.args.
        # If this is not the case, it may be triggered prematurely.
        return 'numpy.select({0}, {1}, default=numpy.nan)'.format(conds, exprs)

    def _print_Relational(self, expr):
        "Relational printer for Equality and Unequality"
        op = {
            '==' :'equal',
            '!=' :'not_equal',
            '<'  :'less',
            '<=' :'less_equal',
            '>'  :'greater',
            '>=' :'greater_equal',
        }
        if expr.rel_op in op:
            lhs = self._print(expr.lhs)
            rhs = self._print(expr.rhs)
            return 'numpy.{op}({lhs}, {rhs})'.format(op=op[expr.rel_op],
                                               lhs=lhs,
                                               rhs=rhs)
        return super(NumPyPrinter, self)._print_Relational(expr)

    def _print_And(self, expr):
        "Logical And printer"
        # We have to override LambdaPrinter because it uses Python 'and' keyword.
        # If LambdaPrinter didn't define it, we could use StrPrinter's
        # version of the function and add 'logical_and' to NUMPY_TRANSLATIONS.
        return '{0}({1})'.format('numpy.logical_and', ','.join(self._print(i) for i in expr.args))

    def _print_Or(self, expr):
        "Logical Or printer"
        # We have to override LambdaPrinter because it uses Python 'or' keyword.
        # If LambdaPrinter didn't define it, we could use StrPrinter's
        # version of the function and add 'logical_or' to NUMPY_TRANSLATIONS.
        return '{0}({1})'.format('numpy.logical_or', ','.join(self._print(i) for i in expr.args))

    def _print_Not(self, expr):
        "Logical Not printer"
        # We have to override LambdaPrinter because it uses Python 'not' keyword.
        # If LambdaPrinter didn't define it, we would still have to define our
        #     own because StrPrinter doesn't define it.
        return '{0}({1})'.format('numpy.logical_not', ','.join(self._print(i) for i in expr.args))

    def _print_Min(self, expr):
        return '{0}(({1}))'.format('numpy.amin', ','.join(self._print(i) for i in expr.args))

    def _print_Max(self, expr):
        return '{0}(({1}))'.format('numpy.amax', ','.join(self._print(i) for i in expr.args))

    def _print_Pow(self, expr):
        if expr.exp == 0.5:
            return '{0}({1})'.format('sqrt', self._print(expr.base))
        else:
            return super(NumPyPrinter, self)._print_Pow(expr)

    def _print_log10(self, expr):  # log10 in C89, but type-generic macro in C99
        return 'numpy.log10({0})'.format(self._print(expr.args[0]))

    def _print_Sqrt(self, expr):
        return 'numpy.sqrt({0})'.format(self._print(expr.args[0]))

    def _print_hypot(self, expr):
        return 'numpy.hypot({0}, {1})'.format(*map(self._print, expr.args))

    def _print_expm1(self, expr):
        return 'numpy.expm1({0})'.format(self._print(expr.args[0]))

    def _print_log1p(self, expr):
        return 'numpy.log1p({0})'.format(self._print(expr.args[0]))

    def _print_exp2(self, expr):
        return 'numpy.exp2({0})'.format(self._print(expr.args[0]))

    def _print_log2(self, expr):
        return 'log2({0})'.format(self._print(expr.args[0]))

    def _print_acos(self, expr):
        return "numpy.arccos(%s)" % self._print(expr.args[0])

    def _print_acosh(self, expr):
        return "numpy.arccosh(%s)" % self._print(expr.args[0])

    def _print_arg(self, expr):
        return "numpy.angle(%s)" % self._print(expr.args[0])

    def _print_asin(self, expr):
        return "numpy.arcsin(%s)" % self._print(expr.args[0])

    def _print_asinh(self, expr):
        return "numpy.arcsinh(%s)" % self._print(expr.args[0])

    def _print_atan(self, expr):
        return "numpy.arctan(%s)" % self._print(expr.args[0])

    def _print_atan2(self, expr):
        return "numpy.arctan2(%s)" % ', '.join(map(self._print, expr.args))

    def _print_atanh(self, expr):
        return "numpy.arctanh(%s)" % self._print(expr.args[0])

    def _print_ceiling(self, expr):
        return "numpy.ceil(%s)" % self._print(expr.args[0])

    def _print_E(self, expr):
        return "numpy.e"

    def _print_im(self, expr):
        return "numpy.imag(%s)" % self._print(expr.args[0])

    def _print_ln(self, expr):
        return "numpy.log(%s)" % self._print(expr.args[0])

    def _print_Mod(self, expr):
        return "numpy.mod(%s)" % ', '.join(map(self._print, expr.args))

    def _print_oo(self, expr):
        return "numpy.inf"

    def _print_pi(self, expr):
        return "numpy.pi"

    def _print_re(self, expr):
        return "numpy.real(%s)" % self._print(expr.args[0])

    def _print_Matrix(self, expr):
        return "numpy.array(%s)" % self._print(expr.tolist())

    _print_SparseMatrix = _print_Matrix
    _print_ImmutableSparseMatrix = _print_Matrix
    _print_MutableDenseMatrix = _print_Matrix
    _print_ImmutableDenseMatrix = _print_Matrix


class SciPyPrinter(NumPyPrinter):

    def _print_SparseMatrix(self, expr):
        i, j, data = [], [], []
        for (r, c), v in expr._smat.items():
            i.append(r)
            j.append(c)
            data.append(v)

        return "scipy.sparse.coo_matrix({data}, ({i}, {j}), shape={shape})".format(
            data=data,
            i=i,
            j=j,
            shape=expr.shape
        )

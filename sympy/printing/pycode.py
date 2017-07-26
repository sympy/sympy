from functools import wraps
from itertools import chain
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
    'Abs': 'abs',
}
_known_functions_math = {
    'acos': 'acos',
    'acosh': 'acosh',
    'asin': 'asin',
    'asinh': 'asinh',
    'atan': 'atan',
    'atan2': 'atan2',
    'atanh': 'atanh',
    'ceiling': 'ceil',
    'cos': 'cos',
    'cosh': 'cosh',
    'erf': 'erf',
    'erfc': 'erfc',
    'exp': 'exp',
    'expm1': 'expm1',
    'factorial': 'factorial',
    'floor': 'floor',
    'gamma': 'gamma',
    'hypot': 'hypot',
    'loggamma': 'lgamma',
    'log': 'log',
    'log10': 'log10',
    'log1p': 'log1p',
    'log2': 'log2',
    'sin': 'sin',
    'sinh': 'sinh',
    'Sqrt': 'sqrt',
    'tan': 'tan',
    'tanh': 'tanh'
}  # Not used from ``math``: [copysign isclose isfinite isinf isnan ldexp frexp pow modf
# radians trunc fmod fsum gcd degrees fabs]

class requires(object):
    """ Decorator for registering requirements on print methods. """
    def __init__(self, **kwargs):
        self._req = kwargs

    def __call__(self, method):
        def _method_wrapper(self_, *args, **kwargs):
            for k, v in self._req.items():
                getattr(self_, k).update(v)
            return method(self_, *args, **kwargs)
        return wraps(method)(_method_wrapper)


def _print_known(self, expr):
    known = self.known_functions[expr.__class__.__name__]
    parts = known.split('.')
    if len(parts) > 1:
        self.modules.update(('.'.join(parts[:-1]),))
    return '{name}({args})'.format(name=known, args=', '.join(map(self._print, expr.args)))


class PythonCodePrinter(CodePrinter):
    printmethod = "_pythoncode"
    language = "Python"
    standard = "python3"
    reserved_words = _kw_py2and3.union(_kw_only_py3)
    modules = None  # initialized to a set in __init__
    tab = '    '
    _kf = dict(chain(
        _known_functions.items(),
        [(k, 'math.' + v) for k, v in _known_functions_math.items()]
    ))
    _operators = {'and': 'and', 'or': 'or', 'not': 'not'}
    _default_settings = dict(CodePrinter._default_settings, precision=15)

    def __init__(self, settings=None):
        super(PythonCodePrinter, self).__init__(settings)
        self.modules = set()
        self.known_functions = dict(self._kf, **(settings or {}).get(
            'user_functions', {}))

    def _format_code(self, lines):
        return lines

    def _get_comment(self, text):
        return "  # {0}".format(text)

    @requires(modules={'math'})
    def _print_Exp1(self, expr):
        return 'math.e'

    @requires(modules={'math'})
    def _print_pi(self, expr):
        return 'math.pi'

    @requires(modules={'math'})
    def _print_Infinity(self, expr):
        return 'math.inf'

    @requires(modules={'math'})
    def _print_NaN(self, expr):
        return 'math.nan'

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

for k in PythonCodePrinter._kf:
    setattr(PythonCodePrinter, '_print_%s' % k, _print_known)


def pycode(expr, **settings):
    return PythonCodePrinter(settings).doprint(expr)


class MpmathPrinter(PythonCodePrinter):
    """
    Lambda printer for mpmath which maintains precision for floats
    """

    @requires(modules={'mpmath'})
    def _print_Integer(self, e):
        return 'mpf(%d)' % e

    def _print_Float(self, e):
        # XXX: This does not handle setting mpmath.mp.dps. It is assumed that
        # the caller of the lambdified function will have set it to sufficient
        # precision to match the Floats in the expression.

        # Remove 'mpz' if gmpy is installed.
        args = str(tuple(map(int, e._mpf_)))
        return 'mpf(%s)' % args

    @requires(modules={'mpmath'})
    def _print_uppergamma(self,e): #printer for the uppergamma function
        return "gammainc({0}, {1}, inf)".format(self._print(e.args[0]), self._print(e.args[1]))

    @requires(modules={'mpmath'})
    def _print_lowergamma(self,e): #printer for the lowergamma functioin
        return "gammainc({0}, 0, {1})".format(self._print(e.args[0]), self._print(e.args[1]))


_not_in_numpy = 'erf erfc factorial gamma lgamma'.split()
_in_numpy = [(k, v) for k, v in _known_functions_math.items() if k not in _not_in_numpy]
_known_functions_numpy = dict(_in_numpy, **{
    'acos': 'arccos',
    'acosh': 'arccosh',
    'asin': 'arcsin',
    'asinh': 'arcsinh',
    'atan': 'arctan',
    'atan2': 'arctan2',
    'atanh': 'arctanh',
})


class NumPyPrinter(PythonCodePrinter):
    """
    Numpy printer which handles vectorized piecewise functions,
    logical operators, etc.
    """
    _kf = dict(chain(
        PythonCodePrinter._kf.items(),
        [(k, 'numpy.' + v) for k, v in _known_functions_numpy.items()]
    ))


    @requires(modules={'numpy'})
    def _print_Exp1(self, expr):
        return "numpy.e"

    @requires(modules={'numpy'})
    def _print_Infinity(self, expr):
        return "numpy.inf"

    @requires(modules={'numpy'})
    def _print_pi(self, expr):
        return "numpy.pi"

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
            return '{0}({1})'.format('numpy.sqrt', self._print(expr.base))
        else:
            return super(NumPyPrinter, self)._print_Pow(expr)

    def _print_arg(self, expr):
        return "numpy.angle(%s)" % self._print(expr.args[0])

    def _print_im(self, expr):
        return "numpy.imag(%s)" % self._print(expr.args[0])

    def _print_Mod(self, expr):
        return "numpy.mod(%s)" % ', '.join(map(self._print, expr.args))

    def _print_re(self, expr):
        return "numpy.real(%s)" % self._print(expr.args[0])

    def _print_Matrix(self, expr):
        return "numpy.array(%s)" % self._print(expr.tolist())

    _print_SparseMatrix = _print_Matrix
    _print_ImmutableSparseMatrix = _print_Matrix
    _print_MutableDenseMatrix = _print_Matrix
    _print_ImmutableDenseMatrix = _print_Matrix

for k in NumPyPrinter._kf:
    setattr(NumPyPrinter, '_print_%s' % k, _print_known)


_known_functions_scipy_special = {
    'erf': 'erf',
    'erfc': 'erfc',
    'gamma': 'gamma',
    'loggamma': 'gammaln'
}

class SciPyPrinter(NumPyPrinter):

    _kf = dict(chain(
        NumPyPrinter._kf.items(),
        [(k, 'scipy.special.' + v) for k, v in _known_functions_scipy_special.items()]
    ))

    @requires(modules={'scipy.sparse'})
    def _print_SparseMatrix(self, expr):
        i, j, data = [], [], []
        for (r, c), v in expr._smat.items():
            i.append(r)
            j.append(c)
            data.append(v)

        return "scipy.sparse.coo_matrix({data}, ({i}, {j}), shape={shape})".format(
            data=data, i=i, j=j, shape=expr.shape
        )

    _print_ImmutableSparseMatrix = _print_SparseMatrix
    _print_MutableSparseMatrix = _print_SparseMatrix


for k in SciPyPrinter._kf:
    setattr(SciPyPrinter, '_print_%s' % k, _print_known)


class SymPyPrinter(PythonCodePrinter):

    _kf = dict([(k, 'sympy.' + v) for k, v in chain(
        _known_functions.items(),
        _known_functions_math.items()
    )])

    @requires(modules={'sympy'})
    def _print_Function(self, expr):
        return 'sympy.{name}({args})'.format(name=expr.func.__name__, args=', '.join(map(self._print, expr.args)))

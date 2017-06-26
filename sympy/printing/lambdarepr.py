from __future__ import print_function, division

from .str import StrPrinter
from sympy.utilities import default_sort_key


class LambdaPrinter(StrPrinter):
    """
    This printer converts expressions into strings that can be used by
    lambdify.
    """

    def _print_MatrixBase(self, expr):
        return "%s(%s)" % (expr.__class__.__name__,
                           self._print((expr.tolist())))

    _print_SparseMatrix = \
        _print_MutableSparseMatrix = \
        _print_ImmutableSparseMatrix = \
        _print_Matrix = \
        _print_DenseMatrix = \
        _print_MutableDenseMatrix = \
        _print_ImmutableMatrix = \
        _print_ImmutableDenseMatrix = \
        _print_MatrixBase

    def _print_Piecewise(self, expr):
        result = []
        i = 0
        for arg in expr.args:
            e = arg.expr
            c = arg.cond
            result.append('((')
            result.append(self._print(e))
            result.append(') if (')
            result.append(self._print(c))
            result.append(') else (')
            i += 1
        result = result[:-1]
        result.append(') else None)')
        result.append(')'*(2*i - 2))
        return ''.join(result)

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

    def _print_And(self, expr):
        result = ['(']
        for arg in sorted(expr.args, key=default_sort_key):
            result.extend(['(', self._print(arg), ')'])
            result.append(' and ')
        result = result[:-1]
        result.append(')')
        return ''.join(result)

    def _print_Or(self, expr):
        result = ['(']
        for arg in sorted(expr.args, key=default_sort_key):
            result.extend(['(', self._print(arg), ')'])
            result.append(' or ')
        result = result[:-1]
        result.append(')')
        return ''.join(result)

    def _print_Not(self, expr):
        result = ['(', 'not (', self._print(expr.args[0]), '))']
        return ''.join(result)

    def _print_BooleanTrue(self, expr):
        return "True"

    def _print_BooleanFalse(self, expr):
        return "False"

    def _print_ITE(self, expr):
        result = [
            '((', self._print(expr.args[1]),
            ') if (', self._print(expr.args[0]),
            ') else (', self._print(expr.args[2]), '))'
        ]
        return ''.join(result)

class TensorflowPrinter(LambdaPrinter):
    """
    Tensorflow printer which handles vectorized piecewise functions,
    logical operators, max/min, and relational operators.
    """

    def _print_And(self, expr):
        "Logical And printer"
        # We have to override LambdaPrinter because it uses Python 'and' keyword.
        # If LambdaPrinter didn't define it, we could use StrPrinter's
        # version of the function and add 'logical_and' to TENSORFLOW_TRANSLATIONS.
        return '{0}({1})'.format('logical_and', ','.join(self._print(i) for i in expr.args))

    def _print_Or(self, expr):
        "Logical Or printer"
        # We have to override LambdaPrinter because it uses Python 'or' keyword.
        # If LambdaPrinter didn't define it, we could use StrPrinter's
        # version of the function and add 'logical_or' to TENSORFLOW_TRANSLATIONS.
        return '{0}({1})'.format('logical_or', ','.join(self._print(i) for i in expr.args))

    def _print_Not(self, expr):
        "Logical Not printer"
        # We have to override LambdaPrinter because it uses Python 'not' keyword.
        # If LambdaPrinter didn't define it, we would still have to define our
        #     own because StrPrinter doesn't define it.
        return '{0}({1})'.format('logical_not', ','.join(self._print(i) for i in expr.args))

    def _print_Min(self, expr, **kwargs):
        from sympy import Min
        if len(expr.args) == 1:
            return self._print(expr.args[0], **kwargs)

        return 'minimum({0}, {1})'.format(
            self._print(expr.args[0], **kwargs),
            self._print(Min(*expr.args[1:]), **kwargs))

    def _print_Max(self, expr, **kwargs):
        from sympy import Max
        if len(expr.args) == 1:
            return self._print(expr.args[0], **kwargs)

        return 'maximum({0}, {1})'.format(
            self._print(expr.args[0], **kwargs),
            self._print(Max(*expr.args[1:]), **kwargs))

    def _print_Piecewise(self, expr, **kwargs):
        from sympy import Piecewise
        e, cond = expr.args[0].args
        if len(expr.args) == 1:
            return 'select({0}, {1}, {2})'.format(
                self._print(cond, **kwargs),
                self._print(e, **kwargs),
                0)

        return 'select({0}, {1}, {2})'.format(
            self._print(cond, **kwargs),
            self._print(e, **kwargs),
            self._print(Piecewise(*expr.args[1:]), **kwargs))

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
            return '{op}({lhs}, {rhs})'.format(op=op[expr.rel_op],
                                               lhs=lhs,
                                               rhs=rhs)
        return super(TensorflowPrinter, self)._print_Relational(expr)


class NumPyPrinter(LambdaPrinter):
    """
    Numpy printer which handles vectorized piecewise functions,
    logical operators, etc.
    """
    _default_settings = {
        "order": "none",
        "full_prec": "auto",
    }

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

        return "dot(%s, %s)" % (self._print(arg1), self._print(arg2))

    def _print_Piecewise(self, expr):
        "Piecewise function printer"
        exprs = '[{0}]'.format(','.join(self._print(arg.expr) for arg in expr.args))
        conds = '[{0}]'.format(','.join(self._print(arg.cond) for arg in expr.args))
        # If [default_value, True] is a (expr, cond) sequence in a Piecewise object
        #     it will behave the same as passing the 'default' kwarg to select()
        #     *as long as* it is the last element in expr.args.
        # If this is not the case, it may be triggered prematurely.
        return 'select({0}, {1}, default=nan)'.format(conds, exprs)

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
            return '{op}({lhs}, {rhs})'.format(op=op[expr.rel_op],
                                               lhs=lhs,
                                               rhs=rhs)
        return super(NumPyPrinter, self)._print_Relational(expr)

    def _print_And(self, expr):
        "Logical And printer"
        # We have to override LambdaPrinter because it uses Python 'and' keyword.
        # If LambdaPrinter didn't define it, we could use StrPrinter's
        # version of the function and add 'logical_and' to NUMPY_TRANSLATIONS.
        return '{0}({1})'.format('logical_and', ','.join(self._print(i) for i in expr.args))

    def _print_Or(self, expr):
        "Logical Or printer"
        # We have to override LambdaPrinter because it uses Python 'or' keyword.
        # If LambdaPrinter didn't define it, we could use StrPrinter's
        # version of the function and add 'logical_or' to NUMPY_TRANSLATIONS.
        return '{0}({1})'.format('logical_or', ','.join(self._print(i) for i in expr.args))

    def _print_Not(self, expr):
        "Logical Not printer"
        # We have to override LambdaPrinter because it uses Python 'not' keyword.
        # If LambdaPrinter didn't define it, we would still have to define our
        #     own because StrPrinter doesn't define it.
        return '{0}({1})'.format('logical_not', ','.join(self._print(i) for i in expr.args))

    def _print_Min(self, expr):
        return '{0}(({1}))'.format('amin', ','.join(self._print(i) for i in expr.args))

    def _print_Max(self, expr):
        return '{0}(({1}))'.format('amax', ','.join(self._print(i) for i in expr.args))

    def _print_Pow(self, expr):
        if expr.exp == 0.5:
            return '{0}({1})'.format('sqrt', self._print(expr.base))
        else:
            return super(NumPyPrinter, self)._print_Pow(expr)

    def _print_log10(self, expr):  # log10 in C89, but type-generic macro in C99
        return 'log10({0})'.format(self._print(expr.args[0]))

    def _print_Sqrt(self, expr):
        return 'sqrt({0})'.format(self._print(expr.args[0]))

    def _print_hypot(self, expr):
        return 'hypot({0}, {1})'.format(*map(self._print, expr.args))

    def _print_expm1(self, expr):
        return 'expm1({0})'.format(self._print(expr.args[0]))

    def _print_log1p(self, expr):
        return 'log1p({0})'.format(self._print(expr.args[0]))

    def _print_exp2(self, expr):
        return 'exp2({0})'.format(self._print(expr.args[0]))

    def _print_log2(self, expr):
        return 'log2({0})'.format(self._print(expr.args[0]))


# numexpr works by altering the string passed to numexpr.evaluate
# rather than by populating a namespace.  Thus a special printer...
class NumExprPrinter(LambdaPrinter):
    # key, value pairs correspond to sympy name and numexpr name
    # functions not appearing in this dict will raise a TypeError
    _numexpr_functions = {
        'sin' : 'sin',
        'cos' : 'cos',
        'tan' : 'tan',
        'asin': 'arcsin',
        'acos': 'arccos',
        'atan': 'arctan',
        'atan2' : 'arctan2',
        'sinh' : 'sinh',
        'cosh' : 'cosh',
        'tanh' : 'tanh',
        'asinh': 'arcsinh',
        'acosh': 'arccosh',
        'atanh': 'arctanh',
        'ln' : 'log',
        'log': 'log',
        'exp': 'exp',
        'sqrt' : 'sqrt',
        'Abs' : 'abs',
        'conjugate' : 'conj',
        'im' : 'imag',
        're' : 'real',
        'where' : 'where',
        'complex' : 'complex',
        'contains' : 'contains',
    }

    def _print_ImaginaryUnit(self, expr):
        return '1j'

    def _print_seq(self, seq, delimiter=', '):
        # simplified _print_seq taken from pretty.py
        s = [self._print(item) for item in seq]
        if s:
            return delimiter.join(s)
        else:
            return ""

    def _print_Function(self, e):
        func_name = e.func.__name__

        nstr = self._numexpr_functions.get(func_name, None)
        if nstr is None:
            # check for implemented_function
            if hasattr(e, '_imp_'):
                return "(%s)" % self._print(e._imp_(*e.args))
            else:
                raise TypeError("numexpr does not support function '%s'" %
                                func_name)
        return "%s(%s)" % (nstr, self._print_seq(e.args))

    def blacklisted(self, expr):
        raise TypeError("numexpr cannot be used with %s" %
                        expr.__class__.__name__)

    # blacklist all Matrix printing
    _print_SparseMatrix = \
    _print_MutableSparseMatrix = \
    _print_ImmutableSparseMatrix = \
    _print_Matrix = \
    _print_DenseMatrix = \
    _print_MutableDenseMatrix = \
    _print_ImmutableMatrix = \
    _print_ImmutableDenseMatrix = \
    blacklisted
    # blacklist some python expressions
    _print_list = \
    _print_tuple = \
    _print_Tuple = \
    _print_dict = \
    _print_Dict = \
    blacklisted

    def doprint(self, expr):
        lstr = super(NumExprPrinter, self).doprint(expr)
        return "evaluate('%s', truediv=True)" % lstr

class MpmathPrinter(LambdaPrinter):
    """
    Lambda printer for mpmath which maintains precision for floats
    """
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

def lambdarepr(expr, **settings):
    """
    Returns a string usable for lambdifying.
    """
    return LambdaPrinter(settings).doprint(expr)

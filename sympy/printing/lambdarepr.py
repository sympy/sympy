from __future__ import print_function, division

from .str import StrPrinter
from sympy.printing.precedence import precedence
from sympy.utilities import default_sort_key


def _find_first_symbol(expr):
    for atom in expr.atoms():
        if atom.is_Symbol:
            return atom
    raise ValueError('expression must contain a Symbol: %r' % expr)


class LambdaPrinter(StrPrinter):
    """
    This printer converts expressions into strings that can be used by
    lambdify.
    """

    def __init__(self, settings=None):
        if settings:
            self.use_numpy = settings.pop('use_numpy', False)
        else:
            self.use_numpy = False
        StrPrinter.__init__(self, settings)

    def _print_MatrixBase(self, expr):
        return "%s(%s)" % (expr.__class__.__name__, str(expr.tolist()))

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
        from sympy.core.sets import Interval
        result = []
        i = 0
        for arg in expr.args:
            e = arg.expr
            c = arg.cond
            result.append('((')
            result.append(self._print(e))
            result.append(') if (')
            if isinstance(c, Interval):
                result.append(self._print(c.contains(_find_first_symbol(e))))
            else:
                result.append(self._print(c))
            result.append(') else (')
            i += 1
        result = result[:-1]
        result.append(') else None)')
        result.append(')'*(2*i - 2))
        return ''.join(result)

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

    def _print_MatMul(self,expr):
        if self.use_numpy:
            from sympy.matrices.expressions import MatrixExpr
            from sympy.matrices.matrices import MatrixBase
            mats = '.dot'.join([self.parenthesize(arg, 1000) for arg
                in expr.args if isinstance(arg, (MatrixBase, MatrixExpr))])
            others = '*'.join([self.parenthesize(arg, precedence(expr))
                for arg in expr.args if not
                isinstance(arg, (MatrixBase, MatrixExpr))])
            if others:
                return others + "*" + mats
            else:
                return mats
        else:
            return '*'.join([self.parenthesize(arg, precedence(expr))
                for arg in expr.args])


def lambdarepr(expr, **settings):
    """
    Returns a string usable for lambdifying.
    """
    return LambdaPrinter(settings).doprint(expr)

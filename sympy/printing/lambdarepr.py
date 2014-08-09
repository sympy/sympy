from __future__ import print_function, division

from .str import StrPrinter
from sympy.utilities import default_sort_key


class LambdaPrinter(StrPrinter):
    """
    This printer converts expressions into strings that can be used by
    lambdify.
    """

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
        from sympy.sets.sets import Interval
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

# numexpr works by altering the string passed to numexpr.evaluate
# rather than by populating a namespace.  Thus a special printer...
# This class cannot be placed in lambdify.py due to circular imports
class NumExprPrinter(LambdaPrinter):
    # strings to substitute for sympy expressions
    str_subs = {
        "Abs": "abs",
        "acos": "arccos",
        "acosh": "arccosh",
        "asin": "arcsin",
        "asinh": "arcsinh",
        "atan": "arctan",
        "atan2": "arctan2",
        "atanh": "arctanh",
        "E": "e",
        "im": "imag",
        "ln": "log",
        "re": "real",
        "I": "1j",
    }
    # numexpr does not support these operations, throw TypeError if found
    blacklisted = ('matrix', 'list', 'tuple')
    def doprint(self, expr):
        lstr = super(NumExprPrinter, self).doprint(expr)
        # check blacklisted
        for b in self.blacklisted:
            if b in lstr.lower():
                raise TypeError("numexpr cannot be used with {}".format(b))
        # substitute strings
        for k, v in sorted(self.str_subs.items(), key=lambda x : -len(x[0])):
            lstr = lstr.replace(k, v)
        return "evaluate('"+lstr+"')"

def lambdarepr(expr, **settings):
    """
    Returns a string usable for lambdifying.
    """
    return LambdaPrinter(settings).doprint(expr)

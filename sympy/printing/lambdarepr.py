from str import StrPrinter
from sympy import S

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

    def _print_Matrix(self, expr):
        return "Matrix([%s])"%expr._format_str(self._print, ",")

    def _print_Piecewise(self, expr):
        from sympy.core.sets import Interval
        result = []
        i = 0
        for e, c in expr.exprcondpairs:
            result.append('((')
            result.append(self._print(e))
            result.append(') if (')
            if isinstance(c, Interval):
                result.append(self._print(c.contains(_find_first_symbol(e))))
            else:
                result.append(self._print(c))
            result.append(') else (')
            i += 1
        if len(result) > 0:
            result = result[:-1]
            result.append(') else ')
        if expr.otherwise is not S.NaN:
            result.append('(')
            result.append(self._print(expr.otherwise))
            result.append(')')
        else:
            result.append('None')
        result.append(')'*(2*i - 1))
        return ''.join(result)

    def _print_And(self, expr):
        result = ['(']
        for arg in expr.args:
            result.extend(['(', self._print(arg), ')'])
            result.append(' and ')
        result = result[:-1]
        result.append(')')
        return ''.join(result)

    def _print_Or(self, expr):
        result = ['(']
        for arg in expr.args:
            result.extend(['(', self._print(arg), ')'])
            result.append(' or ')
        result = result[:-1]
        result.append(')')
        return ''.join(result)

    def _print_Not(self, expr):
        result = ['(', 'not (', self._print(expr.args[0]), '))']
        return ''.join(result)

def lambdarepr(expr, **settings):
    """
    Returns a string usable for lambdifying.
    """
    return LambdaPrinter(settings).doprint(expr)

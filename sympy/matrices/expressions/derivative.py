from matexpr import MatrixExpr, ZeroMatrix
from sympy import Basic, Derivative, Tuple

class MatrixDerivative(MatrixExpr):
    """Matrix Derivative

    """
    is_Derivative = True
    def __new__(cls, expr, *variables, **assumptions):

        evaluate = assumptions.pop('evaluate', True)

        if not variables:
            variables = expr.free_symbols

        if evaluate:
            if hasattr(expr, '_eval_derivative'):
                for var in variables:
                    expr = expr._eval_derivative(var)
                return expr

        return Basic.__new__(cls, expr)

    @property
    def expr(self):
        return self.args[0]

    @property
    def shape(self):
        return self.expr.shape

    def _deep_symbol(self):
        if self.expr.is_Derivative:
            return self.expr._deep_symbol()
        if self.expr.is_MatrixSymbol:
            return self.expr
        raise NotImplementedError()

    def _eval_derivative(self, x):
        if x == self._deep_symbol():
            return MatrixDerivative(self, evaluate=False)
        return ZeroMatrix(*self.shape)

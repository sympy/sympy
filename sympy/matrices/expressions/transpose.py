from matexpr import MatrixExpr
from sympy import Basic

class Transpose(MatrixExpr):
    """Matrix Transpose

    Represents the transpose of a matrix expression.

    Use .T as shorthand

    >>> from sympy import MatrixSymbol, Transpose
    >>> A = MatrixSymbol('A', 3, 5)
    >>> B = MatrixSymbol('B', 5, 3)
    >>> Transpose(A)
    A'
    >>> A.T
    A'
    >>> Transpose(A*B)
    B'*A'
    """
    is_Transpose = True
    def __new__(cls, mat):
        try:
            return mat.transpose()
        except (AttributeError, NotImplementedError):
            return Basic.__new__(cls, mat)

    @property
    def arg(self):
        return self.args[0]

    @property
    def shape(self):
        return self.arg.shape[::-1]

    def _entry(self, i, j):
        return self.arg._entry(j, i)

    def _eval_transpose(self):
        return self.arg

    def _eval_trace(self):
        from trace import Trace
        return Trace(self.arg) # Trace(X.T) => Trace(X)

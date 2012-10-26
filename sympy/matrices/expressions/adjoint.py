from matexpr import MatrixExpr
from sympy import Basic


class Adjoint(MatrixExpr):
    """Matrix Adjoint

    Represents the Adjoint of a matrix expression.

    >>> from sympy import MatrixSymbol, Adjoint
    >>> A = MatrixSymbol('A', 3, 5)
    >>> B = MatrixSymbol('B', 5, 3)
    >>> Adjoint(A*B)
    Adjoint(B)*Adjoint(A)
    """
    is_Adjoint = True

    def __new__(cls, mat):
        try:
            return mat._eval_adjoint()
        except (AttributeError, NotImplementedError):
            pass
        try:
            return mat.adjoint()
        except (AttributeError, NotImplementedError):
            pass
        return Basic.__new__(cls, mat)

    @property
    def arg(self):
        return self.args[0]

    @property
    def shape(self):
        return self.arg.shape[::-1]

    def _entry(self, i, j):
        return conjugate(self.arg._entry(j, i))

    def _eval_adjoint(self):
        return self.arg

    def _eval_trace(self):
        from trace import Trace
        return conjugate(Trace(self.arg))

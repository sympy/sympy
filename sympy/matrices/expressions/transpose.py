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

        if not mat.is_Matrix:
            return mat

        if isinstance(mat, Transpose):
            return mat.arg

        if hasattr(mat, 'transpose'):
            return mat.transpose()

        if mat.is_Mul:
            return MatMul(*[Transpose(arg) for arg in mat.args[::-1]])

        if mat.is_Add:
            return MatAdd(*[Transpose(arg) for arg in mat.args])

        return Basic.__new__(cls, mat)

    @property
    def arg(self):
        return self.args[0]

    @property
    def shape(self):
        return self.arg.shape[::-1]

from matmul import MatMul
from matadd import MatAdd

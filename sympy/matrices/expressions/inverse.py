from sympy.core.sympify import _sympify

from sympy.matrices.expressions.matexpr import ShapeError
from sympy.matrices.expressions.matpow import MatPow


class Inverse(MatPow):
    """Matrix Inverse

    Represents the Inverse of a matrix expression

    Use .I as shorthand

    >>> from sympy import MatrixSymbol, Inverse
    >>> A = MatrixSymbol('A', 3, 3)
    >>> B = MatrixSymbol('B', 3, 3)
    >>> Inverse(A)
    A^-1
    >>> A.I
    A^-1
    >>> Inverse(A*B)
    B^-1*A^-1

    """
    is_Inverse = True

    def __new__(cls, mat):
        mat = _sympify(mat)
        if not mat.is_Matrix:
            return mat**(-1)
        if not mat.is_square:
            raise ShapeError("Inverse of non-square matrix %s" % mat)
        return MatPow.__new__(cls, mat, -1)

    @property
    def arg(self):
        return self.args[0]

    @property
    def shape(self):
        return self.arg.shape

    def _eval_inverse(self):
        return self.arg

    def doit(self, **hints):
        if hints.get('deep', True):
            return self.arg.doit(**hints).inverse()
        else:
            return self.arg.inverse()

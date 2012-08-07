from matexpr import ShapeError
from matpow import MatPow

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

    def __new__(cls, mat, **kwargs):

        if not mat.is_Matrix:
            return mat**(-1)

        if not mat.is_square:
            raise ShapeError("Inverse of non-square matrix %s"%mat)

        try:
            return mat._eval_inverse(**kwargs)
        except (AttributeError, NotImplementedError):
            return MatPow.__new__(cls, mat, -1)

    @property
    def arg(self):
        return self.args[0]

    @property
    def shape(self):
        return self.arg.shape

    def _eval_inverse(self):
        return self.arg

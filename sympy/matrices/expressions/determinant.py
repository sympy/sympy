from sympy import Basic, Expr
from matexpr import ShapeError


class Determinant(Expr):
    """Matrix Determinant

    Represents the determinant of a matrix expression.

    >>> from sympy import MatrixSymbol, Determinant, eye
    >>> A = MatrixSymbol('A', 3, 3)
    >>> Deteterminant(A)
    Determinant(A)

    >>> Determinant(eye(3))
    1
    """

    def __new__(cls, mat):
        if not mat.is_Matrix:
            raise TypeError("Input to Determinant, %s, not a matrix" % str(mat))

        if not mat.is_square:
            raise ShapeError("Det of a non-square matrix")

        try:
            return mat.det()
        except (AttributeError, NotImplementedError):
            return Basic.__new__(cls, mat)

    @property
    def arg(self):
        return self.args[0]

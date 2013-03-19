from sympy import Basic, Expr
from .matexpr import ShapeError


class Det(Expr):
    """Matrix Det

    Represents the Det of a matrix expression.

    >>> from sympy import MatrixSymbol, Det, eye
    >>> A = MatrixSymbol('A', 3, 3)
    >>> Det(A)
    Det(A)

    >>> Det(eye(3))
    1
    """
    is_Det = True

    def __new__(cls, mat):
        if not mat.is_Matrix:
            raise TypeError("input to Det, %s, is not a matrix" % str(mat))

        if not mat.is_square:
            raise ShapeError("Det of a non-square matrix")

        try:
            return mat.det()
        except (AttributeError, NotImplementedError):
            return Basic.__new__(cls, mat)

    @property
    def arg(self):
        return self.args[0]

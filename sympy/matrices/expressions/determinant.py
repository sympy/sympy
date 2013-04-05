from sympy import Basic, Expr, S, Q, ask
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

        return Basic.__new__(cls, mat)

    @property
    def arg(self):
        return self.args[0]

    def doit(self, expand=False):
        if ask(Q.singular(self)):
            return S.Zero
        try:
            return self.arg._eval_determinant()
        except (AttributeError, NotImplementedError):
            return self

def det(matexpr):
    return Determinant(matexpr).doit()

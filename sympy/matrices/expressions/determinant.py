from __future__ import print_function, division

from sympy import Basic, Expr, S, Q
from .matexpr import ShapeError


class Determinant(Expr):
    """Matrix Determinant

    Represents the determinant of a matrix expression.

    >>> from sympy import MatrixSymbol, Determinant, eye
    >>> A = MatrixSymbol('A', 3, 3)
    >>> Determinant(A)
    Determinant(A)

    >>> Determinant(eye(3)).doit()
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
        try:
            return self.arg._eval_determinant()
        except (AttributeError, NotImplementedError):
            return self

def det(matexpr):
    """ Matrix Determinant

    >>> from sympy import MatrixSymbol, det, eye
    >>> A = MatrixSymbol('A', 3, 3)
    >>> det(A)
    Determinant(A)

    >>> det(eye(3))
    1
    """

    return Determinant(matexpr).doit()

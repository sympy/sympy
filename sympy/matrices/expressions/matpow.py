from matexpr import MatrixExpr, ShapeError, Identity
from sympy import Pow, S
from sympy.core.sympify import _sympify

class MatPow(MatrixExpr, Pow):

    def __new__(cls, b, e):
        e = _sympify(e)
        if e is S.One or b.is_ZeroMatrix:
            return b
        elif not b.is_square:
            raise ShapeError("Power of non-square matrix %s"%b)
        elif e is S.Zero:
            return Identity(b.n)
        else:
            return MatrixExpr.__new__(cls, b, e)

    @property
    def shape(self):
        return self.base.shape

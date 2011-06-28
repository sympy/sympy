from matexpr import MatrixExpr, ShapeError
from sympy import Basic

class Inverse(MatrixExpr):

    def __new__(cls, mat):

        if not mat.is_Matrix:
            return mat**(-1)

        if hasattr(mat, 'inv'):
            try:
                return mat.inv()
            except:
                pass

        if isinstance(mat, Inverse):
            return mat.arg

        if not mat.is_square:
            raise ShapeError("Inverse of non-square matrix %s"%mat)

        return Basic.__new__(cls, mat)

    @property
    def arg(self):
        return self.args[0]

    @property
    def shape(self):
        return self.arg.shape[::-1]



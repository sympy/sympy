from matexpr import MatrixExpr, ShapeError
from sympy import Basic

class Transpose(MatrixExpr):
    is_Transpose = True
    def __new__(cls, mat):

        if not mat.is_Matrix:
            return mat

        if isinstance(mat, Transpose):
            return mat.arg

        if hasattr(mat, 'transpose'):
            try:
                return mat.transpose()
            except:
                pass

        return Basic.__new__(cls, mat)

    @property
    def arg(self):
        return self.args[0]

    @property
    def shape(self):
        return self.arg.shape[::-1]



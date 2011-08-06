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

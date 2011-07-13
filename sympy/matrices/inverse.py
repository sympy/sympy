from matexpr import MatrixExpr, ShapeError
from matpow import MatPow
from sympy import Basic

class Inverse(MatPow):
    is_Inverse = True

    def __new__(cls, mat):

        if not mat.is_Matrix:
            return mat**(-1)

        if hasattr(mat, 'inv'):
            try:
                return mat.inv()
            except:
                pass

        if mat.is_Inverse:
            return mat.arg

        if not mat.is_square:
            raise ShapeError("Inverse of non-square matrix %s"%mat)

        if mat.is_Mul:
            try:
                return MatMul(*[Inverse(arg) for arg in mat.args[::-1]])
            except ShapeError:
                pass

        return MatPow.__new__(cls, mat, -1)

    @property
    def arg(self):
        return self.args[0]

    @property
    def shape(self):
        return self.arg.shape

from matmul import MatMul

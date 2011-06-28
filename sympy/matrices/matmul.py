from matexpr import MatrixExpr, ShapeError, matrixify
from sympy.core import Mul

class MatMul(MatrixExpr, Mul):

    def __new__(cls, *args):

        # Check that the shape of the args is consistent
        matrices = [arg for arg in args if arg.is_Matrix]

        for i in range(len(matrices)-1):
            A,B = matrices[i:i+2]
            if A.m != B.n:
                raise ShapeError("Matrices %s and %s are not aligned"%(A, B))

        expr = Mul.__new__(cls, *args)
        return matrixify(expr) # Ensure this is a matrix if it should be

    @property
    def shape(self):
        matrices = [arg for arg in self.args if arg.is_Matrix]
        return (matrices[0].n, matrices[-1].m)




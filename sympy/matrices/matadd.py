from matexpr import MatrixExpr, ShapeError, matrixify
from sympy.core import Add


class MatAdd(MatrixExpr, Add):

    def __new__(cls, *args):
        matrices = [arg for arg in args if arg.is_Matrix]
        if len(matrices) != len(args):
            raise ValueError("Mix of Matrix and Scalar symbols")

        # Check that the shape of the args is consistent
        A = matrices[0]
        for B in matrices[1:]:
            if A.shape != B.shape:
                raise ShapeError("Matrices %s and %s are not aligned"%(A,B))
        expr = Add.__new__(cls, *args)
        return matrixify(expr)

    @property
    def shape(self):
        # Return the shape of the first matrix
        for arg in self.args:
            if arg.is_Matrix:
                return arg.shape




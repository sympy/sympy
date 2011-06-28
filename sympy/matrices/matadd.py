from matexpr import MatrixExpr, ShapeError, matrixify
from sympy.core import Add

class MatAdd(MatrixExpr, Add):

    def __new__(cls, *args):
        if not all(arg.is_Matrix for arg in args):
            raise ValueError("Mix of Matrix and Scalar symbols")

        # Check that the shape of the args is consistent
        A = args[0]
        for B in args[1:]:
            if A.shape != B.shape:
                raise ShapeError("Matrices %s and %s are not aligned"%(A,B))

        expr = Add.__new__(cls, *args)
        return matrixify(expr)

    @property
    def shape(self):
        return self.args[0].shape



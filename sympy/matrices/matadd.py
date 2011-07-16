from matexpr import MatrixExpr, ShapeError, matrixify, ZeroMatrix
from sympy import Add, S

class MatAdd(MatrixExpr, Add):

    def __new__(cls, *args):

        args = map(matrixify, args)

        args = [arg for arg in args if arg!=0]

        if not all(arg.is_Matrix for arg in args):
            raise ValueError("Mix of Matrix and Scalar symbols")


        # Check that the shape of the args is consistent
        A = args[0]
        for B in args[1:]:
            if A.shape != B.shape:
                raise ShapeError("Matrices %s and %s are not aligned"%(A,B))

        expr = Add.__new__(cls, *args)
        if expr == S.Zero:
            return ZeroMatrix(*args[0].shape)
        expr = matrixify(expr)

        if expr.is_Mul:
            return MatMul(*expr.args)

        # Clear out Identities
        # Any zeros around?
        if expr.is_Add and any(M.is_ZeroMatrix for M in expr.args):
            newargs = [M for M in expr.args if not M.is_ZeroMatrix] # clear out
            if len(newargs)==0: # Did we lose everything?
                return ZeroMatrix(*args[0].shape)
            if expr.args != newargs: # Removed some 0's but not everything?
                return MatAdd(*newargs) # Repeat with simpler expr

        return expr

    @property
    def shape(self):
        return self.args[0].shape

    def _check_shape(self):
        return (all(arg.shape == self.args[0].shape for arg in self.args) and
                all(arg._check_shape() for arg in self.args))


from matmul import MatMul

from matexpr import MatrixExpr, ShapeError, matrixify, Identity
from sympy.core import Mul

class MatMul(MatrixExpr, Mul):

    def __new__(cls, *args):

        # Check that the shape of the args is consistent
        matrices = [arg for arg in args if arg.is_Matrix]

        for i in range(len(matrices)-1):
            A,B = matrices[i:i+2]
            if A.m != B.n:
                raise ShapeError("Matrices %s and %s are not aligned"%(A, B))

        expr = matrixify(Mul.__new__(cls, *args))

        # Clear out Identities
        nonmats = [M for M in expr.args if not M.is_Matrix] # scalars
        mats = [M for M in expr.args if M.is_Matrix] # matrices
        if any(M.is_Identity for M in mats): # Any identities around?
            newmats = [M for M in mats if not M.is_Identity] # clear out
            if len(newmats)==0: # Did we lose everything?
                newmats = [Identity(expr.n)] # put just one back in

            if mats != newmats: # Removed some I's but not everything?
                return MatMul(*(nonmats+newmats)) # Repeat with simpler expr

        return expr

    @property
    def shape(self):
        matrices = [arg for arg in self.args if arg.is_Matrix]
        return (matrices[0].n, matrices[-1].m)

    def _check_shape(self):
        matrices = [arg for arg in self.args if arg.is_Matrix]
        for A, B in zip(matrices[:-1], matrices[1:]):
            if A.m != B.n:
                return False
        return all(mat._check_shape() for mat in matrices)

from inverse import Inverse

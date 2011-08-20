from matexpr import MatrixExpr, ShapeError, matrixify, Identity, ZeroMatrix
from sympy.core import Mul

class MatMul(MatrixExpr, Mul):
    """A Product of Matrix Expressions

    MatMul inherits from and operates like SymPy Mul

    >>> from sympy import MatMul, MatrixSymbol
    >>> A = MatrixSymbol('A', 5, 4)
    >>> B = MatrixSymbol('B', 4, 3)
    >>> C = MatrixSymbol('C', 3, 6)
    >>> MatMul(A, B, C)
    A*B*C
    """

    def __new__(cls, *args):

        # Check that the shape of the args is consistent
        matrices = [arg for arg in args if arg.is_Matrix]

        for i in range(len(matrices)-1):
            A,B = matrices[i:i+2]
            if A.m != B.n:
                raise ShapeError("Matrices %s and %s are not aligned"%(A, B))

        if any(arg.is_zero for arg in args):
            return ZeroMatrix(matrices[0].n, matrices[-1].m)

        expr = matrixify(Mul.__new__(cls, *args))
        if expr.is_Add:
            return MatAdd(*expr.args)
        if expr.is_Pow:
            return MatPow(*expr.args)
        if not expr.is_Mul:
            return expr

        if any(arg.is_Matrix and arg.is_ZeroMatrix for arg in expr.args):
            return ZeroMatrix(*expr.shape)

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

from matadd import MatAdd
from matpow import MatPow
from inverse import Inverse

from matexpr import MatrixExpr, ShapeError, matrixify, Identity, ZeroMatrix
from sympy.core import Mul, Add, Basic

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
            if A.cols != B.rows:
                raise ShapeError("Matrices %s and %s are not aligned"%(A, B))

        if any(arg.is_zero for arg in args):
            return ZeroMatrix(matrices[0].rows, matrices[-1].cols)

        expr = matrixify(Mul.__new__(cls, *args))
        if expr.is_Add:
            return MatAdd(*expr.args)
        if expr.is_Pow:
            assert expr.exp.is_Integer
            expr = Basic.__new__(MatMul, *[expr.base for i in range(expr.exp)])
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
                newmats = [Identity(expr.rows)] # put just one back in

            if mats != newmats: # Removed some I's but not everything?
                return MatMul(*(nonmats+newmats)) # Repeat with simpler expr

        return expr

    @property
    def shape(self):
        matrices = [arg for arg in self.args if arg.is_Matrix]
        return (matrices[0].rows, matrices[-1].cols)

    def _entry(self, i, j):
        coeff, matmul = self.as_coeff_mmul()
        if not matmul.is_Mul: # situation like 2*X, matmul is just X
            return coeff * matmul[i,j]

        head, tail = matmul.args[0], matmul.args[1:]
        assert len(tail) != 0

        X = head
        Y = MatMul(*tail)

        if X.shape[1].is_Number:
            # Numeric shape like (3,5)
            return coeff*Add(*[X[i,k]*Y[k,j] for k in range(X.shape[1])])
        else:
            # Symbolic shape like (n, m)
            from sympy import Dummy, summation
            k = Dummy('k', integer=True)
            return summation(coeff*X[i,k]*Y[k,j], (k, 0, X.cols-1))

    def as_coeff_mmul(self):
        scalars = [x for x in self.args if not x.is_Matrix]
        matrices = [x for x in self.args if x.is_Matrix]
        coeff = Mul(*scalars)

        return coeff, MatMul(*matrices)

    def _eval_transpose(self):
        from transpose import Transpose
        return MatMul(*[Transpose(arg) for arg in self.args[::-1]])

    def _eval_trace(self):
        factor = Mul(*[arg for arg in self.args if not arg.is_Matrix])
        matrix = MatMul(*[arg for arg in self.args if arg.is_Matrix])
        if factor != 1:
            from trace import Trace
            return factor * Trace(matrix)
        else:
            raise NotImplementedError("Can't simplify any further")

    def _eval_inverse(self):
        from inverse import Inverse
        try:
            return MatMul(*[Inverse(arg) for arg in self.args[::-1]])
        except ShapeError:
            raise NotImplementedError("Can not decompose this Inverse")


from matadd import MatAdd

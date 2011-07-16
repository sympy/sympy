from matexpr import MatrixExpr
from matmul import MatMul
from matadd import MatAdd
from matpow import MatPow
from transpose import Transpose
from matrices import Matrix
from sympy import Tuple, Basic

class BlockMatrix(MatrixExpr):
    is_BlockMatrix = True
    def __new__(cls, mat):
        if not isinstance(mat, Matrix):
            mat = Matrix(mat)
        data = Tuple(*mat.mat)
        if len(data) == 1:
            return mat[0,0]
        obj = Basic.__new__(cls, data, Tuple(*mat.shape))
        obj.mat = mat
        return obj

    @property
    def shape(self):
        numrows = numcols = 0
        M = self.mat
        for i in range(M.shape[0]):
            numrows += M[i,0].shape[0]
        for i in range(M.shape[1]):
            numcols += M[0,i].shape[1]
        return (numrows, numcols)

    @property
    def blockshape(self):
        return self.mat.shape

    @property
    def rowblocksizes(self):
        return [self.mat[i,0].n for i in range(self.blockshape[0])]

    @property
    def colblocksizes(self):
        return [self.mat[0,i].m for i in range(self.blockshape[1])]

    def _blockmul(self, other):

        if  (other.is_Matrix and other.is_BlockMatrix and
                self.blockshape[1] == other.blockshape[0] and
                self.colblocksizes == other.rowblocksizes):
            return BlockMatrix(self.mat*other.mat)

        return MatrixExpr.__mul__(self, other)

    def _blockadd(self, other):

        if  (other.is_Matrix and other.is_BlockMatrix and
                self.blockshape == other.blockshape and
                self.rowblocksizes == other.rowblocksizes and
                self.colblocksizes == other.colblocksizes):
            return BlockMatrix(self.mat + other.mat)

        return MatrixExpr.__add__(self, other)

    def eval_transpose(self):
        # Flip all the individual matrices
        matrices = [Transpose(matrix) for matrix in self.mat.mat]
        # Flip the block shape
        return BlockMatrix(
                Matrix(self.blockshape[1], self.blockshape[0], matrices))

    def transpose(self):
        return self.eval_transpose()

    def __getitem__(self, *args):
        return self.mat.__getitem__(*args)

def block_collapse(expr):
    if not expr.is_Matrix or (not expr.is_Add and not expr.is_Mul
            and not expr.is_Transpose and not expr.is_Pow):
        return expr

    if expr.is_Transpose:
        return Transpose(block_collapse(expr.arg))

    # Recurse on the subargs
    newargs = map(block_collapse, expr.args)
    if tuple(newargs) != expr.args:
        expr = expr.__class__(*newargs)

    # Turn  -[X, Y] into [-X, -Y]
    if (expr.is_Mul and len(expr.args)==2 and expr.args[0]==-1
            and expr.args[1].is_BlockMatrix):
        return BlockMatrix((-1)*expr.args[1].mat)

    if expr.is_Add:
        nonblocks = [arg for arg in expr.args if not arg.is_BlockMatrix]
        blocks = [arg for arg in expr.args if arg.is_BlockMatrix]
        block = reduce(lambda a,b: a._blockadd(b), blocks[1:], blocks[0])
        return MatAdd(*(nonblocks+[block]))

    if expr.is_Mul:
        nonmatrices = [arg for arg in expr.args if not arg.is_Matrix]
        matrices = [arg for arg in expr.args if arg.is_Matrix]
        newmatrices = []
        i = 0
        while (i+1<len(matrices)):
            A, B = matrices[i:i+2]
            if A.is_BlockMatrix and B.is_BlockMatrix:
                matrices[i] = A._blockmul(B)
                matrices.pop(i+1)
            else:
                i+=1
        return MatMul(*(nonmatrices+matrices))

    if expr.is_Pow:
        rv = expr.base
        for i in range(1, expr.exp):
            rv = rv._blockmul(expr.base)
        return rv




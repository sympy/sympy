from matexpr import MatrixExpr
from matrices import Matrix
from sympy import Tuple, Basic

class BlockMatrix(MatrixExpr):
    is_BlockMatrix = True
    def __new__(cls, mat):
        if not isinstance(mat, Matrix):
            mat = Matrix(mat)
        data = Tuple(*mat.mat)
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

    def __mul__(self, other):

        if  (other.is_Matrix and other.is_BlockMatrix and
                self.blockshape[1] == other.blockshape[0] and
                self.colblocksizes == other.rowblocksizes):
            return BlockMatrix(self.mat*other.mat)

        return MatrixExpr.__mul__(self, other)

    def __add__(self, other):

        if  (other.is_Matrix and other.is_BlockMatrix and
                self.blockshape == other.blockshape and
                self.rowblocksizes == other.rowblocksizes and
                self.colblocksizes == other.colblocksizes):
            return BlockMatrix(self.mat + other.mat)

        return MatrixExpr.__add__(self, other)



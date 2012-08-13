from matexpr import MatrixExpr, ZeroMatrix, Identity
from matmul import MatMul
from matadd import MatAdd
from transpose import Transpose
from trace import Trace
from inverse import Inverse
from sympy.matrices import Matrix, eye
from sympy import Tuple, Basic, sympify, FiniteSet, Add

class BlockMatrix(MatrixExpr):
    """A BlockMatrix is a Matrix composed of other smaller, submatrices

    The submatrices are stored in a SymPy Matrix object but accessed as part of
    a Matrix Expression

    >>> from sympy import MatrixSymbol, BlockMatrix, symbols, Identity, ZeroMatrix, block_collapse
    >>> n,m,l = symbols('n m l')
    >>> X = MatrixSymbol('X', n, n)
    >>> Y = MatrixSymbol('Y', m ,m)
    >>> Z = MatrixSymbol('Z', n, m)
    >>> B = BlockMatrix([[X, Z], [ZeroMatrix(m,n), Y]])
    >>> print B
    [X, Z]
    [0, Y]

    >>> C = BlockMatrix([[Identity(n), Z]])
    >>> print C
    [I, Z]

    >>> print block_collapse(C*B)
    [X, Z + Z*Y]

    """
    is_BlockMatrix = True
    is_BlockDiagMatrix = False
    def __new__(cls, mat):
        if not isinstance(mat, Matrix):
            mat = Matrix(mat)
        data = Tuple(*mat.mat)
        shape = Tuple(*sympify(mat.shape))
        obj = Basic.__new__(cls, data, shape)
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
    def blocks(self):
        return self.mat

    @property
    def rowblocksizes(self):
        return [self.blocks[i,0].rows for i in range(self.blockshape[0])]

    @property
    def colblocksizes(self):
        return [self.blocks[0,i].cols for i in range(self.blockshape[1])]

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

    def _eval_transpose(self):
        # Flip all the individual matrices
        matrices = [Transpose(matrix) for matrix in self.mat.mat]
        # Make a copy
        mat = Matrix(self.blockshape[0], self.blockshape[1], matrices)
        # Transpose the block structure
        mat = mat.transpose()
        return BlockMatrix(mat)

    def _eval_trace(self):
        if self.rowblocksizes == self.colblocksizes:
            return Add(*[Trace(self.mat[i,i])
                        for i in range(self.blockshape[0])])
        raise NotImplementedError("Can't perform trace of irregular blockshape")

    #def transpose(self):
    #    return self.eval_transpose()

    def _eval_inverse(self, expand=False):
        # Inverse of one by one block matrix is easy
        if self.blockshape==(1,1):
            mat = Matrix(1, 1, (Inverse(self.blocks[0]), ))
            return BlockMatrix(mat)
        # Inverse of a two by two block matrix is known
        elif expand and self.blockshape==(2,2):
            # Cite: The Matrix Cookbook Section 9.1.3
            A11, A12, A21, A22 = (self.blocks[0,0], self.blocks[0,1],
                    self.blocks[1,0], self.blocks[1,1])
            C1 = A11 - A12*Inverse(A22)*A21
            C2 = A22 - A21*Inverse(A11)*A12
            mat = Matrix([[Inverse(C1), Inverse(-A11)*A12*Inverse(C2)],
                [-Inverse(C2)*A21*Inverse(A11), Inverse(C2)]])
            return BlockMatrix(mat)
        else:
            raise NotImplementedError()

    def inverse(self):
        return Inverse(self)

    def _entry(self, i, j):
        idx = 0
        # Find row entry
        for row_block, numrows in enumerate(self.rowblocksizes):
            if i < numrows:
                break
            else:
                i -= numrows
        for col_block, numcols in enumerate(self.colblocksizes):
            if j < numcols:
                break
            else:
                j -= numcols
        return self.blocks[row_block, col_block][i,j]

    @property
    def is_Identity(self):
        if self.blockshape[0] != self.blockshape[1]:
            return False
        for i in range(self.blockshape[0]):
            for j in range(self.blockshape[1]):
                if i==j and not self.mat[i,j].is_Identity:
                    return False
                if i!=j and not self.mat[i,j].is_ZeroMatrix:
                    return False
        return True
    @property
    def is_structurally_symmetric(self):
        return self.rowblocksizes == self.colblocksizes

class BlockDiagMatrix(BlockMatrix):
    """A BlockDiagMatrix is a BlockMatrix with matrices only along the diagonal

    >>> from sympy import MatrixSymbol, BlockDiagMatrix, symbols, Identity
    >>> n,m,l = symbols('n m l')
    >>> X = MatrixSymbol('X', n, n)
    >>> Y = MatrixSymbol('Y', m ,m)
    >>> BlockDiagMatrix(X, Y)
    [X, 0]
    [0, Y]

    """
    is_BlockDiagMatrix = True
    def __new__(cls, *mats):
        data_matrix = eye(len(mats))
        for i, mat in enumerate(mats):
            data_matrix[i,i] = mat

        for r in range(len(mats)):
            for c in range(len(mats)):
                if r == c:
                    continue
                n = mats[r].rows
                m = mats[c].cols
                data_matrix[r, c] = ZeroMatrix(n, m)

        shape = Tuple(*sympify(mat.shape))
        data = Tuple(*data_matrix.mat)
        obj = Basic.__new__(cls, data, shape, Tuple(*mats))
        obj.mat = data_matrix
        return obj

    @property
    def diag(self):
        return self.args[2]

    def _eval_inverse(self):
        return BlockDiagMatrix(*[Inverse(mat) for mat in self.diag])

    def _blockmul(self, other):
        if  (other.is_Matrix and other.is_BlockMatrix and
                other.is_BlockDiagMatrix and
                self.blockshape[1] == other.blockshape[0] and
                self.colblocksizes == other.rowblocksizes):
            return BlockDiagMatrix(*[a*b for a, b in zip(self.diag,other.diag)])
        else:
            return BlockMatrix._blockmul(self, other)

    def _blockadd(self, other):

        if  (other.is_Matrix and other.is_BlockMatrix and
                other.is_BlockDiagMatrix and
                self.blockshape == other.blockshape and
                self.rowblocksizes == other.rowblocksizes and
                self.colblocksizes == other.colblocksizes):
            return BlockDiagMatrix(*[a+b for a, b in zip(self.diag,other.diag)])
        else:
            return BlockMatrix._blockadd(self, other)

def block_collapse(expr):
    """Evaluates a block matrix expression

    >>> from sympy import MatrixSymbol, BlockMatrix, symbols, Identity, Matrix, ZeroMatrix, block_collapse
    >>> n,m,l = symbols('n m l')
    >>> X = MatrixSymbol('X', n, n)
    >>> Y = MatrixSymbol('Y', m ,m)
    >>> Z = MatrixSymbol('Z', n, m)
    >>> B = BlockMatrix([[X, Z], [ZeroMatrix(m, n), Y]])
    >>> print B
    [X, Z]
    [0, Y]

    >>> C = BlockMatrix([[Identity(n), Z]])
    >>> print C
    [I, Z]

    >>> print block_collapse(C*B)
    [X, Z + Z*Y]
    """
    if expr.__class__ in [tuple, list, set, frozenset]:
        return expr.__class__([block_collapse(arg) for arg in expr])
    if expr.__class__ in [Tuple, FiniteSet]:
        return expr.__class__(*[block_collapse(arg) for arg in expr])

    if not expr.is_Matrix or (not expr.is_Add and not expr.is_Mul
            and not expr.is_Transpose and not expr.is_Pow
            and not expr.is_Inverse):
        return expr

    if expr.is_Transpose:
        expr = Transpose(block_collapse(expr.arg))
        if expr.is_Transpose and expr.arg.is_BlockMatrix:
            expr = expr.arg._eval_transpose()
        return expr

    if expr.is_Inverse:
        return Inverse(block_collapse(expr.arg))

    # Recurse on the subargs
    args = list(expr.args)
    for i in range(len(args)):
        arg = args[i]
        newarg = block_collapse(arg)
        while(newarg != arg): # Repeat until no new changes
            arg = newarg
            newarg = block_collapse(arg)
        args[i] = newarg

    if tuple(args) != expr.args:
        expr = expr.__class__(*args)

    # Turn  -[X, Y] into [-X, -Y]
    if (expr.is_Mul and len(expr.args)==2 and not expr.args[0].is_Matrix
            and expr.args[1].is_BlockMatrix):
        if expr.args[1].is_BlockDiagMatrix:
            return BlockDiagMatrix(
                    *[expr.args[0]*arg for arg in expr.args[1].diag])
        else:
            return BlockMatrix(expr.args[0]*expr.args[1].mat)

    if expr.is_Add:
        nonblocks = [arg for arg in expr.args if not arg.is_BlockMatrix]
        blocks = [arg for arg in expr.args if arg.is_BlockMatrix]
        if not blocks:
            return MatAdd(*nonblocks)
        block = blocks[0]
        for b in blocks[1:]:
            block = block._blockadd(b)
        if block.blockshape == (1,1):
            # Bring all the non-blocks into the block_matrix
            mat = Matrix(1, 1, (block.blocks[0,0] + MatAdd(*nonblocks), ))
            return BlockMatrix(mat)
        # Add identities to the blocks as block identities
        for i, mat in enumerate(nonblocks):
            c, M = mat.as_coeff_Mul()
            if M.is_Identity and block.is_structurally_symmetric:
                block_id = BlockDiagMatrix(
                        *[c*Identity(k) for k in block.rowblocksizes])
                nonblocks.pop(i)
                block = block._blockadd(block_id)


        return MatAdd(*(nonblocks+[block]))

    if expr.is_Mul:
        nonmatrices = [arg for arg in expr.args if not arg.is_Matrix]
        matrices = [arg for arg in expr.args if arg.is_Matrix]
        i = 0
        while (i+1 < len(matrices)):
            A, B = matrices[i:i+2]
            if A.is_BlockMatrix and B.is_BlockMatrix:
                matrices[i] = A._blockmul(B)
                matrices.pop(i+1)
            else:
                i+=1
        return MatMul(*(nonmatrices + matrices))

    if expr.is_Pow:
        rv = expr.base
        for i in range(1, expr.exp):
            rv = rv._blockmul(expr.base)
        return rv

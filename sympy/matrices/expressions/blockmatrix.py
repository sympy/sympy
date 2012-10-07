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

    >>> from sympy import (MatrixSymbol, BlockMatrix, symbols,
    ...     Identity, ZeroMatrix, block_collapse)
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

    def __new__(cls, *args):
        from sympy.matrices.immutable import ImmutableMatrix
        mat = ImmutableMatrix(*args)

        obj = Basic.__new__(cls, mat)
        return obj


    @property
    def _mat(self):
        return self.args[0]

    @property
    def shape(self):
        numrows = numcols = 0
        M = self._mat
        for i in range(M.shape[0]):
            numrows += M[i, 0].shape[0]
        for i in range(M.shape[1]):
            numcols += M[0, i].shape[1]
        return (numrows, numcols)

    @property
    def blockshape(self):
        return self._mat.shape

    @property
    def blocks(self):
        return self._mat

    @property
    def rowblocksizes(self):
        return [self.blocks[i, 0].rows for i in range(self.blockshape[0])]

    @property
    def colblocksizes(self):
        return [self.blocks[0, i].cols for i in range(self.blockshape[1])]

    def structurally_equal(self, other):
        return (self.shape == other.shape
            and self.blockshape == other.blockshape
            and self.rowblocksizes == other.rowblocksizes
            and self.colblocksizes == other.colblocksizes)

    def _blockmul(self, other):

        if  (other.is_Matrix and other.is_BlockMatrix and
                self.blockshape[1] == other.blockshape[0] and
                self.colblocksizes == other.rowblocksizes):
            return BlockMatrix(self._mat*other._mat)

        return MatrixExpr.__mul__(self, other)

    def _blockadd(self, other):
        if   (other.is_Matrix and other.is_BlockMatrix
                and self.structurally_equal(other)):
            return BlockMatrix(self._mat + other._mat)

        return MatrixExpr.__add__(self, other)

    def _eval_transpose(self):
        # Flip all the individual matrices
        matrices = [Transpose(matrix) for matrix in self._mat]
        # Make a copy
        M = Matrix(self.blockshape[0], self.blockshape[1], matrices)
        # Transpose the block structure
        M = M.transpose()
        return BlockMatrix(M)

    def _eval_trace(self):
        if self.rowblocksizes == self.colblocksizes:
            return Add(*[Trace(self._mat[i, i])
                        for i in range(self.blockshape[0])])
        raise NotImplementedError(
            "Can't perform trace of irregular blockshape")

    def transpose(self):
        """Return transpose of matrix.

        Examples
        ========

        >>> from sympy import MatrixSymbol, BlockMatrix, ZeroMatrix
        >>> from sympy.abc import l, m, n
        >>> X = MatrixSymbol('X', n, n)
        >>> Y = MatrixSymbol('Y', m ,m)
        >>> Z = MatrixSymbol('Z', n, m)
        >>> B = BlockMatrix([[X, Z], [ZeroMatrix(m,n), Y]])
        >>> B.transpose()
        [X',  0]
        [Z', Y']
        >>> _.transpose()
        [X, Z]
        [0, Y]
        """
        return self._eval_transpose()

    def _eval_inverse(self, expand=False):
        # Inverse of one by one block matrix is easy
        if self.blockshape == (1, 1):
            mat = Matrix(1, 1, (Inverse(self.blocks[0]), ))
            return BlockMatrix(mat)
        # Inverse of a two by two block matrix is known
        elif expand and self.blockshape == (2, 2):
            # Cite: The Matrix Cookbook Section 9.1.3
            A11, A12, A21, A22 = (self.blocks[0, 0], self.blocks[0, 1],
                    self.blocks[1, 0], self.blocks[1, 1])
            C1 = A11 - A12*Inverse(A22)*A21
            C2 = A22 - A21*Inverse(A11)*A12
            mat = Matrix([[Inverse(C1), Inverse(-A11)*A12*Inverse(C2)],
                [-Inverse(C2)*A21*Inverse(A11), Inverse(C2)]])
            return BlockMatrix(mat)
        else:
            raise NotImplementedError()

    def inv(self, expand=False):
        """Return inverse of matrix.

        Examples
        ========

        >>> from sympy import MatrixSymbol, BlockMatrix, ZeroMatrix
        >>> from sympy.abc import l, m, n
        >>> X = MatrixSymbol('X', n, n)
        >>> BlockMatrix([[X]]).inv()
        [X^-1]

        >>> Y = MatrixSymbol('Y', m ,m)
        >>> Z = MatrixSymbol('Z', n, m)
        >>> B = BlockMatrix([[X, Z], [ZeroMatrix(m,n), Y]])
        >>> B
        [X, Z]
        [0, Y]
        >>> B.inv(expand=True)
        [X^-1, (-1)*X^-1*Z*Y^-1]
        [   0,             Y^-1]

        """
        return self._eval_inverse(expand)

    def inverse(self):
        # XXX document how this is different than inv
        return Inverse(self)

    def _entry(self, i, j):
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
        return self.blocks[row_block, col_block][i, j]

    @property
    def is_Identity(self):
        if self.blockshape[0] != self.blockshape[1]:
            return False
        for i in range(self.blockshape[0]):
            for j in range(self.blockshape[1]):
                if i == j and not self._mat[i, j].is_Identity:
                    return False
                if i != j and not self._mat[i, j].is_ZeroMatrix:
                    return False
        return True

    @property
    def is_structurally_symmetric(self):
        return self.rowblocksizes == self.colblocksizes

    def equals(self, other):
        if self == other:
            return True
        if (self.is_BlockMatrix and other.is_BlockMatrix and
                self._mat == other._mat):
            return True
        return super(BlockMatrix, self).equals(other)

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
        from sympy.matrices.immutable import ImmutableMatrix
        data = [[mats[i] if i == j else ZeroMatrix(mats[i].rows, mats[j].cols)
                        for j in range(len(mats))]
                        for i in range(len(mats))]
        return Basic.__new__(BlockDiagMatrix, ImmutableMatrix(data))

    @property
    def diag(self):
        return [self.blocks[i,i] for i in range(self.blocks.rows)]

    def _eval_inverse(self):
        return BlockDiagMatrix(*[Inverse(mat) for mat in self.diag])

    def _blockmul(self, other):
        if  (other.is_Matrix and other.is_BlockMatrix and
                other.is_BlockDiagMatrix and
                self.blockshape[1] == other.blockshape[0] and
                self.colblocksizes == other.rowblocksizes):
            return BlockDiagMatrix(*[a*b for a, b in zip(self.diag, other.diag)])
        else:
            return BlockMatrix._blockmul(self, other)

    def _blockadd(self, other):

        if  (other.is_Matrix and other.is_BlockMatrix and
                other.is_BlockDiagMatrix and
                self.blockshape == other.blockshape and
                self.rowblocksizes == other.rowblocksizes and
                self.colblocksizes == other.colblocksizes):
            return BlockDiagMatrix(*[a + b for a, b in zip(self.diag, other.diag)])
        else:
            return BlockMatrix._blockadd(self, other)

from sympy.rr import distribute, typed, canon, chain, try_safe, debug

def block_collapse(expr):
    """Evaluates a block matrix expression

    >>> from sympy import MatrixSymbol, BlockMatrix, symbols, \
                          Identity, Matrix, ZeroMatrix, block_collapse
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
    from sympy import MatAdd, MatMul, MatPow
    rule = canon(typed({MatAdd: chain(bc_matadd, try_safe(bc_block_plus_ident)),
                        MatMul: chain(bc_matmul, try_safe(bc_dist)),
                        MatPow: bc_matpow,
                        BlockMatrix: bc_unpack}))
    result = rule(expr)
    try:                    return result.canonicalize()
    except AttributeError:  return result

#[BlockMatrix]
def bc_unpack(expr):
    if expr.blockshape == (1, 1):
        return expr.blocks[0,0]
    return expr

#[MatAdd]
def bc_matadd(expr):
    nonblocks = [arg for arg in expr.args if not arg.is_BlockMatrix]
    blocks = [arg for arg in expr.args if arg.is_BlockMatrix]
    if not blocks:
        return expr

    block = blocks[0]
    for b in blocks[1:]:
        block = block._blockadd(b)

    return MatAdd(*(nonblocks+[block]))

def bc_block_plus_ident(expr):
    idents = [arg for arg in expr.args if arg.is_Identity]
    if not idents:
        return expr

    blocks = [arg for arg in expr.args if arg.is_BlockMatrix]
    if (blocks and all(b.structurally_equal(blocks[0]) for b in blocks)
               and blocks[0].is_structurally_symmetric):
        block_id = BlockDiagMatrix(*[Identity(k)
                                        for k in blocks[0].rowblocksizes])
        return MatAdd(block_id * len(idents), *blocks)

    return expr

#[MatMul]
def bc_dist(expr):
    """ Turn  -[X, Y] into [-X, -Y] """
    factor, mat = expr.as_coeff_mmul()
    if factor != 1 and mat.is_BlockMatrix:
        B = mat.blocks
        return BlockMatrix([[factor * B[i,j] for j in range(B.cols)]
                                             for i in range(B.rows)])
    return expr

def bc_matmul(expr):
    factor, matrices = expr.as_coeff_matrices()
    i = 0
    while (i+1 < len(matrices)):
        A, B = matrices[i:i+2]
        if A.is_BlockMatrix and B.is_BlockMatrix:
            matrices[i] = A._blockmul(B)
            matrices.pop(i+1)
        else:
            i+=1
    return MatMul(factor, *matrices)

#[MatPow]
def bc_matpow(expr):
    return MatMul(*([expr.base]*expr.exp))

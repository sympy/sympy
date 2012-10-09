from sympy.matrices.expressions.blockmatrix import (bc_matmul,
        bc_block_plus_ident, BlockDiagMatrix, BlockMatrix, bc_dist, bc_matadd)
from sympy.matrices.expressions import MatrixSymbol, Identity
from sympy import symbols, Matrix

n,m,l,k,o,p = symbols('n,m,l,k,o,p')
A = MatrixSymbol('A', n, m)
B = MatrixSymbol('B', n, k)
C = MatrixSymbol('C', l, m)
D = MatrixSymbol('D', l, k)
M = MatrixSymbol('M', m+k, p)
N = MatrixSymbol('N', l+n, k+m)
X = BlockMatrix(Matrix([[A,B],[C,D]]))

G = MatrixSymbol('G', n, n)
H = MatrixSymbol('H', n, n)
b1 = BlockMatrix([[G, H]])
b2 = BlockMatrix([[G], [H]])

def test_bc_matmul():
    assert bc_matmul(H*b1*b2*G) == H*BlockMatrix([[G*G + H*H]])*G

#def test_bc_transpose():
#    assert bc_transpose(BlockMatrix([[G, H]]).T) == BlockMatrix([[G.T, H.T]])

def test_bc_matadd():
    assert bc_matadd(BlockMatrix([[G, H]]) + BlockMatrix([[H, H]])) == \
            BlockMatrix([[G+H, H+H]])

def test_bc_dist_diag():
    A = MatrixSymbol('A', n, n)
    B = MatrixSymbol('B', m, m)
    C = MatrixSymbol('C', l, l)
    X = BlockDiagMatrix(A,B,C)

    assert bc_dist(X+X).equals(BlockDiagMatrix(2*A, 2*B, 2*C))

def test_block_plus_ident():
    A = MatrixSymbol('A', n, n)
    B = MatrixSymbol('B', n, m)
    C = MatrixSymbol('C', m, n)
    D = MatrixSymbol('D', m, m)
    X = BlockMatrix([[A,B],[C,D]])
    assert bc_block_plus_ident(X+Identity(m+n)) == \
            BlockDiagMatrix(Identity(n), Identity(m)) + X

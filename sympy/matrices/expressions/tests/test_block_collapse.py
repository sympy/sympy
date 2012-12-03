from sympy.matrices.expressions.blockmatrix import (bc_matmul,
        bc_block_plus_ident, BlockDiagMatrix, BlockMatrix, bc_dist, bc_matadd,
        bc_matpow)
from sympy.matrices.expressions import MatrixSymbol, Identity, MatMul
from sympy import symbols, Matrix

l, m, n = symbols('l, m, n', integer=True)
G = MatrixSymbol('G', n, n)
H = MatrixSymbol('H', n, n)
b1 = BlockMatrix([[G, H]])
b2 = BlockMatrix([[G], [H]])

def test_bc_matmul():
    assert bc_matmul(H*b1*b2*G) == H*BlockMatrix([[G*G + H*H]])*G

def test_bc_matadd():
    assert bc_matadd(BlockMatrix([[G, H]]) + BlockMatrix([[H, H]])) == \
            BlockMatrix([[G+H, H+H]])

def test_bc_dist_diag():
    A = MatrixSymbol('A', n, n)
    B = MatrixSymbol('B', m, m)
    C = MatrixSymbol('C', l, l)
    X = BlockDiagMatrix(A, B, C)

    assert bc_dist(X+X).equals(BlockDiagMatrix(2*A, 2*B, 2*C))

def test_bc_matpow():
    A = MatrixSymbol('A', n, n)
    assert bc_matpow(A**2) == MatMul(A, A)
    assert bc_matpow(A**n) == A**n

def test_block_plus_ident():
    A = MatrixSymbol('A', n, n)
    B = MatrixSymbol('B', n, m)
    C = MatrixSymbol('C', m, n)
    D = MatrixSymbol('D', m, m)
    X = BlockMatrix([[A, B], [C, D]])
    assert bc_block_plus_ident(X+Identity(m+n)) == \
            BlockDiagMatrix(Identity(n), Identity(m)) + X

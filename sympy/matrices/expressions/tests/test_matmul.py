from sympy.matrices.expressions.matmul import (factor_in_front, remove_ids,
        MatMul, xxinv, any_zeros, unpack)
from sympy import symbols, MatrixSymbol, Identity, Inverse, ZeroMatrix

n, m, l, k = symbols('n m l k', integer=True)
A = MatrixSymbol('A', n, m)
B = MatrixSymbol('B', m, l)
C = MatrixSymbol('C', n, n)
D = MatrixSymbol('D', n, n)
E = MatrixSymbol('E', m, n)


def test_factor_in_front():
    assert factor_in_front(MatMul(A, 2, B, simplify=False)) ==\
                           MatMul(2, A, B, simplify=False)

def test_remove_ids():
    assert remove_ids(MatMul(A, Identity(m), B, simplify=False)) == \
                      MatMul(A, B, simplify=False)
    assert remove_ids(MatMul(Identity(n), simplify=False)) == \
                      MatMul(Identity(n), simplify=False)

def test_xxinv():
    assert xxinv(MatMul(D, Inverse(D), D, simplify=False)) == \
                 MatMul(Identity(n), D, simplify=False)

def test_any_zeros():
    assert any_zeros(MatMul(A, ZeroMatrix(m, k), simplify=False)) == \
                     ZeroMatrix(n, k)

def test_unpack():
    assert unpack(MatMul(A, simplify=False)) == A
    x = MatMul(A, B)
    assert unpack(x) == x

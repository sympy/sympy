from sympy.matrices.expressions.matmul import (factor_in_front, remove_ids,
        MatMul, xxinv, any_zeros, unpack)
from sympy import symbols, MatrixSymbol, Identity, Inverse, ZeroMatrix
from sympy.rules import null_safe

n, m, l, k = symbols('n m l k', integer=True)
A = MatrixSymbol('A', n, m)
B = MatrixSymbol('B', m, l)
C = MatrixSymbol('C', n, n)
D = MatrixSymbol('D', n, n)
E = MatrixSymbol('E', m, n)


def test_factor_in_front():
    assert factor_in_front(MatMul(A, 2, B, evaluate=False)) ==\
                           MatMul(2, A, B, evaluate=False)


def test_remove_ids():
    assert remove_ids(MatMul(A, Identity(m), B, evaluate=False)) == \
                      MatMul(A, B, evaluate=False)
    assert null_safe(remove_ids)(MatMul(Identity(n), evaluate=False)) == \
                                 MatMul(Identity(n), evaluate=False)


def test_xxinv():
    assert xxinv(MatMul(D, Inverse(D), D, evaluate=False)) == \
                 MatMul(Identity(n), D, evaluate=False)


def test_any_zeros():
    assert any_zeros(MatMul(A, ZeroMatrix(m, k), evaluate=False)) == \
                     ZeroMatrix(n, k)


def test_unpack():
    assert unpack(MatMul(A, evaluate=False)) == A
    x = MatMul(A, B)
    assert unpack(x) == x

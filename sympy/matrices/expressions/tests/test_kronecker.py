from sympy import I, symbols, Matrix
from sympy.matrices import ShapeError, MatrixSymbol
from sympy.matrices.expressions import det, trace

from sympy.matrices.expressions.kronecker import KroneckerProduct
from sympy.matrices.expressions.kronecker import kronecker_product
from sympy.core.trace import Tr

mat1 = Matrix([[1, 2 * I], [1 + I, 3]])
mat2 = Matrix([[2 * I, 3], [4 * I, 2]])

n, m, k, x = symbols('n,m,k,x')
Z = MatrixSymbol('Z', n, n)
W = MatrixSymbol('W', m, m)
A = MatrixSymbol('A', n, m)
B = MatrixSymbol('B', n, m)
C = MatrixSymbol('C', m, k)


def test_KroneckerProduct():
    assert isinstance(KroneckerProduct(A, B), KroneckerProduct)
    assert KroneckerProduct(A, B).subs(A, C) == KroneckerProduct(C, B)
    assert KroneckerProduct(A, C).shape == (n*m, m*k)


def test_KroneckerProduct_explicit():
    X = MatrixSymbol('X', 2, 2)
    Y = MatrixSymbol('Y', 2, 2)
    kp = KroneckerProduct(X, Y)
    assert kp.shape == (4, 4)
    assert kp.as_explicit() == Matrix(
        [
            [X[0, 0]*Y[0, 0], X[0, 0]*Y[0, 1], X[0, 1]*Y[0, 0], X[0, 1]*Y[0, 1]],
            [X[0, 0]*Y[1, 0], X[0, 0]*Y[1, 1], X[0, 1]*Y[1, 0], X[0, 1]*Y[1, 1]],
            [X[1, 0]*Y[0, 0], X[1, 0]*Y[0, 1], X[1, 1]*Y[0, 0], X[1, 1]*Y[0, 1]],
            [X[1, 0]*Y[1, 0], X[1, 0]*Y[1, 1], X[1, 1]*Y[1, 0], X[1, 1]*Y[1, 1]]
        ]
    )


def test_tensor_product_adjoint():
    assert KroneckerProduct(I*A, B).adjoint() == \
        -I*KroneckerProduct(A.adjoint(), B.adjoint())
    assert KroneckerProduct(mat1, mat2).adjoint() == \
        kronecker_product(mat1.adjoint(), mat2.adjoint())


def test_tensor_product_conjugate():
    assert KroneckerProduct(I*A, B).conjugate() == \
        -I*KroneckerProduct(A.conjugate(), B.conjugate())
    assert KroneckerProduct(mat1, mat2).conjugate() == \
        kronecker_product(mat1.conjugate(), mat2.conjugate())


def test_tensor_product_transpose():
    assert KroneckerProduct(I*A, B).transpose() == \
        I*KroneckerProduct(A.transpose(), B.transpose())
    assert KroneckerProduct(mat1, mat2).transpose() == \
        kronecker_product(mat1.transpose(), mat2.transpose())


def test_KroneckerProduct_is_associative():
    assert kronecker_product(A, kronecker_product(
        B, C)) == kronecker_product(kronecker_product(A, B), C)
    assert kronecker_product(A, kronecker_product(
        B, C)) == KroneckerProduct(A, B, C)


def test_KroneckerProduct_is_bilinear():
    assert kronecker_product(x*A, B) == x*kronecker_product(A, B)
    assert kronecker_product(A, x*B) == x*kronecker_product(A, B)


def test_KroneckerProduct_determinant():
    kp = kronecker_product(W, Z)
    assert det(kp) == det(W)**n * det(Z)**m


def test_KroneckerProduct_trace():
    kp = kronecker_product(W, Z)
    assert trace(kp) == trace(W)*trace(Z)


def test_KroneckerProduct_isnt_commutative():
    assert KroneckerProduct(A, B) != KroneckerProduct(B, A)
    assert KroneckerProduct(A, B).is_commutative is False


def test_KroneckerProduct_extracts_commutative_part():
    assert kronecker_product(x * A, 2 * B) == x * \
        2 * KroneckerProduct(A, B)


def test_KroneckerProduct_inverse():
    kp = kronecker_product(W, Z)
    assert kp.inverse() == kronecker_product(W.inverse(), Z.inverse())


# def test_kronecker_product_expand():
#     assert KroneckerProduct(A + B, B + C).expand(kroneckerproduct=True) == \
#         KroneckerProduct(A, B) + KroneckerProduct(A, C) + KroneckerProduct(B, B) + KroneckerProduct(B, C)


# def test_kronecker_product_simp():
#     assert kronecker_product_simp(KroneckerProduct(A, B) * KroneckerProduct(B, C)) == KroneckerProduct(A * B, B * C)
#     # tests for Pow-expressions
#     assert kronecker_product_simp(KroneckerProduct(A, B)**x) == KroneckerProduct(A**x, B**x)
#     assert kronecker_product_simp(x * KroneckerProduct(A, B)**2) == x * KroneckerProduct(A**2, B**2)
#     assert kronecker_product_simp(
#         x * (KroneckerProduct(A, B)**2) * KroneckerProduct(C, D)) == x * KroneckerProduct(A**2 * C, B**2 * D)
#     assert kronecker_product_simp(KroneckerProduct(A, B) - KroneckerProduct(C, D)**x) == KroneckerProduct(A, B) - KroneckerProduct(
#         C**x, D**x)

from sympy import PurePoly, Symbol
from sympy.matrices import (Matrix, ImmutableDenseMatrix, SparseMatrix,
    ImmutableSparseMatrix)

from sympy import det
from sympy.matrices.determinant import (
    adjugate, charpoly, cofactor, cofactor_matrix, det_bareiss, det_berkowitz,
    det_LU, minor, minor_submatrix)

x = Symbol('x')


# determinant.py

def test_adjugate():
    M = Matrix(4, 4, [0, 1, 2, 3, 4, 5, 6, 7, 2, 9, 10, 11, 12, 13, 14, 14])
    A = Matrix([
        [   4,  -8,  4,   0],
        [  76, -68, -8,  24],
        [-122, 142,  4, -48],
        [  48, -72,  0,  24]])

    R = adjugate(M)
    assert isinstance(R, ImmutableDenseMatrix)
    assert R == A

    R = adjugate(SparseMatrix(M))
    assert isinstance(R, ImmutableSparseMatrix)
    assert R == A


def test_charpoly():
    M = Matrix([[1, 2], [3, 4]])
    A = PurePoly(x**2 - 5*x - 2, x, domain='RR')

    assert charpoly(M, x='x') == A
    assert charpoly(SparseMatrix(M), x='x') == A


def test_cofactor():
    M = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

    assert cofactor(M, 1, 1) == -12
    assert cofactor(SparseMatrix(M), 1, 1) == -12


def test_cofactor_matrix():
    M = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    A = Matrix([[-3, 6, -3], [6, -12, 6], [-3, 6, -3]])

    R = cofactor_matrix(M)
    assert isinstance(R, ImmutableDenseMatrix)
    assert R == A

    R = cofactor_matrix(SparseMatrix(M))
    assert isinstance(R, ImmutableSparseMatrix)
    assert R == A


def test_det():
    M = Matrix([[1, 2], [3, 4]])
    S = SparseMatrix(M)

    assert det(M) == -2
    assert det_bareiss(M) == -2
    assert det_berkowitz(M) == -2
    assert det_LU(M) == -2

    assert det(S) == -2
    assert det_bareiss(S) == -2
    assert det_berkowitz(S) == -2
    assert det_LU(S) == -2


def test_minor():
    M = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

    assert minor(M, 1, 1) == -12
    assert minor(SparseMatrix(M), 1, 1) == -12


def test_minor_submatrix():
    M = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    A = Matrix([[1, 3], [7, 9]])

    R = minor_submatrix(M, 1, 1)
    assert isinstance(R, ImmutableDenseMatrix)
    assert R == A

    R = minor_submatrix(SparseMatrix(M), 1, 1)
    assert isinstance(R, ImmutableSparseMatrix)
    assert R == A

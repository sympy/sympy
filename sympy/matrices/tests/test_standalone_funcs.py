from sympy import S, PurePoly, Symbol
from sympy.matrices import (Matrix, IMatrix, MMatrix,
    MutableSparseMatrix as MSMatrix, ImmutableSparseMatrix as ISMatrix)
from sympy.utilities.pytest import raises, XFAIL, skip, warns_deprecated_sympy

from sympy.matrices import (
    adjugate, charpoly, cofactor, cofactor_matrix,
    det, det_bareiss, det_berkowitz, det_LU, minor, minor_submatrix,
    is_echelon, echelon_form, rank, rref)

x = Symbol('x')


def test_adjugate():
    M = MMatrix(4, 4, [0, 1, 2, 3, 4, 5, 6, 7, 2, 9, 10, 11, 12, 13, 14, 14])
    A = Matrix([
        [   4,  -8,  4,   0],
        [  76, -68, -8,  24],
        [-122, 142,  4, -48],
        [  48, -72,  0,  24]])

    R = adjugate(M)
    assert isinstance(R, IMatrix)
    assert R == A

    R = adjugate(MSMatrix(M))
    assert isinstance(R, ISMatrix)
    assert R == A


def test_charpoly():
    M = MMatrix([[1, 2], [3, 4]])
    A = PurePoly(x**2 - 5*x - 2, x, domain='RR')

    assert charpoly(M, x='x') == A
    assert charpoly(MSMatrix(M), x='x') == A


def test_cofactor():
    M = MMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

    assert cofactor(M, 1, 1) == -12
    assert cofactor(MSMatrix(M), 1, 1) == -12


def test_cofactor_matrix():
    M = MMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    A = Matrix([[-3, 6, -3], [6, -12, 6], [-3, 6, -3]])

    R = cofactor_matrix(M)
    assert isinstance(R, IMatrix)
    assert R == A

    R = cofactor_matrix(MSMatrix(M))
    assert isinstance(R, ISMatrix)
    assert R == A


def test_det():
    M = MMatrix([[1, 2], [3, 4]])
    S = MSMatrix(M)

    assert det(M) == -2
    assert det_bareiss(M) == -2
    assert det_berkowitz(M) == -2
    assert det_LU(M) == -2

    assert det(S) == -2
    assert det_bareiss(S) == -2
    assert det_berkowitz(S) == -2
    assert det_LU(S) == -2


def test_minor():
    M = MMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

    assert cofactor(M, 1, 1) == -12
    assert cofactor(MSMatrix(M), 1, 1) == -12


def test_minor_submatrix():
    M = MMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    A = Matrix([[1, 3], [7, 9]])

    R = minor_submatrix(M, 1, 1)
    assert isinstance(R, IMatrix)
    assert R == A

    R = minor_submatrix(MSMatrix(M), 1, 1)
    assert isinstance(R, ISMatrix)
    assert R == A


def test_echelon():
    M = MMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    S = MSMatrix(M)
    A = Matrix([[1, 2, 3], [0, -3, -6], [0, 0, 0]])

    assert is_echelon(M) is False
    assert echelon_form(M) == A

    assert is_echelon(S) is False
    assert echelon_form(S) == A


def test_rank():
    M = MMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

    assert rank(M) == 2
    assert rank(MSMatrix(M)) == 2


def test_rref():
    M = MMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    A = (Matrix([[1, 0, -1], [0, 1, 2], [0, 0, 0]]), (0, 1))

    assert rref(M) == A
    assert rref(MSMatrix(M)) == A

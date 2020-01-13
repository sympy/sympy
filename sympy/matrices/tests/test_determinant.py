from sympy import PurePoly, symbols
from sympy.matrices import (Matrix, ImmutableDenseMatrix, SparseMatrix,
    ImmutableSparseMatrix)

from sympy.matrices.determinant import (
    adjugate, charpoly, cofactor, cofactor_matrix, det, det_bareiss,
    det_berkowitz, det_LU, minor, minor_submatrix)

x, y, z = symbols('x y z')


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


def test_determinant():

    for M in [Matrix(), Matrix([[1]])]:
        assert (
            M.det() ==
            M.det_bareiss() ==
            M.det_berkowitz() ==
            M.det_LU() ==
            1)

    M = Matrix(( (-3,  2),
                 ( 8, -5) ))

    assert M.det(method="bareiss") == -1
    assert M.det(method="berkowitz") == -1
    assert M.det(method="lu") == -1

    M = Matrix(( (x,   1),
                 (y, 2*y) ))

    assert M.det(method="bareiss") == 2*x*y - y
    assert M.det(method="berkowitz") == 2*x*y - y
    assert M.det(method="lu") == 2*x*y - y

    M = Matrix(( (1, 1, 1),
                 (1, 2, 3),
                 (1, 3, 6) ))

    assert M.det(method="bareiss") == 1
    assert M.det(method="berkowitz") == 1
    assert M.det(method="lu") == 1

    M = Matrix(( ( 3, -2,  0, 5),
                 (-2,  1, -2, 2),
                 ( 0, -2,  5, 0),
                 ( 5,  0,  3, 4) ))

    assert M.det(method="bareiss") == -289
    assert M.det(method="berkowitz") == -289
    assert M.det(method="lu") == -289

    M = Matrix(( ( 1,  2,  3,  4),
                 ( 5,  6,  7,  8),
                 ( 9, 10, 11, 12),
                 (13, 14, 15, 16) ))

    assert M.det(method="bareiss") == 0
    assert M.det(method="berkowitz") == 0
    assert M.det(method="lu") == 0

    M = Matrix(( (3, 2, 0, 0, 0),
                 (0, 3, 2, 0, 0),
                 (0, 0, 3, 2, 0),
                 (0, 0, 0, 3, 2),
                 (2, 0, 0, 0, 3) ))

    assert M.det(method="bareiss") == 275
    assert M.det(method="berkowitz") == 275
    assert M.det(method="lu") == 275

    M = Matrix(( (1, 0,  1,  2, 12),
                 (2, 0,  1,  1,  4),
                 (2, 1,  1, -1,  3),
                 (3, 2, -1,  1,  8),
                 (1, 1,  1,  0,  6) ))

    assert M.det(method="bareiss") == -55
    assert M.det(method="berkowitz") == -55
    assert M.det(method="lu") == -55

    M = Matrix(( (-5,  2,  3,  4,  5),
                 ( 1, -4,  3,  4,  5),
                 ( 1,  2, -3,  4,  5),
                 ( 1,  2,  3, -2,  5),
                 ( 1,  2,  3,  4, -1) ))

    assert M.det(method="bareiss") == 11664
    assert M.det(method="berkowitz") == 11664
    assert M.det(method="lu") == 11664

    M = Matrix(( ( 2,  7, -1, 3, 2),
                 ( 0,  0,  1, 0, 1),
                 (-2,  0,  7, 0, 2),
                 (-3, -2,  4, 5, 3),
                 ( 1,  0,  0, 0, 1) ))

    assert M.det(method="bareiss") == 123
    assert M.det(method="berkowitz") == 123
    assert M.det(method="lu") == 123

    M = Matrix(( (x, y, z),
                 (1, 0, 0),
                 (y, z, x) ))

    assert M.det(method="bareiss") == z**2 - x*y
    assert M.det(method="berkowitz") == z**2 - x*y
    assert M.det(method="lu") == z**2 - x*y

    # issue 13835
    a = symbols('a')
    M = lambda n: Matrix([[i + a*j for i in range(n)]
                          for j in range(n)])
    assert M(5).det() == 0
    assert M(6).det() == 0
    assert M(7).det() == 0


def test_legacy_det():
    # Minimal support for legacy keys for 'method' in det()
    # Partially copied from test_determinant()

    M = Matrix(( ( 3, -2,  0, 5),
                 (-2,  1, -2, 2),
                 ( 0, -2,  5, 0),
                 ( 5,  0,  3, 4) ))

    assert M.det(method="bareis") == -289
    assert M.det(method="det_lu") == -289
    assert M.det(method="det_LU") == -289

    M = Matrix(( (3, 2, 0, 0, 0),
                 (0, 3, 2, 0, 0),
                 (0, 0, 3, 2, 0),
                 (0, 0, 0, 3, 2),
                 (2, 0, 0, 0, 3) ))

    assert M.det(method="bareis") == 275
    assert M.det(method="det_lu") == 275
    assert M.det(method="Bareis") == 275

    M = Matrix(( (1, 0,  1,  2, 12),
                 (2, 0,  1,  1,  4),
                 (2, 1,  1, -1,  3),
                 (3, 2, -1,  1,  8),
                 (1, 1,  1,  0,  6) ))

    assert M.det(method="bareis") == -55
    assert M.det(method="det_lu") == -55
    assert M.det(method="BAREISS") == -55

    M = Matrix(( (-5,  2,  3,  4,  5),
                 ( 1, -4,  3,  4,  5),
                 ( 1,  2, -3,  4,  5),
                 ( 1,  2,  3, -2,  5),
                 ( 1,  2,  3,  4, -1) ))

    assert M.det(method="bareis") == 11664
    assert M.det(method="det_lu") == 11664
    assert M.det(method="BERKOWITZ") == 11664

    M = Matrix(( ( 2,  7, -1, 3, 2),
                 ( 0,  0,  1, 0, 1),
                 (-2,  0,  7, 0, 2),
                 (-3, -2,  4, 5, 3),
                 ( 1,  0,  0, 0, 1) ))

    assert M.det(method="bareis") == 123
    assert M.det(method="det_lu") == 123
    assert M.det(method="LU") == 123

from sympy import I, PurePoly, S, Symbol
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.matrices import (Matrix, ImmutableDenseMatrix as IMatrix,
    MutableDenseMatrix as MMatrix, MutableSparseMatrix as MSMatrix,
    ImmutableSparseMatrix as ISMatrix)

from sympy.matrices import (
    adjugate, charpoly, cofactor, cofactor_matrix,
    det_matrix, det_bareiss, det_berkowitz, det_LU, minor, minor_submatrix,
    rank,
    columnspace, nullspace, rowspace, orthogonalize,
    eigenvals, eigenvects, diagonalize, jordan_form, left_eigenvects,
    singular_values)

x = Symbol('x')


# determinant.py

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


def test_det_matrix():
    M = MMatrix([[1, 2], [3, 4]])
    S = MSMatrix(M)

    assert det_matrix(M) == -2
    assert det_bareiss(M) == -2
    assert det_berkowitz(M) == -2
    assert det_LU(M) == -2

    assert det_matrix(S) == -2
    assert det_bareiss(S) == -2
    assert det_berkowitz(S) == -2
    assert det_LU(S) == -2


def test_minor():
    M = MMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

    assert minor(M, 1, 1) == -12
    assert minor(MSMatrix(M), 1, 1) == -12


def test_minor_submatrix():
    M = MMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    A = Matrix([[1, 3], [7, 9]])

    R = minor_submatrix(M, 1, 1)
    assert isinstance(R, IMatrix)
    assert R == A

    R = minor_submatrix(MSMatrix(M), 1, 1)
    assert isinstance(R, ISMatrix)
    assert R == A


# reductions.py

def test_rank():
    M = MMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

    assert rank(M) == 2
    assert rank(MSMatrix(M)) == 2


# subspaces.py

def test_columnspace():
    M = MMatrix([
        [ 1,  2,  0,  2,  5],
        [-2, -5,  1, -1, -8],
        [ 0, -3,  3,  4,  1],
        [ 3,  6,  0, -7,  2]])

    A = columnspace(M)
    assert all(isinstance(m, IMatrix) for m in A)
    assert A[0] == Matrix([1, -2, 0, 3])
    assert A[1] == Matrix([2, -5, -3, 6])
    assert A[2] == Matrix([2, -1, 4, -7])

    A = columnspace(MSMatrix(M))
    assert all(isinstance(m, ISMatrix) for m in A)
    assert A[0] == Matrix([1, -2, 0, 3])
    assert A[1] == Matrix([2, -5, -3, 6])
    assert A[2] == Matrix([2, -1, 4, -7])


def test_rowspace():
    M = MMatrix([
        [ 1,  2,  0,  2,  5],
        [-2, -5,  1, -1, -8],
        [ 0, -3,  3,  4,  1],
        [ 3,  6,  0, -7,  2]])

    A = rowspace(M)
    assert all(isinstance(m, IMatrix) for m in A)
    assert A[0] == Matrix([[1, 2, 0, 2, 5]])
    assert A[1] == Matrix([[0, -1, 1, 3, 2]])
    assert A[2] == Matrix([[0, 0, 0, 5, 5]])

    A = rowspace(MSMatrix(M))
    assert all(isinstance(m, ISMatrix) for m in A)
    assert A[0] == Matrix([[1, 2, 0, 2, 5]])
    assert A[1] == Matrix([[0, -1, 1, 3, 2]])
    assert A[2] == Matrix([[0, 0, 0, 5, 5]])


def test_nullspace():
    M = MMatrix([
        [ 1,  3, 0,  2,  6, 3, 1],
        [-2, -6, 0, -2, -8, 3, 1],
        [ 3,  9, 0,  0,  6, 6, 2],
        [-1, -3, 0,  1,  0, 9, 3]])

    A = nullspace(M)
    assert all(isinstance(m, IMatrix) for m in A)
    assert A[0] == Matrix([-3, 1, 0, 0, 0, 0, 0])
    assert A[1] == Matrix([0, 0, 1, 0, 0, 0, 0])
    assert A[2] == Matrix([-2, 0, 0, -2, 1, 0, 0])
    assert A[3] == Matrix([0, 0, 0, 0, 0, -S(1)/3, 1])

    A = nullspace(MSMatrix(M))
    assert all(isinstance(m, ISMatrix) for m in A)
    assert A[0] == Matrix([-3, 1, 0, 0, 0, 0, 0])
    assert A[1] == Matrix([0, 0, 1, 0, 0, 0, 0])
    assert A[2] == Matrix([-2, 0, 0, -2, 1, 0, 0])
    assert A[3] == Matrix([0, 0, 0, 0, 0, -S(1)/3, 1])


def test_orthogonalize():
    v = [MMatrix([1, I]), MMatrix([1, -I])]

    A = orthogonalize(MMatrix, *v)
    assert all(isinstance(m, IMatrix) for m in A)
    assert A == [Matrix([1, I]), Matrix([1, -I])]

    A = orthogonalize(MSMatrix, *v)
    assert all(isinstance(m, ISMatrix) for m in A)
    assert A == [Matrix([1, I]), Matrix([1, -I])]


# eigen.py

def test_eigenvals():
    M = MMatrix([[0, 1, 1],
                [1, 0, 0],
                [1, 1, 1]])
    A = {2*S.One: 1, -S.One: 1, S.Zero: 1}

    R = eigenvals(M)
    assert R == A

    R = eigenvals(MSMatrix(M))
    assert R == A


def test_eigenvects():
    M = MMatrix([[0, 1, 1],
                [1, 0, 0],
                [1, 1, 1]])
    A = (
        [
            (-1, 1, [Matrix([-1, 1, 0])]),
            ( 0, 1, [Matrix([0, -1, 1])]),
            ( 2, 1, [Matrix([S(2)/3, S(1)/3, 1])])
        ])

    R = eigenvects(M)
    assert all(isinstance(ev, IMatrix) for v, m, es in R for ev in es)
    assert R == A

    R = eigenvects(MSMatrix(M))
    assert all(isinstance(ev, ISMatrix) for v, m, es in R for ev in es)
    assert R == A


def test_diagonalize():
    M = MMatrix(3, 3, [1, 2, 0, 0, 3, 0, 2, -4, 2])

    P, D = diagonalize(M)
    assert isinstance(P, IMatrix)
    assert isinstance(D, IMatrix)
    assert P.inv() * M * P == D

    P, D = diagonalize(MSMatrix(M))
    assert isinstance(P, ISMatrix)
    assert isinstance(D, ISMatrix)
    assert P.inv() * M * P == D


def test_jordan_form():
    M = MMatrix(3, 3, [0, 1, 0, -4, 4, 0, -2, 1, 2])
    N = MMatrix(3, 3, [2, 1, 0, 0, 2, 0, 0, 0, 2])

    P, J = jordan_form(M)
    assert isinstance(P, IMatrix)
    assert isinstance(J, IMatrix)
    assert J == N

    P, J = jordan_form(MSMatrix(M))
    assert isinstance(P, ISMatrix)
    assert isinstance(J, ISMatrix)
    assert J == N


def test_left_eigenvects():
    M = MMatrix([[0, 1, 1],
                [1, 0, 0],
                [1, 1, 1]])
    A = (
        [
            (-1, 1, [Matrix([[-2, 1, 1]])]),
            (0, 1, [Matrix([[-1, -1, 1]])]),
            (2, 1, [Matrix([[1, 1, 1]])])
        ])

    R = left_eigenvects(M)
    assert all(isinstance(ev, IMatrix) for v, m, es in R for ev in es)
    assert R == A

    R = left_eigenvects(MSMatrix(M))
    assert all(isinstance(ev, ISMatrix) for v, m, es in R for ev in es)
    assert R == A


def test_singular_values():
    M = MMatrix([
        [2, 4],
        [1, 3],
        [0, 0],
        [0, 0]
        ])
    A = [sqrt(sqrt(221) + 15), sqrt(15 - sqrt(221))]

    assert singular_values(M) == A
    assert singular_values(MSMatrix(M)) == A

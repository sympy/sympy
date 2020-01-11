from sympy import Float, I, PurePoly, S, Symbol
from sympy.matrices import (NonSquareMatrixError, Matrix,
    NumPyMatrix as NPMatrix)
from sympy.utilities.pytest import raises, XFAIL, skip, warns_deprecated_sympy

try:
    import numpy as np

    skipcheck = lambda: None

except ImportError:
    skipcheck = lambda: skip('NumPy must be available to test NumPyMatrix')

x = Symbol('x')


def test_creation():
    skipcheck()

    a = np.array([[1, 2, 3], [4, 5, 6]])
    M = NPMatrix(a)
    assert M == Matrix([[1, 2, 3], [4, 5, 6]])
    M = NPMatrix(3, 2, a)
    assert M == Matrix([[1, 2], [3, 4], [5, 6]])

    m = np.matrix([[1, 2, 3], [4, 5, 6]])
    M = NPMatrix(m)
    assert M == Matrix([[1, 2, 3], [4, 5, 6]])
    M = NPMatrix(3, 2, m)
    assert M == Matrix([[1, 2], [3, 4], [5, 6]])

    a = np.array([[1, 2], [3, 4]], dtype='i4')
    M = NPMatrix(m)
    assert M.dtype == np.int32
    M = NPMatrix(m, dtype='f4')
    assert M.dtype == np.float32


def test_as_wrapper():
    skipcheck()

    a = np.array([[1, 2], [3, 4]])
    M = NPMatrix(a, copy=False)
    assert M[0] == 1
    a[0, 0] = -1
    assert M[0] == -1

    m = np.matrix([[1, 2], [3, 4]])
    M = NPMatrix(m, copy=False)
    assert M[0] == 1
    m[0, 0] = -1
    assert M[0] == -1


def test_as_reshaping_wrapper():
    skipcheck()

    a = np.array([[1, 2, 3], [4, 5, 6]])
    M = NPMatrix(3, 2, a, copy=False)
    assert M == Matrix([[1, 2], [3, 4], [5, 6]])
    a[0, 0] = -1
    assert M[0] == -1

    m = np.matrix([[1, 2, 3], [4, 5, 6]])
    M = NPMatrix(3, 2, m, copy=False)
    assert M == Matrix([[1, 2], [3, 4], [5, 6]])
    m[0, 0] = -1
    assert M[0] == -1


def test_adjugate():
    skipcheck()

    M = NPMatrix(4, 4, [0, 1, 2, 3, 4, 5, 6, 7, 2, 9, 10, 11, 12, 13, 14, 14])
    A = M.adjugate()

    assert isinstance(A, NPMatrix)
    assert all(Float(a).epsilon_eq(b, epsilon='1e-12') for a, b in zip(A,
            (4, -8, 4, 0, 76, -68, -8, 24, -122, 142, 4, -48, 48, -72, 0, 24)))

    M = NPMatrix(8, 2, M)
    raises(NonSquareMatrixError, lambda: M.adjugate())


def test_charpoly():
    skipcheck()

    M = NPMatrix([[1, 2], [3, 4]])
    assert M.charpoly(x='x') == PurePoly(x**2 - 5*x - 2, x, domain='RR')

    M = NPMatrix(4, 1, M)
    raises(NonSquareMatrixError, lambda: M.charpoly())


def test_cofactor():
    skipcheck()

    M = NPMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    assert Float(M.cofactor(1, 1)).epsilon_eq(-12)


def test_cofactor_matrix():
    skipcheck()

    M = NPMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    A = M.cofactor_matrix()

    assert isinstance(A, NPMatrix)
    assert all(Float(a).epsilon_eq(b, epsilon='1e-12') for a, b in zip(A,
            (-3, 6, -3, 6, -12, 6, -3, 6, -3)))


def test_det():
    skipcheck()

    M = NPMatrix([[1, 2], [3, 4]])
    assert Float(M.det()).epsilon_eq(-2)
    assert Float(M.det_bareiss()).epsilon_eq(-2)
    assert Float(M.det_berkowitz()).epsilon_eq(-2)
    assert Float(M.det_LU()).epsilon_eq(-2)

    M = NPMatrix(4, 1, M)
    raises(NonSquareMatrixError, lambda: M.det())
    raises(NonSquareMatrixError, lambda: M.det_bareiss())
    raises(NonSquareMatrixError, lambda: M.det_berkowitz())
    raises(NonSquareMatrixError, lambda: M.det_LU())


def test_minor():
    skipcheck()

    M = NPMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    assert Float(M.cofactor(1, 1)).epsilon_eq(-12)


def test_minor_submatrix():
    skipcheck()

    M = NPMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    A = M.minor_submatrix(1, 1)

    assert isinstance(A, NPMatrix)
    assert A == Matrix([[1, 3], [7, 9]])


def test_echelon():
    skipcheck()

    M = NPMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype='i4') # float gives different result
    A = M.echelon_form()

    assert M.is_echelon is False
    assert isinstance(A, NPMatrix)
    assert A == Matrix([[1, 2, 3], [0, -3, -6], [0, 0, 0]])


def test_rank():
    skipcheck()

    M = NPMatrix([[S('1+I'),S('-19/4+5/4*I'),S('1/2-I'),S('9/4+55/16*I'),S('-3/4'),S('45/32-37/16*I'),S('1/4+1/2*I'),S('-129/64-9/64*I'),S('1/4-5/16*I'),S('65/128+87/64*I'),S('-9/32-1/16*I'),S('183/256-97/128*I'),S('3/64+13/64*I'),S('-23/32-59/256*I'),S('15/128-3/32*I'),S('19/256+551/1024*I')],
        [S('21/8+I'),S('-537/64+143/16*I'),S('-5/8-39/16*I'),S('2473/256+137/64*I'),S('-149/64+49/32*I'),S('-177/128-1369/128*I'),S('125/64+87/64*I'),S('-2063/256+541/128*I'),S('85/256-33/16*I'),S('805/128+2415/512*I'),S('-219/128+115/256*I'),S('6301/4096-6609/1024*I'),S('119/128+143/128*I'),S('-10879/2048+4343/4096*I'),S('129/256-549/512*I'),S('42533/16384+29103/8192*I')],
        [S('-2'),S('17/4-13/2*I'),S('1+I'),S('-19/4+5/4*I'),S('1/2-I'),S('9/4+55/16*I'),S('-3/4'),S('45/32-37/16*I'),S('1/4+1/2*I'),S('-129/64-9/64*I'),S('1/4-5/16*I'),S('65/128+87/64*I'),S('-9/32-1/16*I'),S('183/256-97/128*I'),S('3/64+13/64*I'),S('-23/32-59/256*I')],
        [S('1/4+13/4*I'),S('-825/64-147/32*I'),S('21/8+I'),S('-537/64+143/16*I'),S('-5/8-39/16*I'),S('2473/256+137/64*I'),S('-149/64+49/32*I'),S('-177/128-1369/128*I'),S('125/64+87/64*I'),S('-2063/256+541/128*I'),S('85/256-33/16*I'),S('805/128+2415/512*I'),S('-219/128+115/256*I'),S('6301/4096-6609/1024*I'),S('119/128+143/128*I'),S('-10879/2048+4343/4096*I')],
        [S('-4*I'),S('27/2+6*I'),S('-2'),S('17/4-13/2*I'),S('1+I'),S('-19/4+5/4*I'),S('1/2-I'),S('9/4+55/16*I'),S('-3/4'),S('45/32-37/16*I'),S('1/4+1/2*I'),S('-129/64-9/64*I'),S('1/4-5/16*I'),S('65/128+87/64*I'),S('-9/32-1/16*I'),S('183/256-97/128*I')],
        [S('1/4+5/2*I'),S('-23/8-57/16*I'),S('1/4+13/4*I'),S('-825/64-147/32*I'),S('21/8+I'),S('-537/64+143/16*I'),S('-5/8-39/16*I'),S('2473/256+137/64*I'),S('-149/64+49/32*I'),S('-177/128-1369/128*I'),S('125/64+87/64*I'),S('-2063/256+541/128*I'),S('85/256-33/16*I'),S('805/128+2415/512*I'),S('-219/128+115/256*I'),S('6301/4096-6609/1024*I')],
        [S('-4'),S('9-5*I'),S('-4*I'),S('27/2+6*I'),S('-2'),S('17/4-13/2*I'),S('1+I'),S('-19/4+5/4*I'),S('1/2-I'),S('9/4+55/16*I'),S('-3/4'),S('45/32-37/16*I'),S('1/4+1/2*I'),S('-129/64-9/64*I'),S('1/4-5/16*I'),S('65/128+87/64*I')],
        [S('-2*I'),S('119/8+29/4*I'),S('1/4+5/2*I'),S('-23/8-57/16*I'),S('1/4+13/4*I'),S('-825/64-147/32*I'),S('21/8+I'),S('-537/64+143/16*I'),S('-5/8-39/16*I'),S('2473/256+137/64*I'),S('-149/64+49/32*I'),S('-177/128-1369/128*I'),S('125/64+87/64*I'),S('-2063/256+541/128*I'),S('85/256-33/16*I'),S('805/128+2415/512*I')],
        [S('0'),S('-6'),S('-4'),S('9-5*I'),S('-4*I'),S('27/2+6*I'),S('-2'),S('17/4-13/2*I'),S('1+I'),S('-19/4+5/4*I'),S('1/2-I'),S('9/4+55/16*I'),S('-3/4'),S('45/32-37/16*I'),S('1/4+1/2*I'),S('-129/64-9/64*I')],
        [S('1'),S('-9/4+3*I'),S('-2*I'),S('119/8+29/4*I'),S('1/4+5/2*I'),S('-23/8-57/16*I'),S('1/4+13/4*I'),S('-825/64-147/32*I'),S('21/8+I'),S('-537/64+143/16*I'),S('-5/8-39/16*I'),S('2473/256+137/64*I'),S('-149/64+49/32*I'),S('-177/128-1369/128*I'),S('125/64+87/64*I'),S('-2063/256+541/128*I')],
        [S('0'),S('-4*I'),S('0'),S('-6'),S('-4'),S('9-5*I'),S('-4*I'),S('27/2+6*I'),S('-2'),S('17/4-13/2*I'),S('1+I'),S('-19/4+5/4*I'),S('1/2-I'),S('9/4+55/16*I'),S('-3/4'),S('45/32-37/16*I')],
        [S('0'),S('1/4+1/2*I'),S('1'),S('-9/4+3*I'),S('-2*I'),S('119/8+29/4*I'),S('1/4+5/2*I'),S('-23/8-57/16*I'),S('1/4+13/4*I'),S('-825/64-147/32*I'),S('21/8+I'),S('-537/64+143/16*I'),S('-5/8-39/16*I'),S('2473/256+137/64*I'),S('-149/64+49/32*I'),S('-177/128-1369/128*I')]])
    assert M.rank() == 8


def test_rref():
    skipcheck()

    M = NPMatrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    A = M.rref()

    assert isinstance(A[0], NPMatrix)
    assert A == (Matrix([[1, 0, -1], [0, 1, 2], [0, 0, 0]]), (0, 1))


def test_columnspace():
    skipcheck()

    M = NPMatrix([
        [ 1,  2,  0,  2,  5],
        [-2, -5,  1, -1, -8],
        [ 0, -3,  3,  4,  1],
        [ 3,  6,  0, -7,  2]])
    A = M.columnspace()

    assert all(isinstance(m, NPMatrix) for m in A)
    assert A[0] == Matrix([1, -2, 0, 3])
    assert A[1] == Matrix([2, -5, -3, 6])
    assert A[2] == Matrix([2, -1, 4, -7])


def test_rowspace():
    skipcheck()

    M = NPMatrix([
        [ 1,  2,  0,  2,  5],
        [-2, -5,  1, -1, -8],
        [ 0, -3,  3,  4,  1],
        [ 3,  6,  0, -7,  2]], dtype='i4') # float gives different result
    A = M.rowspace()

    assert all(isinstance(m, NPMatrix) for m in A)
    assert A[0] == Matrix([[1, 2, 0, 2, 5]])
    assert A[1] == Matrix([[0, -1, 1, 3, 2]])
    assert A[2] == Matrix([[0, 0, 0, 5, 5]])


def test_nullspace():
    skipcheck()

    M = NPMatrix([
        [ 1,  3, 0,  2,  6, 3, 1],
        [-2, -6, 0, -2, -8, 3, 1],
        [ 3,  9, 0,  0,  6, 6, 2],
        [-1, -3, 0,  1,  0, 9, 3]])
    A = M.nullspace()

    assert all(isinstance(m, NPMatrix) for m in A)
    assert A[0] == Matrix([-3, 1, 0, 0, 0, 0, 0])
    assert A[1] == Matrix([0, 0, 1, 0, 0, 0, 0])
    assert A[2] == Matrix([-2, 0, 0, -2, 1, 0, 0])
    assert A[3] == Matrix([0, 0, 0, 0, 0, -1.0/3, 1])


def test_orthogonalize():
    skipcheck()

    v = [Matrix([1, I]), Matrix([1, -I])]
    A = NPMatrix.orthogonalize(*v)

    assert all(isinstance(m, NPMatrix) for m in A)
    assert Float(A[0][0]).epsilon_eq(1)
    assert Float(A[0][1] / I).epsilon_eq(1)
    assert Float(A[1][0]).epsilon_eq(1)
    assert Float(A[1][1] / I).epsilon_eq(-1)
